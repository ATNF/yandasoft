/// @copyright (c) 2007 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///

#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/ImageParamsHelper.h>
#include <imageaccess/ImageAccessFactory.h>
#include <utils/PolConverter.h>
#include <fitting/Axes.h>
#include <utils/MultiDimArrayPlaneIter.h>
#include <utils/PaddingUtils.h>
#include <gridding/SupportSearcher.h>
#include <casacore/lattices/LatticeMath/Fit2D.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.synthesisparamshelper");

#include <askap/AskapError.h>
#include <profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Quanta.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/LinearCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>

#include <Common/ParameterSet.h>
#include <Common/Exceptions.h>

#include <casacore/measures/Measures/Stokes.h>

#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <string>

using namespace askap::scimath;
using namespace casa;

namespace askap
{
  namespace synthesis
  {
    /// @brief image accessor
    boost::shared_ptr<accessors::IImageAccess> SynthesisParamsHelper::theirImageAccessor;

    /// @brief default frequency frame
    casa::MFrequency::Ref SynthesisParamsHelper::theirFreqFrame(casa::MFrequency::TOPO);

    void SynthesisParamsHelper::setUpImages(const askap::scimath::Params::ShPtr& params,
				const LOFAR::ParameterSet &parset)
    {
      ASKAPDEBUGTRACE("SynthesisParamsHelper::setUpImages");
      try {
	     vector<string> images=parset.getStringVector("Names");
         std::vector<int> shape=parset.getInt32Vector("shape",std::vector<int>());
         std::vector<std::string> cellsize=parset.getStringVector("cellsize",std::vector<std::string>());

         for (vector<string>::const_iterator it=images.begin();it!=images.end();++it) {
              ASKAPCHECK(it->find("image") == 0, "All image names given in Names are supposed to start from 'image', you have "<<
                         *it);
              int nchan=parset.getInt32(*it+".nchan");
              std::vector<double> freq=parset.getDoubleVector(*it+".frequency");
              ASKAPCHECK(freq.size()>=2, "Parameter "<<*it<<".frequency should have at least 2-elements, you have "<<freq.size());
              if (parset.isDefined(*it+".shape")) {
                  if (shape.size()!=0) {
                      ASKAPLOG_INFO_STR(logger, "Global image shape "<<shape<<
                           " is overridden by specialisation for the parameter "<<*it<<", new shape is "<<
                           parset.getInt32Vector(*it+".shape"));
                  }
                  shape=parset.getInt32Vector(*it+".shape");
              } else {
                  ASKAPCHECK(shape.size()!=0, "If a global image shape is not defined, you should give shape separately for every image");
              }
              if (parset.isDefined(*it+".cellsize")) {
                  if (cellsize.size()!=0) {
                      ASKAPLOG_INFO_STR(logger, "Global cell size "<<cellsize<<
                           " is overridden by specialisation for the parameter "<<*it<<", new cell size is "<<
                            parset.getStringVector(*it+".cellsize"));
                  }
                  cellsize=parset.getStringVector(*it+".cellsize");
              } else {
                  ASKAPCHECK(cellsize.size()!=0, "If a global cell size is not defined, you should give it separately for every image");
              }

              const int nfacets = parset.getInt32(*it+".nfacets",1);
              ASKAPCHECK(nfacets>0, "Number of facets is supposed to be a positive number, you gave "<<nfacets);
              ASKAPCHECK(shape.size()>=2, "Image is supposed to be at least two dimensional. "<<
                          "check shape parameter, you gave "<<shape);

              const std::vector<std::string> direction=parset.getStringVector(*it+".direction");
              const std::vector<std::string> tangent = parset.getStringVector(*it+".tangent",std::vector<std::string>());
              if (tangent.size() != 0) {
                  ASKAPCHECK(nfacets == 1, "Faceting with user-defined tangent point is not supported");
              }

              // required polarisation
              if (!parset.isDefined(*it+".polarisation")) {
                  ASKAPLOG_INFO_STR(logger, "Polarisation frame is not defined for "<<*it
                                    <<", only stokes I will be generated");
              }
              const std::vector<std::string> stokesVec = parset.getStringVector(*it+".polarisation",
                             std::vector<std::string>(1,"I"));
              // there could be many ways to define stokes, e.g. ["XX YY"] or ["XX","YY"] or "XX,YY"
              // to allow some flexibility we have to concatenate all elements first and then
              // allow the parser from PolConverter to take care of extracting the products.
              std::string stokesStr;
              for (size_t i=0; i<stokesVec.size(); ++i) {
                   stokesStr += stokesVec[i];
              }
              casa::Vector<casa::Stokes::StokesTypes> stokes = scimath::PolConverter::fromString(stokesStr);

              const bool ewProj = parset.getBool(*it+".ewprojection", false);
              if (ewProj) {
                  ASKAPLOG_INFO_STR(logger, "Image parameter "<< *it<<" will have SCP/NCP projection");
              } else {
                  ASKAPLOG_INFO_STR(logger, "Image parameter "<< *it<<" will have plain SIN projection");
              }

              const int nTaylorTerms = parset.getInt32(*it+".nterms",1);
              ASKAPCHECK(nTaylorTerms>0, "Number of Taylor terms is supposed to be a positive number, you gave "<<
                         nTaylorTerms);
              //
              ImageParamsHelper iph(*it);
              for (int order = 0; order<2*nTaylorTerms-1; ++order) {
                   if (nTaylorTerms>1) {
                       // this is an MFS case, setup Taylor terms
                       iph.makeTaylorTerm(order);
                       ASKAPLOG_INFO_STR(logger,"Setting up Taylor term "<<order);
                   }
                   if (nfacets == 1) {
                       ASKAPLOG_INFO_STR(logger, "Setting up new empty image "<< iph.paramName());
                       if (tangent.size()) {
                           add(*params, iph.paramName(), tangent, cellsize, shape, ewProj, freq[0], freq[1], nchan,stokes, direction);
                       } else {
		                   add(*params, iph.paramName(), direction, cellsize, shape, ewProj, freq[0], freq[1], nchan,stokes);
		               }
		           } else {
		               // this is a multi-facet case
		               ASKAPLOG_INFO_STR(logger, "Setting up "<<nfacets<<" x "<<nfacets<<
		                                         " new empty facets for image "<< iph.paramName());
		               const int facetstep = parset.getInt32(*it+".facetstep",casa::min(shape[0],shape[1]));
		               ASKAPCHECK(facetstep>0, "facetstep parameter is supposed to be positive, you have "<<facetstep);
		               ASKAPLOG_INFO_STR(logger, "Facet centers will be "<<facetstep<<
		                           " pixels apart, each facet size will be "<<shape[0]<<" x "<<shape[1]);
		               add(*params, iph.paramName(), direction, cellsize, shape, ewProj, freq[0], freq[1], nchan, stokes, nfacets,facetstep);
		           }
		      }
		      ASKAPLOG_INFO_STR(logger, "Number of channels = "<<nchan);
		      ASKAPLOG_INFO_STR(logger, "Polarisation planes correspond to "<<scimath::PolConverter::toString(stokes));
	     }
	  }
	  catch (const LOFAR::APSException &ex) {
	      throw AskapError(ex.what());
	  }
	}

    /// @brief load images according to the parset file
	/// @details This method is somewhat analogous to setUpImages, but it loads the images
	/// from the disk instead of setting them up from the scratch. Encapsulation of all loading
	/// of multiple images in a single method is required to provide a seamless handling of
	/// the faceted image.
	/// @param[in] params Images to be created here
	/// @param[in] parset a parset object to read the parameters from
	void SynthesisParamsHelper::loadImages(const askap::scimath::Params::ShPtr& params, const LOFAR::ParameterSet &parset)
    {
      ASKAPDEBUGTRACE("SynthesisParamsHelper::loadImages");
      ASKAPDEBUGASSERT(params);
      try {
         const vector<string> images=parset.getStringVector("Names");
         for (vector<string>::const_iterator ci = images.begin(); ci != images.end(); ++ci) {
              ASKAPCHECK(ci->find("image") == 0, "All image names given in Names are supposed to start from 'image', you have "<<
                         *ci);
              // @todo add more checking that the image loaded from the disk conforms with the
              // parameters given in the parset file

              const int nfacets = parset.getInt32(*ci+".nfacets",1);
              ASKAPCHECK(nfacets>0, "Number of facets is supposed to be a positive number, you gave "<<nfacets);

              const int nTaylorTerms = parset.getInt32(*ci+".nterms",1);
              ASKAPCHECK(nTaylorTerms>0, "Number of Taylor terms is supposed to be a positive number, you gave "<<
                         nTaylorTerms);
              //
              ImageParamsHelper iph(*ci);
              for (int order = 0; order<2*nTaylorTerms-1; ++order) {
                   if (nTaylorTerms>1) {
                       // this is an MFS case, setup Taylor terms
                       iph.makeTaylorTerm(order);
                       ASKAPLOG_INFO_STR(logger,"Processing Taylor term "<<order);
                   }
                   if (nfacets == 1) {
                       ASKAPLOG_INFO_STR(logger, "Reading image "<<iph.paramName());
                       SynthesisParamsHelper::loadImageParameter(*params,iph.paramName(),iph.paramName());
                   } else {
                       ASKAPLOG_INFO_STR(logger, "Loading multi-facet image image "<<iph.paramName());
                       SynthesisParamsHelper::getMultiFacetImage(*params,iph.paramName(),iph.paramName(), nfacets);
                   }
              }
         }

	  }
	  catch (const LOFAR::APSException &ex) {
	      throw AskapError(ex.what());
	  }
    }

    /// @brief helper method to get projection
    /// @details We support both standard SIN projection and SCP/NCP variants (for East-West arrays).
    /// This method encapsulates the logic and returns a projection class
    /// @param[in] ewprojection true for SCP/NCP variant, false otherwise
    /// @param[in] dec declination in radians (unused for standard SIN projection)
    /// @return casa::Projection class
    casa::Projection SynthesisParamsHelper::getProjection(const bool ewprojection, const double dec)
    {
       if (ewprojection) {
           const double sdec = sin(dec);
           ASKAPCHECK(sdec != 0., "Singular SCP/NCP projection dec="<<dec/casa::C::pi*180.<<
                      " deg, sin(dec)="<<sdec<<". Use plain SIN projection instead!");
           casa::Vector<casa::Double> projectionParameters(2,0.);
           projectionParameters(1) = cos(dec) / sdec;
           return casa::Projection(casa::Projection::SIN, projectionParameters);
       }
       return casa::Projection(casa::Projection::SIN);
    }


    void SynthesisParamsHelper::add(askap::scimath::Params& ip,
				    const string& name, const vector<string>& direction,
				    const vector<string>& cellsize, const vector<int>& shape,
				    const bool ewprojection,
				    const double freqmin, const double freqmax, const int nchan,
				    const casa::Vector<casa::Stokes::StokesTypes> &stokes,
				    const vector<string>& centreDir)
    {
      int nx=shape[0];
      int ny=shape[1];
      ASKAPCHECK(cellsize.size() == 2, "Cell size should have exactly 2 parameters, you have "<<cellsize.size());
      ASKAPCHECK(direction.size() == 3, "Direction should have exactly 3 parameters, you have "<<direction.size());
      ASKAPCHECK(direction[2] == "J2000", "Only J2000 is implemented at the moment, you have requested "<<direction[2]);
      ASKAPCHECK((centreDir.size() == 0) || (centreDir.size() == 3),
                 "centreDir should have exactly 3 parameters or be empty, you have "<<centreDir.size());
      ASKAPCHECK(stokes.nelements()>=1, "At least one polarisation plane should be defined, you have defined none");

      const double ra = convertQuantity(direction[0],"rad");
      const double dec = convertQuantity(direction[1],"rad");

      const double xcellsize =-1.0*convertQuantity(cellsize[0],"rad");
      const double ycellsize = convertQuantity(cellsize[1],"rad");

      /// @todo Do something with the frame info in direction[2]
      Axes axes;
      casa::Matrix<double> xform(2,2,0.);
      xform.diagonal() = 1.;
      // direction coordinate corresponding to the case with tangent point in the image centre
      const casa::DirectionCoordinate dcTangent(casa::MDirection::J2000,
                  getProjection(ewprojection, dec), ra,dec,xcellsize,ycellsize,xform,nx/2,ny/2);
      if (centreDir.size()) {
          ASKAPCHECK(centreDir[2] == "J2000", "Only J2000 is implemented at the moment, you have requested "<<centreDir[2]);
          ASKAPLOG_INFO_STR(logger, "Image parameter "<<name<<" have tangent point "<<direction<<" and image centre "<<centreDir);
          const double raCentre = convertQuantity(centreDir[0],"rad");
          const double decCentre = convertQuantity(centreDir[1],"rad");
          const casa::MVDirection centre(raCentre, decCentre);
          casa::Vector<casa::Double> pix;
          dcTangent.toPixel(pix,centre);
          ASKAPDEBUGASSERT(pix.nelements() == 2);
          axes.addDirectionAxis(casa::DirectionCoordinate(casa::MDirection::J2000,
                  getProjection(ewprojection, dec), ra,dec,xcellsize,ycellsize,xform,double(nx)-pix[0],double(ny)-pix[1]));
      } else {
        ASKAPLOG_INFO_STR(logger, "Image parameter "<<name<<" have tangent point "<<direction<<" at the image centre");
        axes.addDirectionAxis(dcTangent);
      }
      axes.addStokesAxis(stokes);

      casa::Array<double> pixels(casa::IPosition(4, nx, ny, stokes.nelements(), nchan));
      pixels.set(0.0);
      axes.add("FREQUENCY", freqmin, freqmax);
      ASKAPLOG_INFO_STR(logger, "Spectral axis will have startFreq="<<freqmin<<" Hz, endFreq="<<freqmax<<
                                "Hz, nChan="<<nchan);
      ip.add(name, pixels, axes);
    }

    /// @brief Add a parameter as a faceted image
    /// @param[in] ip Parameters
    /// @param[in] name Name of parameter
    /// @param[in] direction Strings containing [ra, dec, frame] (common tangent point)
    /// @param[in] cellsize Cellsize as a string e.g. [12arcsec, 12arcsec]
    /// @param[in] shape Number of pixels in RA and DEC for each facet e.g. [256, 256]
    /// @param[in] ewprojection If true, SCP or NCP variant of SIN projection will be used
    /// @param[in] freqmin Minimum frequency (Hz)
    /// @param[in] freqmax Maximum frequency (Hz)
    /// @param[in] nchan Number of spectral channels
    /// @param[in] stokes Polarisation frame (vector of stokes enums)
    /// @param[in] nfacets Number of facets in each axis (assumed the same for both axes)
    /// @param[in] facetstep Offset in pixels between facet centres (equal to shape to
    ///            have no overlap between adjacent facets), assumed the same for both axes
    void SynthesisParamsHelper::add(askap::scimath::Params& ip, const string& name,
       const vector<string>& direction,
       const vector<string>& cellsize,
       const vector<int>& shape,
       const bool ewprojection,
       const double freqmin, const double freqmax, const int nchan,
       const casa::Vector<casa::Stokes::StokesTypes> &stokes,
       const int nfacets, const int facetstep)
    {
      ASKAPDEBUGASSERT(nfacets>0);
      ASKAPDEBUGASSERT(facetstep>0);
      const int nx=shape[0];
      const int ny=shape[1];
      ASKAPCHECK(cellsize.size() == 2, "Cell size should have exactly 2 parameters, you have "<<cellsize.size());
      ASKAPCHECK(direction.size() == 3, "Direction should have exactly 3 parameters, you have "<<direction.size());
      ASKAPCHECK(direction[2] == "J2000", "Only J2000 is implemented at the moment, you have requested "<<direction[2]);
      ASKAPCHECK(stokes.nelements()>=1, "At least one polarisation plane should be defined, you have defined none");

      const double ra = convertQuantity(direction[0],"rad");
      const double dec = convertQuantity(direction[1],"rad");

      const double xcellsize =-1.0*convertQuantity(cellsize[0],"rad");
      const double ycellsize = convertQuantity(cellsize[1],"rad");

      // zero-filled array is the same for all facets as it is copied inside Params
      // class
      casa::Array<double> pixels(casa::IPosition(4, nx, ny, stokes.nelements(), nchan));
      pixels.set(0.0);

      // void linear transform used to set up coordinate system
      casa::Matrix<double> xform(2,2,0.);
      xform.diagonal() = 1.;

      // have to create facet parameter in two steps as it could be
      // a Taylor decomposition
      ImageParamsHelper iph(name);
      // a loop over facets
      const double facetFactor = (nfacets-1)/2.;
      for (int ix=0;ix<nfacets;++ix) {
           // reference position (tangent point) on the first axis assuming the pixels of this
           // facet go from 0 to nx-1.
           const double xrefpix = double(nx)/2.-facetstep*(ix-facetFactor);
           for (int iy=0;iy<nfacets;++iy) {

                // for debugging to avoid running out of memory
                //if (ix!=1 || iy!=1) continue;

                // reference position (tangent point) on the second axis assuming the pixels of this
                // facet go from 0 to ny-1.
                const double yrefpix = double(ny)/2.-facetstep*(iy-facetFactor);

                /// @todo Do something with the frame info in direction[2]
                Axes axes;
                axes.addDirectionAxis(casa::DirectionCoordinate(casa::MDirection::J2000,
                  getProjection(ewprojection, dec), ra,dec,xcellsize,ycellsize,xform,xrefpix,yrefpix));

                // a fake axis to know which part of the image actually contains useful
                // information. Otherwise, this parameter is impossible to derive from a
                // single facet only (and we may need, e.g., to clip the outer edges in each
                // major cycle)
                axes.add("FACETSTEP",double(facetstep),double(facetstep));

                axes.addStokesAxis(stokes);

                axes.add("FREQUENCY", freqmin, freqmax);

                // add/change facet indices
                iph.makeFacet(ix,iy);
                ip.add(iph.paramName(), pixels, axes);

                // for debigging
                //if (ix!=0 || iy!=0) ip.fix(iph.paramName());
           }
      }

    }

    /// @brief helper method to clip the outer edges of the image
    /// @details For experiments with faceting we want to be able to clip the outer
    /// edges of each model image (beyond the facet step) to zero. This is one way to
    /// reduce cross-talk problem (when facets overlap). This method encapsulates all
    /// the required operations. It takes facet step from the fake image axis FACETSTEP
    /// and does nothing if such a parameter doesn't exist or is larger than the shape
    /// along the directional axes.
    /// @param[in] ip parameters
    /// @param[in] name full name of the image (i.e. with .facet.x.y for facets)
    void SynthesisParamsHelper::clipImage(const askap::scimath::Params &ip, const string &name)
    {
       const askap::scimath::Axes axes(ip.axes(name));
       if (!axes.has("FACETSTEP")) {
           // it is not a facet image, do nothing.
           return;
       }
       const int facetStep = int(axes.start("FACETSTEP"));
       ASKAPDEBUGASSERT(facetStep>0);
       casa::Array<double> pixels = ip.value(name);
       const casa::IPosition shape = pixels.shape();
       ASKAPDEBUGASSERT(shape.nelements()>=2);
       casa::IPosition end(shape);
       for (uint index=0;index<end.nelements();++index) {
            ASKAPDEBUGASSERT(end[index]>=1);
            end[index]--;
       }

       if (shape[0]>facetStep+1) {
           // need clipping along the first axis
           casa::IPosition start(shape.nelements(),0);
           end[0] = (shape[0]-facetStep)/2-1;
           end[1] = shape[1]-1; // although this step is strictly speaking unnecessary
           pixels(start,end).set(0.);

           end[0] = shape[0]-1;
           start[0] = (shape[0]+facetStep)/2;
           pixels(start,end).set(0.);
       }

       if (shape[1]>facetStep+1) {
           // need clipping along the second axis
           casa::IPosition start(shape.nelements(),0);
           start[0]=(shape[0]-facetStep)/2;
           end[0]=(shape[0]+facetStep)/2;
           if (start[0]<0) {
               start[0] = 0;
           }
           if (end[0]+1 > shape[0]) {
               end[0] = shape[0] - 1;
           }
           start[1] = 0;
           end[1] = (shape[1]-facetStep)/2-1;
           pixels(start,end).set(0.);

           start[1] = (shape[1]+facetStep)/2;
           end[1] = shape[1]-1;
           pixels(start,end).set(0.);
       }
    }

    /// @brief helper method to store restoring beam for an image
    /// @details We have to carry restore beam parameters together with the image.
    /// This is done by creating 2 fake axes MAJMIN (with start = maj and end = min)
    /// and PA with position angle. All angles are given in radians. The presence of
    /// this fake axes distinguishes a restored image from model image. Restored image
    /// will have units Jy/beam instead of Jy/pixel and beam info will be added to the
    /// image (in saveImageParamter).
    /// @param[in] ip parameters
    /// @param[in] name full name of the parameter representing this image
    /// @param[in] beam major, minor axes and position anlge as quantities
    void SynthesisParamsHelper::setBeam(askap::scimath::Params &ip, const string &name,
                            const casa::Vector<casa::Quantum<double> > &beam)
    {
       askap::scimath::Axes &axes = ip.axes(name);
       ASKAPDEBUGASSERT(beam.nelements()>=3);
       if (axes.has("MAJMIN")) {
           axes.update("MAJMIN",beam[0].getValue("rad"),beam[1].getValue("rad"));
       } else {
           axes.add("MAJMIN",beam[0].getValue("rad"),beam[1].getValue("rad"));
       }

       if (axes.has("PA")) {
           axes.update("PA",beam[2].getValue("rad"),0.);
       } else {
           axes.add("PA",beam[2].getValue("rad"),0.);
       }
    }

    /// @brief add a parameter as a merged faceted image
    /// @details Each facet is represented by a number of independent parameters with
    /// the appropriate names. This method looks at the coordinate systems of all
    /// subimages and forms a parameter representing merged image. It can then be
    /// populated with the data from the appropriate slices.
    /// Note that this function is called for the first facet only, thus
    /// the facet center is that of facet 0.0.
    /// @param[in] ip parameters
    /// @param[in] name Base name of the parameter (i.e. without .facet.0.0)
    /// @param[in] nfacets number of facets defined
    void SynthesisParamsHelper::add(askap::scimath::Params& ip, const string &name,
              const int nfacets)
    {
       ASKAPDEBUGASSERT(nfacets>1);
       // no consistency check of the coordinate systems of individual patches at this stage

       // create image handler in two steps because name may contain a taylor-order suffix
       ImageParamsHelper iph(name);
       iph.makeFacet(0,0);
       const askap::scimath::Axes axes(ip.axes(iph.paramName()));
       ASKAPDEBUGASSERT(axes.has("FACETSTEP") && axes.has("STOKES") && axes.has("FREQUENCY")
                        && axes.hasDirection());
       const casa::IPosition shape = ip.value(iph.paramName()).shape();
       ASKAPDEBUGASSERT(shape.nelements()>=2);

       const int facetStep = int(axes.start("FACETSTEP"));
       ASKAPDEBUGASSERT(facetStep>0);

       // Determine the shape of the total image from the facetStep.
       // Note it is not done from the facet size because facets can overlap.
       // Normally facetStep and facet size will be equal.
       casa::IPosition newShape(shape);
       newShape[0]=facetStep*nfacets;
       newShape[1]=facetStep*nfacets;

       Axes newAxes(axes);
       casa::DirectionCoordinate dc(axes.directionAxis());
       const casa::Vector<casa::Double> refPix(2,double(newShape[0])/2.);
       dc.setReferencePixel(refPix);
       newAxes.addDirectionAxis(dc);

       casa::Array<double> pixels(newShape);
       pixels.set(0.0);
       ip.add(iph.taylorName(), pixels, newAxes);
    }

    /// @brief obtain an array corresponding to a single facet of a merged faceted image
    /// @details Each facet is represented by a number of independent parameters with
    /// the names containing .facet.x.y at the end. One of the add methods can add a
    /// parameter representing merged image (with the name without any suffixes). This
    /// method allows to translate the name of the facet (with suffixes) into a slice of
    /// the merged array corresponding to this particular facet. The suffixes are removed
    /// automatically to locate the merged image. This is the core method necessary for
    /// merging individual facets together (which happens inside ImageRestoreSolver).
    /// @param[in] ip parameters
    /// @param[in] name name of the facet parameter (with suffix like .facet.0.0)
    /// @return an array of doubles representing a subimage of the merged image
    casa::Array<double> SynthesisParamsHelper::getFacet(askap::scimath::Params &ip, const string &name)
    {
      ASKAPDEBUGASSERT(ip.has(name));
      // parse the name
      ImageParamsHelper iph(name);
      // name with the suffixes related to facets removed (and taylor suffix preserved if present)
      const std::string mergedName = iph.taylorName();
      ASKAPCHECK(ip.has(mergedName), "Merged image ("<<mergedName<<") doesn't exist");
      // there is no consistency check that the given facet corresponds to this particular
      // merged image and coordinate systems match.

      // now find blc and trc of the patch inside the big image
      const askap::scimath::Axes axes(ip.axes(mergedName));
      ASKAPDEBUGASSERT(axes.has("FACETSTEP"));
      ASKAPCHECK(casa::abs(axes.start("FACETSTEP")-axes.end("FACETSTEP"))<0.5, "facet steps extracted from "<<
                 iph.name()<<" are notably different for ra and dec axes. Should be the same integer number");
      const int facetStep = int(axes.start("FACETSTEP")+0.5);

      casa::Array<double> mergedImage = ip.value(mergedName);
      casa::IPosition blc(mergedImage.shape());
      casa::IPosition trc(mergedImage.shape());
      ASKAPDEBUGASSERT(blc.nelements()>=2);
      // adjust extra dimensions
      for (size_t i=2;i<blc.nelements();++i) {
           blc[i] = 0;
           ASKAPDEBUGASSERT(trc[i]!=0);
           trc[i] -= 1;
      }

      casa::IPosition patchShape = ip.value(name).shape();
      ASKAPDEBUGASSERT(patchShape.nelements()>=2);
      ASKAPDEBUGASSERT((facetStep<=patchShape[0]) && (facetStep<=patchShape[1]));

      ASKAPDEBUGASSERT(facetStep>=1);
      blc[0] = iph.facetX()*facetStep;
      trc[0] = blc[0]+facetStep-1;
      blc[1] = iph.facetY()*facetStep;
      trc[1] = blc[1]+facetStep-1;

      /*
      const casa::DirectionCoordinate csPatch = directionCoordinate(ip,name);
      const casa::DirectionCoordinate csFull = directionCoordinate(ip,iph.name());
      casa::Vector<double> world(2);
      // first get blc
      casa::Vector<double> blcPixel(2);
      blcPixel(0)=double((patchShape[0]-facetStep)/2);
      blcPixel(1)=double((patchShape[1]-facetStep)/2);
      std::cout<<blcPixel<<endl;
      csPatch.toWorld(world,blcPixel);
      csFull.toPixel(blcPixel,world);
      std::cout<<blcPixel<<endl;

      // now get trc
      casa::Vector<double> trcPixel(2);
      trcPixel[0]=double((patchShape[0]+facetStep)/2-1);
      trcPixel[1]=double((patchShape[1]+facetStep)/2-1);
      ASKAPDEBUGASSERT((trcPixel[0]>0) && (trcPixel[1]>0));
      std::cout<<trcPixel<<endl;
      csPatch.toWorld(world,trcPixel);
      csFull.toPixel(trcPixel,world);
      std::cout<<trcPixel<<endl;
      for (size_t dim=0;dim<2;++dim) {
           const int pix1 = int(blcPixel[dim]);
           const int pix2 = int(trcPixel[dim]);
           blc[dim] = pix1>pix2 ? pix2 : pix1;
           trc[dim] = pix1>pix2 ? pix1 : pix2;
      }
      */
      // ready to make a slice
      //std::cout<<blc<<" "<<trc<<" "<<facetStep<<" "<<mergedImage.shape()<<std::endl;
      ASKAPDEBUGASSERT((trc[0]-blc[0]+1 == facetStep) && (trc[1]-blc[1]+1 == facetStep));
      return mergedImage(blc,trc);
    }

    /// @brief A helper method to parse string of quantities
    /// @details Many parameters in parset file are given as quantities or
    /// vectors of quantities, e.g. [8.0arcsec,8.0arcsec]. This method allows
    /// to parse vector of strings corresponding to such parameter and return
    /// a vector of double values in the required units.
    /// @param[in] strval input vector of strings
    /// @param[in] unit required units (given as a string)
    /// @return vector of doubles with converted values
    std::vector<double> SynthesisParamsHelper::convertQuantity(const std::vector<std::string> &strval,
                       const std::string &unit)
    {
       std::vector<double> result(strval.size());
       std::vector<std::string>::const_iterator inIt = strval.begin();
       for (std::vector<double>::iterator outIt = result.begin(); inIt != strval.end();
                                                ++inIt,++outIt) {
            ASKAPDEBUGASSERT(outIt != result.end());
            *outIt = convertQuantity(*inIt,unit);
       }
       return result;
    }

    /// @brief A helper method to parse string of quantities
    /// @details Many parameters in parset file are given as quantities or
    /// vectors of quantities, e.g. 8.0arcsec. This method allows
    /// to parse a single string corresponding to such a parameter and return
    /// a double value converted to the requested units.
    /// @param[in] strval input string
    /// @param[in] unit required units (given as a string)
    /// @return converted value
    double SynthesisParamsHelper::convertQuantity(const std::string &strval,
                       const std::string &unit)
    {
       casa::Quantity q;

       casa::Quantity::read(q, strval);
       return q.getValue(casa::Unit(unit));
    }

    void SynthesisParamsHelper::saveImageParameter(const askap::scimath::Params& ip, const string& name,
						const string& imagename)
    {
      ASKAPTRACE("SynthesisParamsHelper::saveImageParameter");
      const casa::Array<double> imagePixels(ip.value(name));
      ASKAPDEBUGASSERT(imagePixels.ndim()!=0);
      const casa::CoordinateSystem imageCoords(coordinateSystem(ip,name));

      casa::Array<float> floatImagePixels(imagePixels.shape());
      casa::convertArray<float, double>(floatImagePixels, imagePixels);
      ASKAPLOG_DEBUG_STR(logger, "Data of "<<name<<" parameter peak at "<<casa::max(floatImagePixels));

      imageHandler().create(imagename, floatImagePixels.shape(), imageCoords);
      imageHandler().write(imagename, floatImagePixels);

      const Axes &axes = ip.axes(name);
      if (axes.has("MAJMIN")) {
          // this is a restored image with beam parameters set
          ASKAPCHECK(axes.has("PA"),"PA axis should always accompany MAJMIN");
          imageHandler().setUnits(imagename, "Jy/beam");
          imageHandler().setBeamInfo(imagename, axes.start("MAJMIN"), axes.end("MAJMIN"),
                                     axes.start("PA"));
      } else {
          if (imagename.find("sensitivity") == 0) {
               // sensitivity image is a special case, it has units of Jy/beam, but no beam info
              imageHandler().setUnits(imagename, "Jy/beam");
          }
          else if (imagename.find("residual") == 0) {
              // residual image is also a special case - perhaps the beam info should be from the fitted
              // psf
              imageHandler().setUnits(imagename, "Jy/beam");
          }
          else {
             imageHandler().setUnits(imagename, "Jy/pixel");
          }
      }
    }

    /// @brief obtain image handler
    /// @details For some operations it may be necessary to access the (global) instance of the
    /// image handler. This method allows that. An exception is thrown if no image handler has
    /// been previously set up.
    /// @return a reference to image handler
    accessors::IImageAccess& SynthesisParamsHelper::imageHandler()
    {
      ASKAPCHECK(theirImageAccessor, "setUpImageHandler has to be called before any read/write operation");
      return *theirImageAccessor;
    }

    /// @brief configure default frequency frame
    /// @details All code workes in a single frequency frame (convertions are done, if
    /// necessary when the data are read (using conversion mechanism provided by the accessor).
    /// A call to this method sets up new default.
    /// @param[in] frame reference frame to use for all created images
    void SynthesisParamsHelper::setDefaultFreqFrame(const casa::MFrequency::Ref &frame)
    {
      theirFreqFrame = frame;
    }

    /// @brief setup image handler
    /// @details This method uses the factory to setup a helper class handling the
    /// operations with images (default is casa). It is necessary to call this method
    /// at least once before any read or write operation can happen.
    /// @param[in] parset a parset file containing parameters describing which image handler to use
    /// @note The key parameter describing the image handler is "imagetype". By default, the
    /// casa image handler is created (however, a call to this method is still required)
    void SynthesisParamsHelper::setUpImageHandler(const LOFAR::ParameterSet &parset)
    {
      theirImageAccessor = accessors::imageAccessFactory(parset);
    }


    void SynthesisParamsHelper::loadImageParameter(askap::scimath::Params& ip, const string& name,
						 const string& imagename)
    {
      ASKAPTRACE("SynthesisParamsHelper::loadImageParameter");
      casa::Array<float> pixels = imageHandler().read(imagename);
      casa::Array<double> imagePixels(pixels.shape());
      casa::convertArray<double, float>(imagePixels, pixels);

      casa::CoordinateSystem imageCoords = imageHandler().coordSys(imagename);

      /// Fill in the axes information
      Axes axes;
      /// First do the direction
      int whichDir=imageCoords.findCoordinate(Coordinate::DIRECTION);
      ASKAPCHECK(whichDir>-1, "No direction coordinate present in the image "<<imagename);
      casa::DirectionCoordinate radec(imageCoords.directionCoordinate(whichDir));
      casa::Vector<casa::Int> axesDir = imageCoords.pixelAxes(whichDir);
      ASKAPCHECK(axesDir.nelements() == 2, "Direction axis "<<whichDir<<
                 " is expected to correspond to just two pixel axes, you have "<<axesDir);
      ASKAPCHECK((axesDir[0] == 0) && (axesDir[1] == 1),
               "At present we support only images with first axes being the direction pixel axes, image "<<name<<
               " has "<< axesDir);

      casa::Vector<casa::String> units(2);
      units.set("rad");
      radec.setWorldAxisUnits(units);

      /*
	  casa::Vector<double> start(2);
      casa::Vector<double> end(2);
      casa::Vector<double> pixelbuf(2);
      const casa::Vector<double> refPix(radec.referencePixel());

      pixelbuf(0)=0;
      pixelbuf(1)=refPix(1);
      ASKAPCHECK(radec.toWorld(start,pixelbuf), "Pixel to world conversion error: "<<radec.errorMessage());
      pixelbuf(0)=imagePixels.shape()[0]-1;
      ASKAPCHECK(radec.toWorld(end,pixelbuf), "Pixel to world conversion error: "<<radec.errorMessage());
      axes.add("RA", start(0), end(0));

      pixelbuf(0)=refPix(0);
      pixelbuf(1)=0;
      ASKAPCHECK(radec.toWorld(start,pixelbuf), "Pixel to world conversion error: "<<radec.errorMessage());
      pixelbuf(1)=imagePixels.shape()[1]-1;
      ASKAPCHECK(radec.toWorld(end,pixelbuf), "Pixel to world conversion error: "<<radec.errorMessage());
      axes.add("DEC", start(1), end(1));
      */
      axes.addDirectionAxis(radec);

      int whichStokes = imageCoords.findCoordinate(Coordinate::STOKES);
      int nPol = 1;
      if (whichStokes<0) {
          const casa::Vector<casa::Stokes::StokesTypes> dummyStokes(1,casa::Stokes::I);
          axes.addStokesAxis(dummyStokes);
      } else {
          casa::StokesCoordinate sc(imageCoords.stokesCoordinate(whichStokes));
          const casa::Vector<casa::Int> stokesAsInt = sc.stokes();
          casa::Vector<casa::Stokes::StokesTypes> stokes(stokesAsInt.nelements());
          for (casa::uInt pol=0; pol<stokes.nelements(); ++pol) {
               stokes[pol] = casa::Stokes::StokesTypes(stokesAsInt[pol]);
          }
          axes.addStokesAxis(stokes);

          casa::Vector<casa::Int> axesStokes = imageCoords.pixelAxes(whichStokes);
          ASKAPCHECK(axesStokes.nelements() == 1, "Stokes axis "<<whichStokes<<
                 " is expected to correspond to just one pixel axes, you have "<<axesStokes);
          ASKAPASSERT(casa::uInt(axesStokes[0])<imagePixels.shape().nelements());
          nPol = imagePixels.shape()(axesStokes[0]);
          ASKAPASSERT(uInt(nPol) == stokesAsInt.nelements());
      }

      int whichSpectral=imageCoords.findCoordinate(Coordinate::SPECTRAL);
      ASKAPCHECK(whichSpectral>-1, "No spectral coordinate present in model");
      casa::Vector<casa::Int> axesSpectral = imageCoords.pixelAxes(whichSpectral);
      ASKAPCHECK(axesSpectral.nelements() == 1, "Spectral axis "<<whichSpectral<<
                 " is expected to correspond to just one pixel axes, you have "<<axesSpectral);
      ASKAPASSERT(casa::uInt(axesSpectral[0])<imagePixels.shape().nelements());
      const int nChan = imagePixels.shape()(axesSpectral[0]);
      casa::SpectralCoordinate freq(imageCoords.spectralCoordinate(whichSpectral));
      double startFreq, endFreq;
      freq.toWorld(startFreq, 0.0);
      freq.toWorld(endFreq, double(nChan-1));
      axes.add("FREQUENCY", startFreq, endFreq);
      const casa::IPosition targetShape(4, imagePixels.shape()(0), imagePixels.shape()(1), nPol, nChan);
      ASKAPLOG_INFO_STR(logger, "About to add new image parameter with name "<<name<<
                  " reshaped to "<<targetShape<<" from original image shape "<<imagePixels.shape());
      ASKAPLOG_INFO_STR(logger, "Spectral axis will have startFreq="<<startFreq<<" Hz, endFreq="<<endFreq<<
                                "Hz, nChan="<<nChan);
      ASKAPDEBUGASSERT(targetShape.product() == imagePixels.shape().product());
      ip.add(name, imagePixels.reform(targetShape),axes);

    }

    /// @brief Get parameters corresponding to all facets from a CASA image
    /// @param[in] ip Parameters
    /// @param[in] name Base name of the parameter (.facet.x.y will be added)
    /// @param[in] fileName Base name of the image file (.facet.x.y will be added)
    /// @param[in] nfacets Number of facets on each axis (assumed the same for both axes)
    void SynthesisParamsHelper::getMultiFacetImage(askap::scimath::Params &ip, const string &name,
           const string &fileName, const int nfacets)
    {
      ASKAPCHECK(nfacets>0, "The number of facets is supposed to be positive, you have "<<nfacets);
      // create helper in two steps because the name may represent a Taylor term
      ImageParamsHelper iph(name);
      for (int ix=0; ix<nfacets; ++ix) {
           for (int iy=0; iy<nfacets; ++iy) {
                // assign facet indices to the helper
                iph.makeFacet(ix,iy);
                const std::string paramName = iph.paramName();
                loadImageParameter(ip,paramName,fileName);
           }
      }
    }

    boost::shared_ptr<casa::TempImage<float> >
    SynthesisParamsHelper::tempImage(const askap::scimath::Params& ip,
				     const string& name)
    {
      const casa::Array<double> imagePixels(ip.value(name));

      casa::CoordinateSystem imageCoords(coordinateSystem(ip, name));

      boost::shared_ptr<casa::TempImage<float> >
	im(new casa::TempImage<float> (TiledShape(imagePixels.shape()),
				       imageCoords));

      im->setUnits("Jy/pixel");

      casa::Array<float> floatImagePixels(imagePixels.shape());
      casa::convertArray<float, double>(floatImagePixels, imagePixels);
      casa::ArrayLattice<float> latImagePixels(floatImagePixels);
      im->copyData(latImagePixels);
      return im;
    }

    casa::CoordinateSystem
    SynthesisParamsHelper::coordinateSystem(const askap::scimath::Params& ip,
					    const string& name)
    {
      const Axes axes(ip.axes(name));

      casa::DirectionCoordinate radec(directionCoordinate(ip, name));

      casa::CoordinateSystem imageCoords;
      imageCoords.addCoordinate(radec);

      // default is a dummy stokes coordinate with only stokes I present
      casa::Vector<int> iquv(1);
      iquv(0) = Stokes::I;
      if (axes.has("STOKES")) {
          casa::Vector<casa::Stokes::StokesTypes> stokes = axes.stokesAxis();
          ASKAPDEBUGASSERT(stokes.nelements()>=1);
          iquv.resize(stokes.nelements());
          for (size_t pol=0; pol<stokes.nelements(); ++pol) {
               iquv[pol] = int(stokes[pol]);
          }
      }

      casa::StokesCoordinate stokes(iquv);
      imageCoords.addCoordinate(stokes);

      const casa::IPosition shape = ip.value(name).shape();
      const int nchan = shape.nelements() >= 4 ? shape(3) : 1;
      const double restfreq = 0.0;
      const double crpix = double(nchan-1)/2.;
      const double crval = (axes.start("FREQUENCY")+axes.end("FREQUENCY"))/2.0;
      // we can't estimate increment if there is only one channel and start=stop
      const double cdelt = nchan>1 ? (axes.end("FREQUENCY")-axes.start("FREQUENCY"))/double(nchan-1) : 1.;
      const casa::SpectralCoordinate freq(casa::MFrequency::castType(theirFreqFrame.getType()), crval, cdelt, crpix, restfreq);
      imageCoords.addCoordinate(freq);

      return imageCoords;
    }

    casa::DirectionCoordinate
    SynthesisParamsHelper::directionCoordinate(const askap::scimath::Params& ip,
					       const string& name)
    {
      const Axes axes(ip.axes(name));
      ASKAPCHECK(!axes.has("RA-TANGENT") && !axes.has("DEC-TANGENT"),
          "Detected an obsolete way to specify tangent point!");

      ASKAPCHECK(axes.hasDirection(), "Direction coordinate is missing. axes["<<name<<"]:"<<axes);
      return axes.directionAxis();

    }

    void SynthesisParamsHelper::update(askap::scimath::Params& ip, const string& name,
				       const casa::ImageInterface<float>& im)
    {
      ASKAPDEBUGTRACE("SynthesisParamsHelper::update");
      /// This next copy should be a reference unless it is too big
      casa::Array<float> floatImagePixels(im.shape());
      casa::ArrayLattice<float> latImagePixels(floatImagePixels);
      latImagePixels.copyData(im);

      casa::Array<double> imagePixels(im.shape());
      casa::convertArray<double, float>(imagePixels, floatImagePixels);
      ip.update(name, imagePixels);
    }

    /// @brief check whether parameter list defines at least one component
    /// @details Parameter lists can have a mixture of components and
    /// images defined. This method checks whether the given parameter
    /// list defines at least one component.
    /// @param[in] params a shared pointer to the parameter container
    /// @return true, if at least one component is defined
    bool SynthesisParamsHelper::hasComponent(const askap::scimath::Params::ShPtr &params)
    {
       ASKAPDEBUGASSERT(params);
       return (params->completions("flux.i").size()!=0) || (params->completions("calibrator.").size()!=0);
    }

    /// @brief check whether parameter list defines at least one image
    /// @details Parameter lists can have a mixture of components and
    /// images defined. This method checks whether the given parameter
    /// list defines at least one image.
    /// @param[in] params a shared pointer to the parameter container
    /// @return true, if at least one image is defined
    bool SynthesisParamsHelper::hasImage(const askap::scimath::Params::ShPtr &params)
    {
       ASKAPDEBUGASSERT(params);
       return params->completions("image").size()!=0;
    }

    /// @brief find a slicer matching another direction coordinate
    /// @details This method builds BLC and TRC of an image corresponding
    /// to the provided direction coordinate which cover the same area as
    /// the given image parameter (the functionality is similar to that of
    /// getFacet, although this interface allows more flexibility).
    /// @param[in] ip parameters
    /// @param[in] name image parameter to map
    /// @param[in] dc direction coordinate
    /// @return slicer encapsulating BLC and TRC
    /// @note The dimensionality of the output corresponds to the input image
    casa::Slicer SynthesisParamsHelper::facetSlicer(const askap::scimath::Params& ip,
           const std::string &name, const casa::DirectionCoordinate &dc)
    {
       ASKAPDEBUGASSERT(ip.has(name));
       const casa::DirectionCoordinate inDC = directionCoordinate(ip,name);
       const casa::IPosition inShape = ip.value(name).shape();
       ASKAPDEBUGASSERT(inShape.nelements() >= 2);
       casa::IPosition blc(inShape.nelements(),0);
       casa::IPosition trc(inShape);
       for (casa::uInt dim=0; dim<inShape.nelements(); ++dim) {
            --trc(dim);
       }
       // currently blc,trc describe the whole input image; convert coordinates
       casa::Vector<casa::Double> pix(2);

       // first process BLC
       pix[0] = casa::Double(blc[0]);
       pix[1] = casa::Double(blc[1]);
       casa::MDirection tempDir;
       casa::Bool success = inDC.toWorld(tempDir, pix);
       ASKAPCHECK(success, "Pixel to world coordinate conversion failed for input BLC: "<<inDC.errorMessage());
       success = dc.toPixel(pix,tempDir);
       ASKAPCHECK(success, "World to pixel coordinate conversion failed for output BLC: "<<dc.errorMessage());
       // converting to Int without rounding appears to be changing the image size.
       //blc[0] = casa::Int(pix[0]);
       //blc[1] = casa::Int(pix[1]);
       blc[0] = casa::Int(round(pix[0]));
       blc[1] = casa::Int(round(pix[1]));

       // now process TRC
       pix[0] = casa::Double(trc[0]);
       pix[1] = casa::Double(trc[1]);
       success = inDC.toWorld(tempDir, pix);
       ASKAPCHECK(success, "Pixel to world coordinate conversion failed for input TRC: "<<inDC.errorMessage());
       success = dc.toPixel(pix,tempDir);
       ASKAPCHECK(success, "World to pixel coordinate conversion failed for output TRC: "<<dc.errorMessage());
       // converting to Int without rounding appears to be changing the image size.
       //trc[0] = casa::Int(pix[0]);
       //trc[1] = casa::Int(pix[1]);
       trc[0] = casa::Int(round(pix[0]));
       trc[1] = casa::Int(round(pix[1]));

       return casa::Slicer(blc,trc,casa::Slicer::endIsLast);
    }

    /// @brief make a merged image parameter covering all given facets
    /// @details This method is very similar to another version of the add method which creates
    /// an image paramters covering all facets named in the appropriate fashion. Although doing
    /// essentially the same job, this method works with any images, i.e. they are not necessarily
    /// regularly spaced and appropriately named facets. With time we can probably change how we
    /// do faceting and retire the old methods.
    /// @param[in] ip parameters
    /// @param[in] names names of all images to merge
    /// @param[in] mergedName name of the image to create
    void SynthesisParamsHelper::add(askap::scimath::Params& ip, const std::vector<std::string> &names,
              const std::string &mergedName)
    {
       ASKAPCHECK(names.size()>0, "At least one input image is expected by SynthesisParamsHelper::add");
       const casa::DirectionCoordinate templateDC = directionCoordinate(ip,names[0]);
       // we could've generated initial slicer from the shape and avoid one unnecessary call to facetSlicer
       const casa::Slicer tempSlicer = facetSlicer(ip,names[0],templateDC);
       casa::IPosition tempBLC = tempSlicer.start();
       casa::IPosition tempTRC = tempSlicer.end();
       ASKAPDEBUGASSERT(tempBLC.nelements() >= 2);
       ASKAPDEBUGASSERT(tempTRC.nelements() >= 2);
       for (size_t i = 1; i<names.size(); ++i) {
            const casa::Slicer newSlicer = facetSlicer(ip,names[i],templateDC);
            const casa::IPosition newBLC = newSlicer.start();
            const casa::IPosition newTRC = newSlicer.end();
            ASKAPDEBUGASSERT(newBLC.nelements() >= 2);
            ASKAPDEBUGASSERT(newTRC.nelements() >= 2);
            for (casa::uInt dim=0; dim<2; ++dim) {
                 if (newBLC(dim) < tempBLC(dim)) {
                     tempBLC(dim) = newBLC(dim);
                 }
                 if (newTRC(dim) > tempTRC(dim)) {
                     tempTRC(dim) = newTRC(dim);
                 }
            }
       }
       casa::IPosition newShape = ip.value(names[0]).shape();
       ASKAPDEBUGASSERT(newShape.nelements() >= 2);
       newShape(0) = tempTRC(0) - tempBLC(0) + 1;
       newShape(1) = tempTRC(1) - tempBLC(1) + 1;
       ASKAPDEBUGASSERT(newShape(0) > 0);
       ASKAPDEBUGASSERT(newShape(1) > 0);
       casa::Vector<casa::Double> refPix = templateDC.referencePixel();
       refPix[0] -= casa::Double(tempBLC(0) - tempSlicer.start()(0));
       refPix[1] -= casa::Double(tempBLC(1) - tempSlicer.start()(1));
       casa::DirectionCoordinate newDC(templateDC);
       newDC.setReferencePixel(refPix);

       askap::scimath::Axes newAxes(ip.axes(names[0]));
       newAxes.addDirectionAxis(newDC);

       casa::Array<double> pixels(newShape);
       pixels.set(0.0);
       ip.add(mergedName, pixels, newAxes);

    }

    /// @brief A helper method to build a list of faceted images
    /// @details All multi-facet images are split between a number of
    /// parameters named like "image.i.fieldname.facet.0.0". Single
    /// facet images correspond to parameters named like "image.i.fieldname".
    /// This method reads a supplied vector of names (may be either all names
    /// or just free parameters extracted from Params object) and builds a map
    /// of the image name (up to and including fieldname) and the number of
    /// facets. It also does the necessary checks that all required facets are
    /// defined and throws an exception if it is not the case.
    /// @param[in] names parameter names to work with
    /// @param[out] facetmap a map of (possibly truncated names) and the number of facets
    /// @note 1. facetmap.size()<=names.size() after the call to this method
    /// 2. This method just adds the content to the facet map without erasing the
    /// existing information.
    void SynthesisParamsHelper::listFacets(const std::vector<std::string> &names,
                          std::map<std::string, int> &facetmap)
    {
       // temporary maps, just to check that no facets were missed
       std::map<std::string, std::set<int> >  tempMapX;
       std::map<std::string, std::set<int> >  tempMapY;

       for (std::vector<std::string>::const_iterator ci = names.begin(); ci!=names.end(); ++ci) {
            ImageParamsHelper iph(*ci);
            // name with the facet-related suffixes removed (and taylor suffix preserved, if present)
            const std::string baseName = iph.taylorName();
            if (iph.isFacet()) {
                tempMapX[baseName].insert(iph.facetX());
                tempMapY[baseName].insert(iph.facetY());
                facetmap[baseName] = 0; // a flag that we need to figure out the exact number later
            } else {
                // this is not a faceted image, just add it to the final list
                facetmap[baseName] = 1; // one facet
            }
       }
       for (std::map<std::string, int>::iterator it = facetmap.begin(); it!=facetmap.end(); ++it) {
            if (it->second == 0) {
                ASKAPDEBUGASSERT(hasValue(tempMapX, it->first));
                ASKAPDEBUGASSERT(hasValue(tempMapY, it->first));

                // the code below assumes equal number of facets in both axes. It should be
                // modified slightly to lift this restriction.

                ASKAPDEBUGASSERT(tempMapX[it->first].size());
                ASKAPDEBUGASSERT(tempMapY[it->first].size());


                const int maxFacetX = *(std::max_element(tempMapX[it->first].begin(),
                                      tempMapX[it->first].end()));
                const int maxFacetY = *(std::max_element(tempMapY[it->first].begin(),
                                      tempMapY[it->first].end()));
                const int nFacets = (maxFacetX > maxFacetY ? maxFacetX : maxFacetY)+1;

                // doing checks
                for (int facet = 0; facet<nFacets; ++facet) {
                     ASKAPCHECK(hasValue(tempMapX[it->first],facet), "Facet "<<facet<<
                          " is missing for the first axis");
                     ASKAPCHECK(hasValue(tempMapY[it->first],facet), "Facet "<<facet<<
                          " is missing for the second axis");
                }

                it->second = nFacets;
            }
       }
    }

    /// @brief A helper method to build a list of images representing Taylor terms
    /// @details Different Taylor terms in the multi-frequency algorithm are
    /// represented by parameters named like "image.fieldname.taylor.0". This method
    /// reads a supplied vector of names (may be just free parameters or all names
    /// available) and builds a map of the actual image name (without suffixes) and
    /// the number of Taylor orders encountered. It also does check that all orders
    /// starting from 0 are present and throws an exception if it is not the case.
    /// To some extent this method is similar to listFacets, but is intended for
    /// Taylor terms.
    /// @param[in] names parameter names to work with
    /// @param[out] taylormap a map of (possibly truncated names) and the number of
    ///             Taylor terms (1 means no decomposition into Taylor series, i.e. no MFS)
    /// @note 1. taylormap.size()<=names.size() after a call to this method, if it was originally
    ///          empty
    ///       2. This method just adds new elements to the taylormap without erasing the
    ///          existing information.
    void SynthesisParamsHelper::listTaylor(const std::vector<std::string> &names,
                           std::map<std::string, int> &taylormap)
    {
       // temporary map, just to check that no Taylor terms are missed
       // (parameters may not come in any particular order)
       std::map<std::string, std::set<int> >  tempMap;

       for (std::vector<std::string>::const_iterator ci = names.begin(); ci!=names.end(); ++ci) {
            ImageParamsHelper iph(*ci);
            // name with the taylor-related suffixes removed (and facet suffixes preserved, if present)
            const std::string baseName = iph.facetName();
            if (iph.isTaylorTerm()) {
                // this is a Taylor term, we need to remember all orders sited for this base name
                tempMap[baseName].insert(iph.order());
                taylormap[baseName] = 0; // just a flag, we need to figure out the exact number later
            } else {
                // this is not an MFS'ed image, add it to the final list
                taylormap[baseName] = 1; // single order
            }
       }

       for (std::map<std::string, int>::iterator it = taylormap.begin(); it!=taylormap.end(); ++it) {
            if (it->second == 0) {
                // This is the MFS case, need to figure out the exact number of Taylor terms
                ASKAPDEBUGASSERT(hasValue(tempMap, it->first));
                ASKAPDEBUGASSERT(tempMap[it->first].size());

                const int nTaylorTerms = *(std::max_element(tempMap[it->first].begin(),
                                      tempMap[it->first].end())) + 1;

                // doing checks
                for (int order = 0; order<nTaylorTerms; ++order) {
                     ASKAPCHECK(hasValue(tempMap[it->first],order), "Taylor term "<<order<<
                          " is missing for the image "<<it->first);
                }

                it->second = nTaylorTerms;
            }
       }
    }
    /// @brief find a parameter representing a PSF
    /// @details If multiple PSF parameters are present, the first encountered is returned
    /// @param[in] ip parameters
    /// @return full name of some PSF parameter or an empty string if it is not found
    std::string SynthesisParamsHelper::findPSF(const askap::scimath::Params &ip)
    {
        std::string psfName;
        // Find all the free parameters beginning with image
        vector<string> names(ip.completions("image"));
        for (vector<string>::const_iterator it = names.begin(); it!=names.end(); ++it) {
             std::string curName = "psf.image";
             if (ip.has(curName + *it)) {
                 if (psfName != "") {
                     ASKAPLOG_WARN_STR(logger, "Multiple PSF parameters are present, using "<<
                              psfName<<" which was first encountered");
                     break;
                 }
                 psfName = curName + *it;
             } else if (ip.has(string("psf") + *it)) {
                 if (psfName != "") {
                     ASKAPLOG_WARN_STR(logger, "Multiple PSF parameters are present, using "<<
                             psfName<<" which was first encountered");
                     break;
                 }
                 psfName = string("psf") + *it;
             }
        }
        return psfName;
    }

    /// @brief fit gaussian beam into PSF
    /// @details This method fits a 2D Gaussian into the given PSF image. If no parameter
    /// name is given (i.e. an empty string is passed to this method), the most appropriate
    /// parameter is automatically selected (i.e. psf.image.something if preconditioning is
    /// done and psf.something if not). First match is always used. If the image is
    /// multi-dimensional, only first plane is used. A warning is given in the case of
    /// a potential ambiguity.
    /// @param[in] ip parameters
    /// @param[in] cutoff cutoff defining the support size where the fitting is done (default
    ///            is 0.05, i.e. fitting is done to pixels enclosed in a rectangular support
    ///            defined by 5% cutoff from the peak)
    /// @param[in] name full name of the parameter representing the PSF (default is to figure this out)
    casa::Vector<casa::Quantum<double> > SynthesisParamsHelper::fitBeam(const askap::scimath::Params &ip,
                                     const double cutoff, const std::string &name)
    {
       ASKAPTRACE("SynthesisParamsHelper::fitBeam");

       std::string psfName = name;
       if (name == "") {
           // we have to figure out the name
           psfName = findPSF(ip);
       }
       ASKAPCHECK(psfName != "", "Unable to find psf paramter to fit, params="<<ip);
       ASKAPLOG_INFO_STR(logger, "Fitting 2D Gaussian into PSF parameter "<<psfName);

       casa::Array<double> psfArray = ip.value(psfName);
       return fitBeam(psfArray, ip.axes(psfName), cutoff);

    }

    casa::Vector<casa::Quantum<double> > SynthesisParamsHelper::fitBeam(casa::Array<double> &psfArray,
                                                                        const scimath::Axes &axes,
                                                                        const double cutoff) {

       const casa::IPosition shape = psfArray.shape();
       ASKAPCHECK(shape.nelements()>=2,"PSF image is supposed to be at least 2-dimensional, shape="<<psfArray.shape());
       if (shape.product() != shape[0]*shape[1]) {
           ASKAPLOG_WARN_STR(logger, "Multi-dimensional PSF is present (shape="<<shape<<
                             "), using the first 2D plane only to fit the beam");
       }

       casa::Array<double> psfSlice = MultiDimArrayPlaneIter::getFirstPlane(psfArray).nonDegenerate();
       ASKAPDEBUGASSERT(psfSlice.shape().nelements() == 2);
       casa::Matrix<double> psfSliceMatrix = psfSlice;

       // search for support to speed up beam fitting
       ASKAPLOG_INFO_STR(logger, "Searching for support with the relative cutoff of "<<cutoff<<" to speed fitting up");
       SupportSearcher ss(cutoff);
       ss.search(psfSliceMatrix);
       casa::uInt support = ss.symmetricalSupport(psfSlice.shape());
       support = 2 * (support/2) + 1; // if even, move up to next odd.
       
       ASKAPLOG_INFO_STR(logger, "Extracting support of "<<support<<" pixels for 2D gaussian fitting");
       const casa::IPosition newShape(2,support,support);
       for (int dim=0; dim<2; ++dim) {
            ASKAPCHECK(psfSlice.shape()[dim] >= int(support), "Support is greater than the original size, shape="<<
                       psfSlice.shape());
       }
       //
       casa::Array<float> floatPSFSlice(newShape);
       casa::convertArray<float, double>(floatPSFSlice, scimath::PaddingUtils::centeredSubArray(psfSlice,newShape));

       // hack for debugging only
       //floatPSFSlice = imageHandler().read("tmp.img").nonDegenerate();
       //

       // normalise to 1
       const float maxPSF = casa::max(floatPSFSlice);
       if (fabs(maxPSF-1.)>1e-6) {
           floatPSFSlice /= maxPSF;
       }
       //
       
       // actual fitting
       casa::Vector<casa::Double> initialEstimate(6,0.);
       initialEstimate[0]=1.; // PSF peak is always 1
       initialEstimate[1]=newShape[0]/2; // centre
       initialEstimate[2]=newShape[1]/2; // centre
       initialEstimate[3]=1;  // 1 pixel wide
       initialEstimate[4]=0.9;  // 1 pixel wide
       initialEstimate[5]=casa::C::pi/4.; // quire arbitrary  pa.
       casa::Vector<casa::Bool> parameterMask(6,casa::False);
       parameterMask[3] = casa::True; // fit maj
       parameterMask[4] = casa::True; // fit min
       parameterMask[5] = casa::True; // fit pa

       casa::LogIO os;
       casa::Fit2D fitter(os);
       fitter.addModel(casa::Fit2D::GAUSSIAN,initialEstimate,parameterMask);
       casa::Array<casa::Float> sigma(floatPSFSlice.shape(),1.);
       const casa::Fit2D::ErrorTypes fitError = fitter.fit(floatPSFSlice,sigma);
       ASKAPCHECK(fitError == casa::Fit2D::OK, "Error fitting the beam. fitError="<<fitError<<
                  " message: "<<fitter.errorMessage());
       casa::Vector<casa::Double> result = fitter.availableSolution();
       ASKAPLOG_INFO_STR(logger, "Got fit result (in pixels) "<<result<<" and uncertainties "<<fitter.availableErrors());
       ASKAPCHECK(result.nelements() == 6, "Expect 6 parameters for 2D gaussian, result vector has "<<result.nelements());
       ASKAPCHECK(axes.hasDirection(), "Direction axes are missing from the PSF parameter, unable to convert pixels to angular units");
       const casa::Vector<casa::Double> increments = axes.directionAxis().increment();
       ASKAPCHECK(increments.nelements() == 2, "Expect just two elements for increments of the direction axis, you have "<<
                  increments);
       ASKAPCHECK(increments[1]>0, "Expect positive increment on the declination axis. increments="<<increments);
       ASKAPCHECK(fabs(fabs(increments[0])-fabs(increments[1]))<1e-6,
                  "Different cell sizes mean that the current beam fitting code would give a wrong position angle. increments="
                  <<increments);
       casa::Vector<casa::Quantum<double> > beam(3);
       beam[0] = casa::Quantum<double>(fabs(increments[0])*result[3],"rad");
       beam[1] = casa::Quantum<double>(fabs(increments[1])*result[4],"rad");
       // position angle in radians
       double pa = increments[0]<0 ? result[5] - casa::C::pi/2 : casa::C::pi/2 - result[5];
       if (pa < -casa::C::pi/2) {
           pa += casa::C::pi;
       }
       beam[2] = casa::Quantum<double>(pa,"rad");
       return beam;
    }


    /// @brief zero all free model images
    /// @details I (MV) hope that this method is temporary. In the current design of the code we need to
    /// discard the solution of the model updates in the case of a dirty image. Otherwise, the restored image
    /// is wrong. This method iterates over all free model image parameters and sets them to 0.
    /// @param[in] params collection of parameters
    void SynthesisParamsHelper::zeroAllModelImages(const askap::scimath::Params::ShPtr& params) {
         ASKAPDEBUGASSERT(params);
         // Find all the free parameters beginning with image
         vector<string> names(params->completions("image"));
         for (vector<string>::const_iterator it=names.begin(); it!=names.end(); ++it) {
              const std::string name="image"+*it;
              if (params->isFree(name)) {
                  params->value(name).set(0.);
              }
         }
    }


    /// @brief load component-related parameters from a parset file
    /// @details Parameter layout is different in scimath::Params and
    /// parset files for some reason. Typically a source is defined with
    /// parameters like flux.i.name, direction.ra.name, ... within the
    /// scimath::Params, but in the parset file the names of the parameters
    /// are sources.name.flux.i, sources.name.direction.ra, etc). This
    /// method translates the parameter names and copies the values accross.
    /// @param[in] params a shared pointer to the parameter container
    /// @param[in] parset a parset object to read the data from
    /// @param[in] srcName name of the source
    /// @param[in] baseKey a prefix added to parset parameter names (default
    /// is "sources.", wich matches the current layout of the parset file)
    void SynthesisParamsHelper::copyComponent(const askap::scimath::Params::ShPtr &params,
           const LOFAR::ParameterSet &parset,
           const std::string &srcName, const std::string &baseKey)
    {
       ASKAPDEBUGASSERT(params);
       // check the special case of predefined calibrators
       const std::string calParamName = baseKey + srcName + ".calibrator";
       if (parset.isDefined(calParamName)) {
           params->add("calibrator."+parset.getString(calParamName));
           ASKAPCHECK(params->completions("calibrator.").size() == 1,
              "It is not intended to have two pre-defined calibrators in the same model simultaneously. params: "<<*params);
           return; // no need to load individual parameters as they're not used in this case
       }

       // load explicitly specified component

       // first, create a list of parameters describing the component
       // if the value of the map is true, the parameter is mandatory
       // (in the future we may have a more flexible code here filling this map)
       std::map<std::string, bool>  parameterList;
       parameterList["flux.i"] = true;
       parameterList["direction.ra"] = true;
       parameterList["direction.dec"] = true;
       parameterList["shape.bmaj"] = false;
       parameterList["shape.bmin"] = false;
       parameterList["shape.bpa"] = false;

       // now iterate through all parameters
       for (std::map<std::string, bool>::const_iterator ci = parameterList.begin();
            ci!=parameterList.end(); ++ci) {
            const std::string parName = baseKey+srcName+"."+ci->first;
            if (parset.isDefined(parName)) {
                const double val = parset.getDouble(parName);
                params->add(ci->first+"."+srcName, val);
            } else {
                if (ci->second) {
                    ASKAPTHROW(AskapError, "Parameter "<<parName<<
                           " is required to define the source "<<srcName<<
                           ", baseKey="<<baseKey<<" or "<<baseKey + srcName +".calibrator parameter should be present");
                }
            }
       }
    }
  }
}
