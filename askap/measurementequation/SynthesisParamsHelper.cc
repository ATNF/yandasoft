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

#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/imageaccess/ImageAccessFactory.h>
#include <askap/scimath/utils/PolConverter.h>
#include <askap/scimath/fitting/Axes.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/imagemath/linmos/LinmosImageRegrid.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/gridding/SupportSearcher.h>
#include <askap/scimath/fft/FFTUtils.h>
#include <casacore/lattices/LatticeMath/Fit2D.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.synthesisparamshelper");

#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Quanta.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/scimath/Mathematics/Interpolate2D.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/LinearCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>

#include <askap/scimath/fft/FFT2DWrapper.h>

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
    boost::shared_ptr<accessors::IImageAccess<casacore::Float> > SynthesisParamsHelper::theirImageAccessor;

    /// @brief default frequency frame
    casacore::MFrequency::Ref SynthesisParamsHelper::theirFreqFrame(casacore::MFrequency::TOPO);

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
              casacore::Vector<casacore::Stokes::StokesTypes> stokes = scimath::PolConverter::fromString(stokesStr);

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
		               const int facetstep = parset.getInt32(*it+".facetstep",casacore::min(shape[0],shape[1]));
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

      // Check whether or not images will need to be downsampled after loading
      boost::optional<float> extraOS;
      if (parset.isDefined("extraoversampling")) {
          extraOS = parset.getFloat("extraoversampling");
          ASKAPDEBUGASSERT(*extraOS > 1.);
      }

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
                       if (extraOS) {
                           SynthesisParamsHelper::loadImageParameter(*params,iph.paramName(),iph.paramName(),*extraOS);
                       } else {
                           SynthesisParamsHelper::loadImageParameter(*params,iph.paramName(),iph.paramName());
                       }
                   } else {
                       ASKAPLOG_INFO_STR(logger, "Loading multi-facet image image "<<iph.paramName());
                       SynthesisParamsHelper::getMultiFacetImage(*params,iph.paramName(),iph.paramName(),nfacets);
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
    /// @return casacore::Projection class
    casacore::Projection SynthesisParamsHelper::getProjection(const bool ewprojection, const double dec)
    {
       if (ewprojection) {
           const double sdec = sin(dec);
           ASKAPCHECK(sdec != 0., "Singular SCP/NCP projection dec="<<dec/casacore::C::pi*180.<<
                      " deg, sin(dec)="<<sdec<<". Use plain SIN projection instead!");
           casacore::Vector<casacore::Double> projectionParameters(2,0.);
           projectionParameters(1) = cos(dec) / sdec;
           return casacore::Projection(casacore::Projection::SIN, projectionParameters);
       }
       return casacore::Projection(casacore::Projection::SIN);
    }


    void SynthesisParamsHelper::add(askap::scimath::Params& ip,
				    const string& name, const vector<string>& direction,
				    const vector<string>& cellsize, const vector<int>& shape,
				    const bool ewprojection,
				    const double freqmin, const double freqmax, const int nchan,
				    const casacore::Vector<casacore::Stokes::StokesTypes> &stokes,
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

      const MDirection radec = asMDirection(direction);
      const double ra = radec.getValue().getLong();
      const double dec = radec.getValue().getLat();

      const double xcellsize =-1.0*convertQuantity(cellsize[0],"rad");
      const double ycellsize = convertQuantity(cellsize[1],"rad");

      /// @todo Do something with the frame info in direction[2]
      Axes axes;
      casacore::Matrix<double> xform(2,2,0.);
      xform.diagonal() = 1.;
      // direction coordinate corresponding to the case with tangent point in the image centre
      const casacore::DirectionCoordinate dcTangent(casacore::MDirection::J2000,
                  getProjection(ewprojection, dec), ra,dec,xcellsize,ycellsize,xform,nx/2,ny/2);
      if (centreDir.size()) {
          ASKAPCHECK(centreDir[2] == "J2000", "Only J2000 is implemented at the moment, you have requested "<<centreDir[2]);
          ASKAPLOG_INFO_STR(logger, "Image parameter "<<name<<" have tangent point "<<direction<<" and image centre "<<centreDir);
          const casacore::MVDirection centre = asMDirection(centreDir).getValue();
          casacore::Vector<casacore::Double> pix;
          dcTangent.toPixel(pix,centre);
          ASKAPDEBUGASSERT(pix.nelements() == 2);
          axes.addDirectionAxis(casacore::DirectionCoordinate(casacore::MDirection::J2000,
                  getProjection(ewprojection, dec), ra,dec,xcellsize,ycellsize,xform,double(nx)-pix[0],double(ny)-pix[1]));
      } else {
        ASKAPLOG_INFO_STR(logger, "Image parameter "<<name<<" have tangent point "<<direction<<" at the image centre");
        axes.addDirectionAxis(dcTangent);
      }
      axes.addStokesAxis(stokes);

      casacore::Array<imtype> pixels(casacore::IPosition(4, nx, ny, stokes.nelements(), nchan));
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
       const casacore::Vector<casacore::Stokes::StokesTypes> &stokes,
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

      const MDirection radec = asMDirection(direction);
      const double ra = radec.getValue().getLong();
      const double dec = radec.getValue().getLat();

      const double xcellsize =-1.0*convertQuantity(cellsize[0],"rad");
      const double ycellsize = convertQuantity(cellsize[1],"rad");

      // zero-filled array is the same for all facets as it is copied inside Params
      // class
      casacore::Array<imtype> pixels(casacore::IPosition(4, nx, ny, stokes.nelements(), nchan));
      pixels.set(0.0);

      // void linear transform used to set up coordinate system
      casacore::Matrix<double> xform(2,2,0.);
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
                axes.addDirectionAxis(casacore::DirectionCoordinate(casacore::MDirection::J2000,
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

                // for debugging
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
       casacore::Array<imtype> pixels = ip.valueT(name);
       const casacore::IPosition shape = pixels.shape();
       ASKAPDEBUGASSERT(shape.nelements()>=2);
       casacore::IPosition end(shape);
       for (uint index=0;index<end.nelements();++index) {
            ASKAPDEBUGASSERT(end[index]>=1);
            end[index]--;
       }

       if (shape[0]>facetStep+1) {
           // need clipping along the first axis
           casacore::IPosition start(shape.nelements(),0);
           end[0] = (shape[0]-facetStep)/2-1;
           end[1] = shape[1]-1; // although this step is strictly speaking unnecessary
           pixels(start,end).set(0.);

           end[0] = shape[0]-1;
           start[0] = (shape[0]+facetStep)/2;
           pixels(start,end).set(0.);
       }

       if (shape[1]>facetStep+1) {
           // need clipping along the second axis
           casacore::IPosition start(shape.nelements(),0);
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
                            const casacore::Vector<casacore::Quantum<double> > &beam)
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
       const casacore::IPosition shape = ip.shape(iph.paramName());
       ASKAPDEBUGASSERT(shape.nelements()>=2);

       const int facetStep = int(axes.start("FACETSTEP"));
       ASKAPDEBUGASSERT(facetStep>0);

       // Determine the shape of the total image from the facetStep.
       // Note it is not done from the facet size because facets can overlap.
       // Normally facetStep and facet size will be equal.
       casacore::IPosition newShape(shape);
       newShape[0]=facetStep*nfacets;
       newShape[1]=facetStep*nfacets;

       Axes newAxes(axes);
       casacore::DirectionCoordinate dc(axes.directionAxis());
       const casacore::Vector<casacore::Double> refPix(2,double(newShape[0])/2.);
       dc.setReferencePixel(refPix);
       newAxes.addDirectionAxis(dc);

       casacore::Array<imtype> pixels(newShape);
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
    casacore::Array<imtype> SynthesisParamsHelper::getFacet(askap::scimath::Params &ip, const string &name)
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
      ASKAPCHECK(casacore::abs(axes.start("FACETSTEP")-axes.end("FACETSTEP"))<0.5, "facet steps extracted from "<<
                 iph.name()<<" are notably different for ra and dec axes. Should be the same integer number");
      const int facetStep = int(axes.start("FACETSTEP")+0.5);

      casacore::Array<imtype> mergedImage = ip.valueT(mergedName);
      casacore::IPosition blc(mergedImage.shape());
      casacore::IPosition trc(mergedImage.shape());
      ASKAPDEBUGASSERT(blc.nelements()>=2);
      // adjust extra dimensions
      for (size_t i=2;i<blc.nelements();++i) {
           blc[i] = 0;
           ASKAPDEBUGASSERT(trc[i]!=0);
           trc[i] -= 1;
      }

      casacore::IPosition patchShape = ip.shape(name);
      ASKAPDEBUGASSERT(patchShape.nelements()>=2);
      ASKAPDEBUGASSERT((facetStep<=patchShape[0]) && (facetStep<=patchShape[1]));

      ASKAPDEBUGASSERT(facetStep>=1);
      blc[0] = iph.facetX()*facetStep;
      trc[0] = blc[0]+facetStep-1;
      blc[1] = iph.facetY()*facetStep;
      trc[1] = blc[1]+facetStep-1;

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
       casacore::Quantity q;

       casacore::Quantity::read(q, strval);
       return q.getValue(casacore::Unit(unit));
    }

    void SynthesisParamsHelper::saveImageParameter(const askap::scimath::Params& ip, const string& name,
					 const string& imagename, const boost::optional<float> extraOversampleFactor,
                     const LOFAR::ParameterSet & keywords, const std::vector<std::string>& historyLines)
    {
      ASKAPTRACE("SynthesisParamsHelper::saveImageParameter");

      casacore::Array<float> imagePixels;
      casacore::CoordinateSystem imageCoords(coordinateSystem(ip,name));

      if (extraOversampleFactor) {
          ASKAPCHECK(*extraOversampleFactor > 1.,"Oversampling factor should be > 1, not "<<*extraOversampleFactor);
          imagePixels.resize(scimath::PaddingUtils::paddedShape(ip.shape(name),*extraOversampleFactor));
          #ifdef ASKAP_FLOAT_IMAGE_PARAMS
          scimath::PaddingUtils::fftPad(ip.valueF(name),imagePixels);
          #else
          casa::Array<imtype> imageTmp(imagePixels.shape());
          scimath::PaddingUtils::fftPad(ip.valueT(name),imageTmp);
          casa::convertArray<float, imtype>(imagePixels, imageTmp);
          #endif

          // update the coordinate system for the new resolution
          const int whichDir = imageCoords.findCoordinate(Coordinate::DIRECTION);
          ASKAPCHECK(whichDir>-1, "No direction coordinate present in the image "<<name);
          casacore::DirectionCoordinate radec(imageCoords.directionCoordinate(whichDir));
          casacore::Vector<casacore::Int> axesDir = imageCoords.pixelAxes(whichDir);
          ASKAPCHECK(axesDir.nelements() == 2, "Direction axis "<<whichDir<<
                     " is expected to correspond to just two pixel axes, you have "<<axesDir);
          ASKAPCHECK((axesDir[0] == 0) && (axesDir[1] == 1),
                   "At present we support only images with first axes being the direction pixel axes, image "<<name<<
                   " has "<< axesDir);
          /// @todo double check that the rounding is correct for the ref pixel
          radec.setReferencePixel(radec.referencePixel()*double(*extraOversampleFactor));
          radec.setIncrement(radec.increment()/double(*extraOversampleFactor));
          imageCoords.replaceCoordinate(radec, whichDir);
      }
      else {
          imagePixels.reference(ip.valueF(name));
      }

      ASKAPDEBUGASSERT(imagePixels.ndim()!=0);
      ASKAPLOG_DEBUG_STR(logger, "Data of "<<name<<" parameter peak at "<<casacore::max(imagePixels));

      imageHandler().create(imagename, imagePixels.shape(), imageCoords);
      imageHandler().write(imagename, imagePixels);

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
      imageHandler().setMetadataKeywords(imagename, keywords);
      imageHandler().addHistory(imagename, historyLines);
    }

    /// @brief obtain image handler
    /// @details For some operations it may be necessary to access the (global) instance of the
    /// image handler. This method allows that. An exception is thrown if no image handler has
    /// been previously set up.
    /// @return a reference to image handler
    accessors::IImageAccess<casacore::Float>& SynthesisParamsHelper::imageHandler()
    {
      ASKAPCHECK(theirImageAccessor, "setUpImageHandler has to be called before any read/write operation");
      return *theirImageAccessor;
    }

    /// @brief configure default frequency frame
    /// @details All code workes in a single frequency frame (convertions are done, if
    /// necessary when the data are read (using conversion mechanism provided by the accessor).
    /// A call to this method sets up new default.
    /// @param[in] frame reference frame to use for all created images
    void SynthesisParamsHelper::setDefaultFreqFrame(const casacore::MFrequency::Ref &frame)
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

    /// @brief setup image handler
    /// @details This method uses the factory to setup a helper class handling the
    /// operations with images (default is casa). It is necessary to call this method
    /// at least once before any read or write operation can happen.
    /// @param[in] parset a parset file containing parameters describing which image handler to use
    /// @param[in] comms, MPI communicator
    /// @note The key parameter describing the image handler is "imagetype". By default, the
    /// casa image handler is created (however, a call to this method is still required)
    void SynthesisParamsHelper::setUpImageHandler(const LOFAR::ParameterSet &parset,
                                                  askap::askapparallel::AskapParallel &comms)
    {
      theirImageAccessor = accessors::imageAccessFactory(parset, comms);
    }


    void SynthesisParamsHelper::loadImageParameter(askap::scimath::Params& ip, const string& name,
						 const string& imagename, const boost::optional<float> extraOversampleFactor)
    {
      ASKAPTRACE("SynthesisParamsHelper::loadImageParameter");
      casacore::Array<float> imagePixels = imageHandler().read(imagename);
      casacore::CoordinateSystem imageCoords = imageHandler().coordSys(imagename);
      loadImageParameter(ip, name, imagePixels, imageCoords, extraOversampleFactor);
    }

    void SynthesisParamsHelper::loadImageParameter(askap::scimath::Params& ip, const string& name,
		casacore::Array<float>& imagePixels, const casacore::CoordinateSystem& imageCoords,
        const boost::optional<float> extraOversampleFactor, int channel)
    {
      /// Fill in the axes information
      Axes axes;
      /// First do the direction
      int whichDir=imageCoords.findCoordinate(Coordinate::DIRECTION);
      ASKAPCHECK(whichDir>-1, "No direction coordinate present in the supplied image coordinates");
      casacore::DirectionCoordinate radec(imageCoords.directionCoordinate(whichDir));
      casacore::Vector<casacore::Int> axesDir = imageCoords.pixelAxes(whichDir);
      ASKAPCHECK(axesDir.nelements() == 2, "Direction axis "<<whichDir<<
                 " is expected to correspond to just two pixel axes, you have "<<axesDir);
      ASKAPCHECK((axesDir[0] == 0) && (axesDir[1] == 1),
               "At present we support only images with first axes being the direction pixel axes, image "<<name<<
               " has "<< axesDir);

      casacore::Vector<casacore::String> units(2);
      units.set("rad");
      radec.setWorldAxisUnits(units);
      axes.addDirectionAxis(radec);

      int whichStokes = imageCoords.findCoordinate(Coordinate::STOKES);
      int nPol = 1;
      if (whichStokes<0) {
          const casacore::Vector<casacore::Stokes::StokesTypes> dummyStokes(1,casacore::Stokes::I);
          axes.addStokesAxis(dummyStokes);
      } else {
          casacore::StokesCoordinate sc(imageCoords.stokesCoordinate(whichStokes));
          const casacore::Vector<casacore::Int> stokesAsInt = sc.stokes();
          casacore::Vector<casacore::Stokes::StokesTypes> stokes(stokesAsInt.nelements());
          for (casacore::uInt pol=0; pol<stokes.nelements(); ++pol) {
               stokes[pol] = casacore::Stokes::StokesTypes(stokesAsInt[pol]);
          }
          axes.addStokesAxis(stokes);

          casacore::Vector<casacore::Int> axesStokes = imageCoords.pixelAxes(whichStokes);
          ASKAPCHECK(axesStokes.nelements() == 1, "Stokes axis "<<whichStokes<<
                 " is expected to correspond to just one pixel axes, you have "<<axesStokes);
          ASKAPASSERT(casacore::uInt(axesStokes[0])<imagePixels.shape().nelements());
          nPol = imagePixels.shape()(axesStokes[0]);
          ASKAPASSERT(uInt(nPol) == stokesAsInt.nelements());
      }

      int whichSpectral=imageCoords.findCoordinate(Coordinate::SPECTRAL);
      ASKAPCHECK(whichSpectral>-1, "No spectral coordinate present in model");
      casacore::Vector<casacore::Int> axesSpectral = imageCoords.pixelAxes(whichSpectral);
      ASKAPCHECK(axesSpectral.nelements() == 1, "Spectral axis "<<whichSpectral<<
                 " is expected to correspond to just one pixel axes, you have "<<axesSpectral);
      ASKAPASSERT(casacore::uInt(axesSpectral[0])<imagePixels.shape().nelements());
      const int nChan = (channel < 0 ? imagePixels.shape()(axesSpectral[0]) : 1);
      casacore::SpectralCoordinate freq(imageCoords.spectralCoordinate(whichSpectral));
      double startFreq, endFreq;
      if (channel < 0) {
        freq.toWorld(startFreq, 0.0);
        freq.toWorld(endFreq, double(nChan-1));
      } else {
        freq.toWorld(startFreq, double(channel));
        freq.toWorld(endFreq, double(channel));
      }
      axes.add("FREQUENCY", startFreq, endFreq);

      casacore::IPosition targetShape(4, imagePixels.shape()(0), imagePixels.shape()(1), nPol, nChan);
      ASKAPDEBUGASSERT(targetShape.product() == imagePixels.shape().product());

      if (extraOversampleFactor) {
          ASKAPCHECK(*extraOversampleFactor > 1.,"Oversampling factor should be > 1, not "<<*extraOversampleFactor);

          // need to add two params: the loaded full-res image and the downsampled image

          // first save the full-res image
          std::string fullresname = name;
          const size_t index = fullresname.find("image");
          ASKAPCHECK(index == 0, "Trying to swap to full-resolution param name but something is wrong");
          fullresname.replace(index,5,"fullres");

          ASKAPLOG_INFO_STR(logger, "About to add new image parameters with name "<<fullresname<<
                      " reshaped to "<<targetShape<<" from original image shape "<<imagePixels.shape());
          if (!ip.has(fullresname)) {
            ip.add(fullresname, imagePixels.reform(targetShape), axes);
          } else {
            ip.update(fullresname, imagePixels.reform(targetShape));
          }

          // now downsample the input sky model to the working resolution
          /// @todo should this be done at higher precision ifndef ASKAP_FLOAT_IMAGE_PARAMS?
          downsample(imagePixels,*extraOversampleFactor);
          // update the target shape for the new resolution
          targetShape(0) = imagePixels.shape()(0);
          targetShape(1) = imagePixels.shape()(1);
          // update the coordinate system for the new resolution
          /// @todo double check that the rounding is correct for the ref pixel
          radec.setReferencePixel(radec.referencePixel()/double(*extraOversampleFactor));
          radec.setIncrement(radec.increment()*double(*extraOversampleFactor));
          axes.addDirectionAxis(radec);

          ASKAPLOG_INFO_STR(logger, "Also adding downsampled image parameters with name "<<name<<
                      " reshaped to "<<targetShape<<" from original image shape "<<imagePixels.shape());
          if (!ip.has(name)) {
            ip.add(name, imagePixels.reform(targetShape), axes);
          } else {
            ip.update(name, imagePixels.reform(targetShape));
          }

      }
      else {
          ASKAPLOG_INFO_STR(logger, "About to add new image parameter with name "<<name<<
                      " reshaped to "<<targetShape<<" from original image shape "<<imagePixels.shape());
          if (!ip.has(name)) {
            ip.add(name, imagePixels.reform(targetShape), axes);
          } else {
            ip.update(name, imagePixels.reform(targetShape));
          }
      }
      ASKAPLOG_INFO_STR(logger, "Spectral axis will have startFreq="<<startFreq<<" Hz, endFreq="<<endFreq<<
                        "Hz, nChan="<<nChan);


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

    /// @brief zero-pad in the Fourier domain to increase resolution before cleaning
    /// @param[in] image the array to oversample
    /// @param[in] osfactor extra oversampling factor
    /// @param[in] norm bool true if we want to renormalise the array (e.g., keep psf peak at 1)
    /// @return oversampled image
    /// @todo move osfactor to itsOsFactor to enforce consistency between oversample() & downsample()?
    /// @todo use scimath::PaddingUtils::fftPad? Works with imtype rather than float so template there or here?
    void SynthesisParamsHelper::oversample(casacore::Array<float> &image, const float osfactor, const bool norm)
    {
        ASKAPCHECK(osfactor >= 1.0,
            "Oversampling factor in the solver is supposed to be greater than or equal to 1.0, you have "<<osfactor);
        if (osfactor == 1.) return;
        const uint ndim = image.ndim();
        ASKAPASSERT(ndim>1);
        const IPosition matShape = image.shape().nonDegenerate();
        const IPosition osShape = scimath::PaddingUtils::paddedShape(matShape,osfactor);
        // Limit number of fft threads to 8 (more is slower for our fft sizes)
        scimath::FFT2DWrapper<casacore::Complex> fft2d(true,8);

        // renormalise based on the imminent padding
        if (norm) {
            image *= static_cast<float>(osfactor*osfactor);
        }

        // Set up scratch array for changing resolution
        casacore::Matrix<casacore::Complex> Agrid(matShape,casacore::Complex(0.f));

        // copy image into complex scratch space & free memory
        casacore::setReal(Agrid, image.nonDegenerate());
        image.resize();

        // fft to uv
        fft2d(Agrid, true);

        // "extract" original Fourier grid into the central portion of the new Fourier grid
        casacore::Matrix<casacore::Complex> AgridOS(osShape,casacore::Complex(0.f));
        casacore::Array<casacore::Complex> subGrid = scimath::PaddingUtils::extract(AgridOS,osfactor);
        subGrid = Agrid;
        Agrid.resize();

        // ifft back to image and return the real part
        fft2d(AgridOS,false);

        image.reference(real(AgridOS).addDegenerate(ndim-2));
    }

    /// @brief remove Fourier zero-padding region to re-establish original resolution after cleaning
    /// @brief zero-pad in the Fourier domain to increase resolution before cleaning
    /// @param[in] image the array to oversample
    /// @param[in] osfactor extra oversampling factor
    /// @param[in] norm bool true if we want to renormalise the array (e.g., keep psf peak at 1)
    /// @return downsampled image
    /// @todo move osfactor to itsOsFactor to enforce consistency between oversample() & downsample()?
    /// @todo use scimath::PaddingUtils::fftPad? Downsampling may not be supported at this stage
    void SynthesisParamsHelper::downsample(casacore::Array<float> &image, const float osfactor, const bool norm)
    {
        ASKAPCHECK(osfactor >= 1.0,
            "Oversampling factor in the solver is supposed to be greater than or equal to 1.0, you have "<<osfactor);
        if (osfactor == 1.) return;

        // Limit number of fft threads to 8 (more is slower for our fft sizes)
        scimath::FFT2DWrapper<casacore::Complex> fft2d(true,8);
        const uint ndim = image.ndim();
        ASKAPASSERT(ndim>1);

        // Set up scratch arrays for changing resolution
        casacore::Matrix<casacore::Complex> AgridOS(image.shape().nonDegenerate(),casacore::Complex(0.f));

        // copy image into a complex scratch space & free memory
        casacore::setReal(AgridOS, image.nonDegenerate());
        image.resize();

        // fft to uv. This is what we need to degrid from, so no need to renormalise
        fft2d(AgridOS, true);

        // extract the central portion of the Fourier grid
        // We want Agrid to be contiguous, so force a copy by using separate assignment
        casacore::Matrix<casacore::Complex> Agrid;
        Agrid = scimath::PaddingUtils::extract(AgridOS,osfactor);
        // free memory
        AgridOS.resize();

        if (norm) {
            // renormalise based on unpadding
            Agrid /= static_cast<float>(osfactor*osfactor);
        }

        // ifft back to image and return the real part
        fft2d(Agrid, false);
        image.reference(real(Agrid).addDegenerate(ndim-2));
    }

    /// @brief determine sampling to use for nyquistgridding
    /// @param[in] advice VisMetaDataStats object used to get max U, V, W
    /// @param[in/out] parset to read and modify (setting extraoversampling and changing cellsize and shape/subshape)
    // @todo should probably throw an exception if individual cell or image sizes are given
    void SynthesisParamsHelper::setNyquistSampling(const VisMetaDataStats& advice, LOFAR::ParameterSet& parset) {
        // add Nyquist gridding parameters if needed. Wait until after doing others requiring VisMetaDataStats.
        ASKAPCHECK(!parset.isDefined("Images.extraoversampling"), "Images.extraoversampling cannot be set by user");
        if (parset.getBool("Images.nyquistgridding",false) || parset.isDefined("Images.griddingcellsize")) {

            ASKAPCHECK(parset.isDefined("Images.cellsize"),
            "The global image cellsize is currently required with Nyquist gridding");
            const std::vector<double> cellSize = convertQuantity(parset.getStringVector("Images.cellsize"),"arcsec");
            ASKAPCHECK(cellSize.size() == 2, "nyquistgridding requires a cellsize vector of length 2");
            ASKAPCHECK(cellSize[0] == cellSize[1], "nyquistgridding only set up for square pixels");

            const std::vector<std::string> names = parset.getStringVector("Images.Names");
            const int n = names.size();
            ASKAPCHECK(n > 0,"At least one image name needs to be specified");
            Vector<int> sizes(n);
            int minSize;
            if (parset.isDefined("Images.shape")) {
                auto shape = parset.getInt32Vector("Images.shape");
                ASKAPCHECK(shape.size()==2,"The image shape needs two dimensions");
                ASKAPCHECK(shape[0]==shape[1],"The image shape needs to be square for Nyqyuist gridding");
                sizes = shape[0];
                minSize = sizes[0];
            }
            // check if we have image specific image shapes - use smallest one to set gridding cellsize
            // and make sure image shapes only differ by integer factors from smallest image
            // find smallest image and make sure size is sensible
            for (int i=0; i<n; i++) {
                if (parset.isDefined("Images."+names[i]+".shape")) {
                    auto shape = parset.getInt32Vector("Images."+names[i]+".shape");
                    ASKAPCHECK(shape.size()==2,"The image shape needs two dimensions");
                    ASKAPCHECK(shape[0]==shape[1],"The image shape needs to be square for Nyqyuist gridding");
                    sizes[i] = parset.getInt32Vector("Images."+names[i]+".shape")[0];
                }
            }
            minSize = scimath::goodFFTSize(*std::min_element(sizes.begin(),sizes.end()));

            // now make other images multiples of this
            for (int i=0; i<n; i++) {
                const int factor = static_cast<int>(ceil(static_cast<double>(sizes[i])/minSize));
                if (sizes[i] % minSize != 0) {
                    // adjust each size to be a multiple of minSize
                    sizes[i] = minSize * factor;
                    ASKAPLOG_INFO_STR(logger, "Adjusting image size for "<<names[i]<<
                    " to be a multiple of smallest size ("<<minSize<<") - factor: "<<factor);
                }
            }

            // need to make sure that extraOsFactor results in an integer number of pixels,
            // which could get complicated for rectangular grids and pixels.

            std::vector<double> gCellSize(2);
            if (parset.isDefined("Images.griddingcellsize")) {
                const std::vector<string> gParam = parset.getStringVector("Images.griddingcellsize");
                ASKAPCHECK(gParam.size() == 2, "nyquistgridding requires a griddingcellsize vector of length 2");
                gCellSize = SynthesisParamsHelper::convertQuantity(gParam,"arcsec");
            }
            else {
                const double uv_max = casacore::max(advice.maxU(), advice.maxV());
                const double fov = cellSize[0] * sizes[0] * casacore::C::arcsec;
                const double wk_max = 6/fov + advice.maxW()*fov;
                ASKAPASSERT(uv_max > 0);
                // calculate the resolution in arcsec corresponding to the smallest grid that will fit all the data
                //  - the reciprical of twice the longest baseline plus the largest w support
                //  - doubled wk_max because it seemed too small
                gCellSize[0] = 0.5 / (uv_max + 2*wk_max) / casacore::C::arcsec;
                gCellSize[1] = gCellSize[0];
            }
            ASKAPCHECK(gCellSize[0]==gCellSize[1], "nyquistgridding only set up for square pixels");
            ASKAPCHECK(gCellSize[0]>=cellSize[0], "griddingcellsize must not be less than cellsize");

            // nominal ratio between gridding resolution and cleaning resolution
            ASKAPDEBUGASSERT(cellSize[0] > 0);
            double extraOsFactor = gCellSize[0]/cellSize[0];
            // now tweak the ratio to result in an integer number of pixels and reset the gridding cell size
            ASKAPDEBUGASSERT(extraOsFactor >= 1);
            int nPix = static_cast<int>(ceil(minSize/extraOsFactor));
            // ensure we only use factors of 2, 3, 5 and 7 to avoid slow FFT issues AXA-2430
            nPix = scimath::goodFFTSize(nPix);

            int nSubPix = 0;
            if (parset.isDefined("Images.subshape")) {
                ASKAPCHECK(sizes.size()==1,"Cannot specify subshape if you have multiple images");
                const std::vector<int> subSize = parset.getInt32Vector("Images.subshape");
                ASKAPCHECK(subSize.size() == 2, "nyquistgridding requires a subshape vector of length 2");
                ASKAPCHECK(subSize[0] > 1, "subshape image size too small: "<<subSize[0]);
                const int factor = static_cast<int>(ceil(static_cast<double>(sizes[0])/subSize[0]));
                if (sizes[0] % subSize[0] != 0) {
                    // adjust imSize to be a multiple of subSize
                    sizes[0] = subSize[0] * factor;
                    ASKAPLOG_INFO_STR(logger, "Adjusting image shape to be a multiple of subshape - factor: "<<factor);
                }
                nSubPix = static_cast<int>(ceil(subSize[0]/extraOsFactor));
                ASKAPDEBUGASSERT(nSubPix > 1);
                // ensure we only use factors of 2, 3, 5 and 7 to avoid slow FFT issues AXA-2430
                nSubPix = scimath::goodFFTSize(nSubPix);
                // reset the extra multiplicative factor
                extraOsFactor = static_cast<double>(subSize[0]) / nSubPix;
                nPix = nSubPix * factor;
            } else {
                // reset the extra multiplicative factor
                extraOsFactor = static_cast<double>(minSize) / nPix;
            }

            ASKAPDEBUGASSERT(extraOsFactor >= 1);

            gCellSize[0] = cellSize[0] * extraOsFactor;
            gCellSize[1] = gCellSize[0];

            // Only set or reset these parameters if they are needed
            //  - could require a minimum increase factor (20%, 50%, 100%, etc.)
            if (extraOsFactor > 1.) {
                ASKAPLOG_INFO_STR(logger, "  Adding new parameter extraoversampling = "<<extraOsFactor);
                parset.add("Images.extraoversampling", utility::toString(extraOsFactor));
                ASKAPLOG_INFO_STR(logger, "  Changing cellsize from "<<parset.getStringVector("Images.cellsize")<<
                " to "<<"["<<gCellSize[0]<<"arcsec,"<<gCellSize[1]<<"arcsec]");
                parset.replace("Images.cellsize", "["+utility::toString(gCellSize[0])+"arcsec,"+utility::toString(gCellSize[1])+"arcsec]");
                const std::string imShape0("Images.shape");
                bool first = true;
                for (int i=0; i<n; i++) {
                    const int factor = static_cast<int>(ceil(static_cast<double>(sizes[i])/minSize));
                    sizes[i] = nPix * factor;
                    const std::string imShape("Images."+names[i]+".shape");
                    if (parset.isDefined(imShape)) {
                        ASKAPLOG_INFO_STR(logger, "  Changing shape for "<<names[i]<<" from "<<parset.getInt32Vector(imShape)<<
                        " to "<<"["<<sizes[i]<<","<<sizes[i]<<"]");
                        parset.replace(imShape, "["+utility::toString(sizes[i])+","+utility::toString(sizes[i])+"]");
                    } else {
                        ASKAPLOG_INFO_STR(logger, "  Changing shape for "<<names[i]<<" from "<<parset.getInt32Vector(imShape0)<<
                        " to "<<"["<<sizes[i]<<","<<sizes[i]<<"]");
                        if (first) {
                            parset.replace(imShape0, "["+utility::toString(sizes[i])+","+utility::toString(sizes[i])+"]");
                            first = false;
                        }
                    }
                }
                if (nSubPix > 0) {
                    ASKAPLOG_INFO_STR(logger, "  Changing subshape from "<<parset.getInt32Vector("Images.subshape")<<
                    " to "<<"["<<nSubPix<<","<<nSubPix<<"]");
                    parset.replace("Images.subshape", "["+utility::toString(nSubPix)+","+utility::toString(nSubPix)+"]");
                }
            }
        }
    }

    void SynthesisParamsHelper::adjustCellsize(Array<float> & pixels, double inputInc, double outputInc, 
        uInt inSize, uInt outSize)
    {
        const double tol = 1e-5;
        ASKAPDEBUGASSERT(inputInc != 0);
        const double fac = outputInc / inputInc;
        ASKAPCHECK(fac > 1-tol,"Input cellsize should be smaller than or equal to output cellsize");
        if (abs(fac-1) > tol) {
            // try to downsample - need to use exact cellsize ratio but
            // need to make sure image sizes are good FFT sizes for downsampling
            // We try padding or cropping the input image
            bool found = false;
            uInt n1 = inSize;
            for (int k = 0; k < inSize/2 && !found ; k++) {
                // first try padding
                n1 = scimath::goodFFTSize(inSize + 2 * k);
                const uInt n2 = scimath::goodFFTSize(static_cast<int>(lround(n1/fac)));
                found = abs(double(n1)/n2 - fac) < tol;
                ASKAPLOG_DEBUG_STR(logger,"tried n1="<<n1<<", n2="<<n2<<", fac="<<fac<<", n1/n2="<<double(n1)/n2<<", match="<<found);
                // try cropping if the input image covers a larger area than the output
                if (!found && k > 0) {
                    n1 = scimath::goodFFTSize(inSize - 2 * k);
                    if (inputInc*n1 > outputInc*outSize) {
                        const uInt n2 = scimath::goodFFTSize(static_cast<int>(lround(n1/fac)));
                        found = abs(double(n1)/n2 - fac) < tol;
                        ASKAPLOG_DEBUG_STR(logger,"tried n1="<<n1<<", n2="<<n2<<", fac="<<fac<<", n1/n2="<<double(n1)/n2<<", match="<<found);
                    } 
                }
            }
            ASKAPCHECK(found,"No valid pixel ratio found for given cell sizes "<<inputInc<<" and "<<outputInc);
            // pad or crop the pixels array
            if (n1 > inSize) {
                // pad
                Array<float> imagePixels(IPosition(4,n1,n1,1,1),0.0f);
                // make reference to central area
                IPosition blc(4, n1/2 - inSize/2, n1/2 - inSize/2, 0, 0);
                IPosition trc(4, n1/2 + inSize/2 - 1, n1/2 + inSize/2 - 1, 0, 0);
                Array<float> subImage = imagePixels(blc, trc);
                // copy in the pixels
                subImage = pixels;
                SynthesisParamsHelper::downsample(imagePixels,fac,false);
                // replace original array with the downsampled array
                pixels.reference(imagePixels);
            } else if (n1 < inSize) {
                // crop
                IPosition blc(4, inSize/2 - n1/2, inSize/2 - n1/2, 0, 0);
                IPosition trc(4, inSize/2 + n1/2 - 1, inSize/2 + n1/2 - 1, 0, 0);
                Array<float> imagePixels;
                // make copy of central area
                imagePixels = pixels(blc, trc);
                SynthesisParamsHelper::downsample(imagePixels,fac,false);
                // replace original array with the downsampled array
                pixels.reference(imagePixels);
            } else {
                SynthesisParamsHelper::downsample(pixels,fac,true);
            }
        }
    }


    boost::shared_ptr<casacore::TempImage<float> >
    SynthesisParamsHelper::tempImage(const askap::scimath::Params& ip,
				     const string& name)
    {
      const casacore::Array<float> imagePixels(ip.valueF(name));

      casacore::CoordinateSystem imageCoords(coordinateSystem(ip, name));

      boost::shared_ptr<casacore::TempImage<float> >
      im(new casacore::TempImage<float> (TiledShape(imagePixels.shape()),
				       imageCoords));

      im->setUnits("Jy/pixel");

      casacore::ArrayLattice<float> latImagePixels(imagePixels);
      im->copyData(latImagePixels);
      return im;
    }

    casacore::CoordinateSystem
    SynthesisParamsHelper::coordinateSystem(const askap::scimath::Params& ip,
					    const string& name)
    {
      const Axes axes(ip.axes(name));

      casacore::DirectionCoordinate radec(directionCoordinate(ip, name));

      casacore::CoordinateSystem imageCoords;
      imageCoords.addCoordinate(radec);


      const casacore::IPosition shape = ip.shape(name);
      const int nchan = shape.nelements() >= 4 ? shape(3) : 1;
      const double restfreq = 0.0;
      const double crpix = double(nchan-1)/2.;
      const double crval = (axes.start("FREQUENCY")+axes.end("FREQUENCY"))/2.0;
      // we can't estimate increment if there is only one channel and start=stop
      const double cdelt = nchan>1 ? (axes.end("FREQUENCY")-axes.start("FREQUENCY"))/double(nchan-1) : 1.;
      const casacore::SpectralCoordinate freq(casacore::MFrequency::castType(theirFreqFrame.getType()), crval, cdelt, crpix, restfreq);
      imageCoords.addCoordinate(freq);

      // default is a dummy stokes coordinate with only stokes I present
      casacore::Vector<int> iquv(1);
      iquv(0) = Stokes::I;
      if (axes.has("STOKES")) {
          casacore::Vector<casacore::Stokes::StokesTypes> stokes = axes.stokesAxis();
          ASKAPDEBUGASSERT(stokes.nelements()>=1);
          iquv.resize(stokes.nelements());
          for (size_t pol=0; pol<stokes.nelements(); ++pol) {
               iquv[pol] = int(stokes[pol]);
          }
      }

      casacore::StokesCoordinate stokes(iquv);
      imageCoords.addCoordinate(stokes);

      return imageCoords;
    }

    casacore::DirectionCoordinate
    SynthesisParamsHelper::directionCoordinate(const askap::scimath::Params& ip,
					       const string& name)
    {
      const Axes axes(ip.axes(name));
      ASKAPCHECK(!axes.has("RA-TANGENT") && !axes.has("DEC-TANGENT"),
          "Detected an obsolete way to specify tangent point!");

      ASKAPCHECK(axes.hasDirection(), "Direction coordinate is missing. axes["<<name<<"]:"<<axes);
      return axes.directionAxis();

    }

    void SynthesisParamsHelper::copyImageParameter(askap::scimath::Params& sourceParam,  askap::scimath::Params& sinkParam, const string& name) {

      ASKAPDEBUGTRACE("SynthesisParamsHelper::copyImageParameter");
      ASKAPDEBUGASSERT(sourceParam.has(name));

      // Stephen Ord 2018.

      // the purpose of this is to produce a local sky model from a global one.
      // in general we will be replacing the model in the local Params - with
      // a chunk from the global model.

      // we need to get a reference to the input and output Images
      // this tempImage method is a copy of the param.value. It is not the value itself.
      // we need to call update on this for the operation to be complete

      boost::shared_ptr<casacore::TempImage<float> > inRef = tempImage(sourceParam,name); // inRef

      // regrid the input plane into the output one.

      boost::shared_ptr<casacore::TempImage<float> > outRef = tempImage(sinkParam,name); // outRef

      // regridder
      askap::imagemath::LinmosImageRegrid<float> regridder;

      const askap::imagemath::Interpolate2D::Method method = askap::imagemath::Interpolate2D::CUBIC;
      const casacore::uInt decimate = 1;
      casacore::IPosition itsAxes = IPosition::makeAxisPath(outRef->shape().nonDegenerate().nelements());

      // the regridding is taking the input array and resampling it into the coorrdinate system of the
      // output array. I dont think it combines - I think it just replaces. The output can be empty.

      regridder.regrid(*outRef, method,
                          itsAxes, *inRef, false, decimate,true,true);

      // now we have to insert the regridded image into the pixels of the output ...
      update(sinkParam,name,*outRef);

      // so this should pretty much do the job. I will put together a unit test to check that this is
      // doing what I think it is doing.



    }
    void SynthesisParamsHelper::update(askap::scimath::Params& ip, const string& name,
				       const casacore::ImageInterface<float>& im)
    {
      ASKAPDEBUGTRACE("SynthesisParamsHelper::update");
      /// This next copy should be a reference unless it is too big
      casacore::Array<float> floatImagePixels(im.shape());
      casacore::ArrayLattice<float> latImagePixels(floatImagePixels);
      latImagePixels.copyData(im);

      ip.update(name, floatImagePixels);
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
    casacore::Slicer SynthesisParamsHelper::facetSlicer(const askap::scimath::Params& ip,
           const std::string &name, const casacore::DirectionCoordinate &dc)
    {
       ASKAPDEBUGASSERT(ip.has(name));
       const casacore::DirectionCoordinate inDC = directionCoordinate(ip,name);
       const casacore::IPosition inShape = ip.shape(name);
       ASKAPDEBUGASSERT(inShape.nelements() >= 2);
       casacore::IPosition blc(inShape.nelements(),0);
       casacore::IPosition trc(inShape);
       for (casacore::uInt dim=0; dim<inShape.nelements(); ++dim) {
            --trc(dim);
       }
       // currently blc,trc describe the whole input image; convert coordinates
       casacore::Vector<casacore::Double> pix(2);

       // first process BLC
       pix[0] = casacore::Double(blc[0]);
       pix[1] = casacore::Double(blc[1]);
       casacore::MDirection tempDir;
       casacore::Bool success = inDC.toWorld(tempDir, pix);
       ASKAPCHECK(success, "Pixel to world coordinate conversion failed for input BLC: "<<inDC.errorMessage());
       success = dc.toPixel(pix,tempDir);
       ASKAPCHECK(success, "World to pixel coordinate conversion failed for output BLC: "<<dc.errorMessage());
       // converting to Int without rounding appears to be changing the image size.
       //blc[0] = casacore::Int(pix[0]);
       //blc[1] = casacore::Int(pix[1]);
       blc[0] = casacore::Int(round(pix[0]));
       blc[1] = casacore::Int(round(pix[1]));

       // now process TRC
       pix[0] = casacore::Double(trc[0]);
       pix[1] = casacore::Double(trc[1]);
       success = inDC.toWorld(tempDir, pix);
       ASKAPCHECK(success, "Pixel to world coordinate conversion failed for input TRC: "<<inDC.errorMessage());
       success = dc.toPixel(pix,tempDir);
       ASKAPCHECK(success, "World to pixel coordinate conversion failed for output TRC: "<<dc.errorMessage());
       // converting to Int without rounding appears to be changing the image size.
       //trc[0] = casacore::Int(pix[0]);
       //trc[1] = casacore::Int(pix[1]);
       trc[0] = casacore::Int(round(pix[0]));
       trc[1] = casacore::Int(round(pix[1]));

       return casacore::Slicer(blc,trc,casacore::Slicer::endIsLast);
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
       const casacore::DirectionCoordinate templateDC = directionCoordinate(ip,names[0]);
       // we could've generated initial slicer from the shape and avoid one unnecessary call to facetSlicer
       const casacore::Slicer tempSlicer = facetSlicer(ip,names[0],templateDC);
       casacore::IPosition tempBLC = tempSlicer.start();
       casacore::IPosition tempTRC = tempSlicer.end();
       ASKAPDEBUGASSERT(tempBLC.nelements() >= 2);
       ASKAPDEBUGASSERT(tempTRC.nelements() >= 2);
       for (size_t i = 1; i<names.size(); ++i) {
            const casacore::Slicer newSlicer = facetSlicer(ip,names[i],templateDC);
            const casacore::IPosition newBLC = newSlicer.start();
            const casacore::IPosition newTRC = newSlicer.end();
            ASKAPDEBUGASSERT(newBLC.nelements() >= 2);
            ASKAPDEBUGASSERT(newTRC.nelements() >= 2);
            for (casacore::uInt dim=0; dim<2; ++dim) {
                 if (newBLC(dim) < tempBLC(dim)) {
                     tempBLC(dim) = newBLC(dim);
                 }
                 if (newTRC(dim) > tempTRC(dim)) {
                     tempTRC(dim) = newTRC(dim);
                 }
            }
       }
       casacore::IPosition newShape = ip.shape(names[0]);
       ASKAPDEBUGASSERT(newShape.nelements() >= 2);
       newShape(0) = tempTRC(0) - tempBLC(0) + 1;
       newShape(1) = tempTRC(1) - tempBLC(1) + 1;
       ASKAPDEBUGASSERT(newShape(0) > 0);
       ASKAPDEBUGASSERT(newShape(1) > 0);
       casacore::Vector<casacore::Double> refPix = templateDC.referencePixel();
       refPix[0] -= casacore::Double(tempBLC(0) - tempSlicer.start()(0));
       refPix[1] -= casacore::Double(tempBLC(1) - tempSlicer.start()(1));
       casacore::DirectionCoordinate newDC(templateDC);
       newDC.setReferencePixel(refPix);

       askap::scimath::Axes newAxes(ip.axes(names[0]));
       newAxes.addDirectionAxis(newDC);

       casacore::Array<imtype> pixels(newShape);
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
    casacore::Vector<casacore::Quantum<double> > SynthesisParamsHelper::fitBeam(const askap::scimath::Params &ip,
                        const double cutoff, const int maxsupport, const std::string &name)
    {
       ASKAPTRACE("SynthesisParamsHelper::fitBeam");

       std::string psfName = name;
       if (name == "") {
           // we have to figure out the name
           psfName = findPSF(ip);
       }
       ASKAPCHECK(psfName != "", "Unable to find psf parameter to fit, params="<<ip);
       ASKAPLOG_INFO_STR(logger, "Fitting 2D Gaussian into PSF parameter "<<psfName);

       casacore::Array<imtype> psfArray = ip.valueT(psfName);
       return fitBeam(psfArray, ip.axes(psfName), cutoff, maxsupport);

    }

    casacore::Vector<casacore::Quantum<double> > SynthesisParamsHelper::fitBeam(casacore::Array<imtype> &psfArray,
       const scimath::Axes &axes, const double cutoff, const int maxsupport)
    {
        ASKAPCHECK(axes.hasDirection(), "Direction axes are missing from the PSF parameter, unable to convert pixels to angular units");
        const casacore::Vector<casacore::Double> increments = axes.directionAxis().increment();
        ASKAPCHECK(increments.nelements() == 2, "Expect just two elements for increments of the direction axis, you have "<<
                   increments);
        ASKAPCHECK(increments[1]>0, "Expect positive increment on the declination axis. increments="<<increments);
        ASKAPCHECK(fabs(fabs(increments[0])-fabs(increments[1]))<1e-6,
                   "Different cell sizes mean that the current beam fitting code would give a wrong position angle. increments="
                   <<increments);

       casacore::Vector<double> result = fitBeam(psfArray, cutoff, maxsupport);

       casacore::Vector<casacore::Quantum<double> > beam(3);
       beam[0] = casacore::Quantum<double>(fabs(increments[0])*result[0],"rad");
       beam[1] = casacore::Quantum<double>(fabs(increments[1])*result[1],"rad");
       // position angle in radians
       // fitBeam call above already does the pi/2 offset
       //double pa = increments[0]<0 ? result[2] - casacore::C::pi/2 : casacore::C::pi/2 - result[2];
       double pa = increments[0]<0 ? result[2] : -result[2];
       if (pa < -casacore::C::pi/2) {
           pa += casacore::C::pi;
       }
       beam[2] = casacore::Quantum<double>(pa,"rad");
       return beam;
    }

    casacore::Vector<double> SynthesisParamsHelper::fitBeam(casacore::Array<imtype> &psfArray,
        const double cutoff, const int maxsupport) {

       const casacore::IPosition shape = psfArray.shape();
       ASKAPCHECK(shape.nelements()>=2,"PSF image is supposed to be at least 2-dimensional, shape="<<psfArray.shape());
       ASKAPCHECK(cutoff>0 && cutoff<1,"beam cutoff level should be between 0 and 1");
       // make maxsupport a valid odd value
       int maxSupport = maxsupport;
       maxSupport = casacore::min(maxSupport,shape(0)-1);
       maxSupport = casacore::min(maxSupport,shape(1)-1);
       maxSupport = casacore::max(3,maxSupport);
       maxSupport = 2 * (maxSupport/2) + 1;

       if (shape.product() != shape[0]*shape[1]) {
           ASKAPLOG_WARN_STR(logger, "Multi-dimensional PSF is present (shape="<<shape<<
                             "), using the first 2D plane only to fit the beam");
       }
       ASKAPCHECK(shape[0] >= 3 && shape[1] >= 3, "Expect at least 3x3 pixel images, you have shape = "<<shape);

       casacore::Matrix<imtype> psfSlice = imagemath::MultiDimArrayPlaneIter::getFirstPlane(psfArray).nonDegenerate();

       // search for support to speed up beam fitting
       ASKAPLOG_INFO_STR(logger, "Searching for support with the relative cutoff of "<<cutoff<<" to speed fitting up");
       SupportSearcher ss(cutoff);
       ss.search(psfSlice);
       // search only looks in x and y direction, now extend the support to diagonals
       ss.extendedSupport(psfSlice);
       casacore::uInt support = ss.symmetricalSupport(psfSlice.shape());

       // limit support size to roughly 100 pixels per beam
       // some mostly flagged channels can have very large support, fit will take too long
       if (support > maxSupport) {
           // maybe we should just fail here and abort processing of the channel
           ASKAPLOG_DEBUG_STR(logger, "Reducing support size for beam fit from "<<support<<" to "<<maxSupport);
           support = maxSupport;
       }

       if (support < 3) {
           support = 3;
       } else {
          if (support % 2 == 0) {
              // if even, move up to next odd.
              ++support;
          }
       }
       ASKAPDEBUGASSERT(support % 2 == 1);

       ASKAPLOG_INFO_STR(logger, "Extracting support of "<<support<<" pixels for 2D gaussian fitting");
       const casa::IPosition newShape(2,support,support);
       for (int dim=0; dim<2; ++dim) {
            ASKAPCHECK(psfSlice.shape()[dim] >= int(support), "Support is greater than the original size, shape="<<
                       psfSlice.shape());
       }
       //
       #ifdef ASKAP_FLOAT_IMAGE_PARAMS
       // avoid making a reference copy since we rescale in next step
       casa::Array<float> floatPSFSlice;
       floatPSFSlice = scimath::PaddingUtils::centeredSubArray(psfSlice,newShape);
       #else
       casa::Array<float> floatPSFSlice(newShape);
       casa::convertArray<float, imtype>(floatPSFSlice, scimath::PaddingUtils::centeredSubArray(psfSlice,newShape));
       #endif
       // hack for debugging only
       //floatPSFSlice = imageHandler().read("tmp.img").nonDegenerate();
       //
       // normalise to 1 - technically unnecessary as we should have it already normalised
       const float maxPSF = casa::max(floatPSFSlice);
       ASKAPLOG_DEBUG_STR(logger," cutoff: "<<cutoff<<" maxPSF="<<maxPSF<<" newShape: "<<newShape);
       if (fabs(maxPSF-1.)>1e-6) {
            floatPSFSlice /= maxPSF;
       }
       //

       // actual fitting
       // the beam fitter fails at times - producing no or a bad solution
       // We've tried a few approaches:
       //  - retry with different support - often works, but still some failures
       //  - use setIncludeRange to exclude pixels below the cutoff - works sometimes, but worse at times
       //  - use estimate to set the initial guess - this seems the best solution so far
       casa::LogIO os;
       casa::Fit2D fitter(os);
       // do not set cutoff for small support, i.e. all pixels will be used for the fit
       // note, we can't get support less than 3 given the code above
       if (support > 3) {
           // for large support case ensure that we get at least 4 neighbouring pixels near the peak above the cutoff
           double newCutoff = cutoff;
           const int centrePixel = static_cast<int>(support) / 2;
           for (int pix = 0; pix < 4; ++pix) {
                // this gives +/- 1 pixel from centre on each axis and all combinations, assumes integer math.
                const IPosition cursor(2, centrePixel + (pix % 2) * 2 - 1, centrePixel + (pix / 2) * 2 - 1);
                const float curVal = floatPSFSlice(cursor);
                if (curVal < newCutoff) {
                    newCutoff = curVal;
                }
           }
           if (newCutoff < cutoff) {
               ASKAPLOG_DEBUG_STR(logger, "Reducing cutoff for fitter to ensure enough pixels are taken into account, new cutoff "<<newCutoff<<
                         ", it probably means the beam is bad");
           }
           fitter.setIncludeRange(newCutoff,1.0);
       }
       // Using estimate to set the initial guess seems to make the fit more robust
       casa::Vector<casa::Double> initialEstimate = fitter.estimate(casa::Fit2D::GAUSSIAN,floatPSFSlice);
       ASKAPLOG_DEBUG_STR(logger,"Initial beam fit: "<<initialEstimate);
       initialEstimate[0]=1.; // PSF peak is always 1
       initialEstimate[1]=newShape[0]/2; // centre
       initialEstimate[2]=newShape[1]/2; // centre
       casa::Vector<casa::Bool> parameterMask(6,casa::False);
       parameterMask[3] = casa::True; // fit maj
       parameterMask[4] = casa::True; // fit min
       parameterMask[5] = casa::True; // fit pa

       fitter.addModel(casa::Fit2D::GAUSSIAN,initialEstimate,parameterMask);
       const casa::Array<casa::Float> sigma(floatPSFSlice.shape(),1.f);
       const casa::Fit2D::ErrorTypes fitError = fitter.fit(floatPSFSlice,sigma);
       ASKAPCHECK(fitError == casa::Fit2D::OK, "Error fitting the beam. fitError="<<fitError<<
                  " message: "<<fitter.errorMessage());
       ASKAPCHECK(fitter.numberPoints() > 4, "The number of points available for fitting the beam ("<<fitter.numberPoints()<<
                  ") is too small. Either cell size or cutoff are too large");

       const casa::Vector<casa::Double> result = fitter.availableSolution();
       const casa::Vector<casa::Double> errors = fitter.availableErrors();

       ASKAPLOG_DEBUG_STR(logger, "Got fit result (in pixels) "<<result<<" and uncertainties "<< errors<<
                          " number of points used for the fit: "<<fitter.numberPoints());
       ASKAPCHECK(result.nelements() == 6, "Expect 6 parameters for 2D gaussian, result vector has "<<result.nelements());
       // Check if fit is reasonable
       for (casacore::uInt i = 0; i < result.nelements(); ++i) {
            ASKAPCHECK(!isnan(result[i]), "Beam fitting produces NaN in the result ("<<result<<"), most likely beam is not sampled properly, decrease the cell size");
            ASKAPCHECK(!isinf(result[i]), "Beam fitting produces infinity in the result ("<<result<<"), most likely beam is not sampled properly, decrease the cell size");
       }
       ASKAPCHECK(result[3] <= 2.0 * support && result[4] <= 2.0 * support,
                  "Error fitting the beam. Gaussian width >2x support");
       casacore::Vector<double> beam(3);
       beam[0] = result[3];
       beam[1] = result[4];
       // position angle on the pixel grid in radians
       double pa = result[5] - casacore::C::pi/2;
       if (pa < -casacore::C::pi/2) {
           pa += casacore::C::pi;
       }
       beam[2] = pa;
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
              //const std::string name="image"+*it;
              std::string name="image"+*it;
              if (params->isFree(name)) {
                  params->valueT(name).set(0.);
              }
              name="fullres"+*it;
              if (params->isFree(name)) {
                  params->valueT(name).set(0.);
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
    /// @param[in] compName name of the component
    /// @param[in] baseKey a prefix added to parset parameter names (default
    /// is "sources.", wich matches the current layout of the parset file)
    void SynthesisParamsHelper::copyComponent(const askap::scimath::Params::ShPtr &params,
           const LOFAR::ParameterSet &parset, const std::string &srcName,
           const std::string &compName, const std::string &baseKey)
    {
       ASKAPDEBUGASSERT(params);
       // check the special case of predefined calibrators
       const std::string calParamName = baseKey + compName + ".calibrator";
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
       parameterList["flux.spectral_index"] = false;
       parameterList["flux.ref_freq"] = false;

       // DDCALTAG COMPTAG -- keep a list of the sources and give each a unique ID
       const std::vector<std::string> sourceList = params->completions("sourceID.");
       if (std::find(sourceList.begin(), sourceList.end(), srcName) == sourceList.end()) {
           params->add("sourceID."+srcName, sourceList.size());
       }

       // DDCALTAG COMPTAG -- link this component to its source ID
       //ASKAPLOG_INFO_STR(logger, "DDCALTAG linking component "<< compName<<
       //   " to source "<<params->scalarValue("sourceID."+srcName));
       params->add("source."+compName, params->scalarValue("sourceID."+srcName));

       // now iterate through all parameters
       for (std::map<std::string, bool>::const_iterator ci = parameterList.begin();
            ci!=parameterList.end(); ++ci) {
            const std::string parName = baseKey+compName+"."+ci->first;
            if (parset.isDefined(parName)) {
                const double val = parset.getDouble(parName);
                params->add(ci->first+"."+compName, val);
            } else {
                if (ci->second) {
                    ASKAPTHROW(AskapError, "Parameter "<<parName<<
                           " is required to define the source "<<compName<<
                           ", baseKey="<<baseKey<<" or "<<baseKey + compName +".calibrator parameter should be present");
                }
            }
       }
    }
  }
}
