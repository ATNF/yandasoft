/// @file
/// @brief Specialised gridder just for uv-weight construction
/// @details We don't need everything from the gridder (i.e. the actual gridding, CF generation, etc) for
/// the weight construction. This is essentially a cut down version of the Box gridder trimmed specifically
/// so it can only construct weight.
///
/// @copyright (c) 2023 CSIRO
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

// own includes
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>

#include <askap/gridding/UVWeightGridder.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/AskapError.h>
#include <askap/scimath/utils/PaddingUtils.h>

ASKAP_LOGGER(logger, ".gridding.uvweightgridder");

namespace askap {

namespace synthesis {

/// @brief default constructor
/// @note this class constructed via the default constructor will be useless without the builder set (via setUVWeightBuilder call)
UVWeightGridder::UVWeightGridder() : itsPaddingFactor(1.f), itsUCellSize(0.), itsVCellSize(0.), itsMaxPointingSeparation(-1.), 
       itsFirstAccumulatedVis(false), itsDoBeamAndFieldSelection(true), itsSourceIndex(0u), itsCurrentField(0u),
       itsPointingTolerance(0.0001), itsOversample(1)
{}

/// @brief constructor setting the weight builder up front
/// @details Equivalent to the default constructor followed by a call to setUVWeightBuilder
/// @param[in] wtBuilder shared pointer to the weight builder to use
UVWeightGridder::UVWeightGridder(const boost::shared_ptr<IUVWeightBuilder> &wtBuilder) : itsPaddingFactor(1.f), itsUCellSize(0.), itsVCellSize(0.), 
       itsUVWeightBuilder(wtBuilder), itsMaxPointingSeparation(-1.),
       itsFirstAccumulatedVis(false), itsDoBeamAndFieldSelection(true), itsSourceIndex(0u), itsCurrentField(0u),
       itsPointingTolerance(0.0001), itsOversample(1)
{}


/// @brief Initialise the gridding and the associated builder class
/// @details This method is supposed to be called before gridding first data. For convenience parameters resemble those
/// the proper gridders from the IVisGridder class hierarchy are using. In particular, the shape parameter is 4-dimensional
/// (as used for the gridders) with uSize, vSize, nPol and nChan as opposed to the 3-dimensional shape used for weight grids
/// (uSize, vSize, nChan - i.e. it is assumed that we always have the same weight for all polarisation products).
/// @param axes axes specifications
/// @param shape desired shape of the weight grid, same as passed to the proper gridder for image creation, i.e. u, v, pol, chan
/// @note this method plays the role of initialiseGrid in the gridder hierarchy
void UVWeightGridder::initialise(const scimath::Axes& axes, const casacore::IPosition& shape)
{
   ASKAPDEBUGASSERT(shape.nelements()>=2);
   itsShape = scimath::PaddingUtils::paddedShape(shape,paddingFactor());
   // the following section has some code duplication with TableVisGridder::initialiseCellSize, perhaps some of it should be moved to Axes
   itsAxes = axes;
   ASKAPCHECK(itsAxes.hasDirection(), "Direction axis is missing. itsAxes:"<<itsAxes);
   const casacore::Vector<casacore::Double> increments = itsAxes.directionAxis().increment();
   ASKAPCHECK(increments.nelements() == 2, "Expect 2 elements in the increment vector, you have "<<increments);
   itsUCellSize = 1./(increments[0]*double(itsShape[0]));
   itsVCellSize = 1./(increments[1]*double(itsShape[1]));

   // now initialise the builder class
   ASKAPCHECK(itsUVWeightBuilder, "Weight builder class is supposed to be set before initialise call!");
   itsUVWeightBuilder->initialise(itsShape[0], itsShape[1], itsShape.nelements() > 3 ? itsShape[3] : 1u);

   // to setup accumulation of the first encountered field
   itsFirstAccumulatedVis = true;
   // I (MV) am not sure we need to reset the list of known pointings which is used to index fields. However, this only matters if we're
   // going to reuse the weight gridder and run it multiple times for different data.
   itsKnownPointings.clear();

   // frequency mapping (similar to TableVisGridder::initialiseFreqMapping)
   if (itsAxes.has("FREQUENCY") && itsShape.nelements()>=4) {
      itsFreqMapper.setupImage(itsAxes, itsShape[3]);
   } else {
      ASKAPLOG_DEBUG_STR(logger, "Forced to use single spectral plane weight gridding (either "
               "FREQUENCY axis or the number of channels are missing");
      itsFreqMapper.setupSinglePlaneGridding();
   }

   ASKAPLOG_DEBUG_STR(logger, "UV Weight building is initialised assuming the tangent centre is "<<
        printDirection(getTangentPoint())<<" and the image centre "<<
        printDirection(getImageCentre()));
}

/// @brief process the visibility data.
/// @param acc const data accessor to work with
/// @note this method plays the role of 'generic' or 'grid' methods in the gridder hierarchy. I (MV) not sure at this stage whether
/// we need some selection methods to control what actually contributes to weights or should use the accessor selector instead 
/// (as this would be a separate iteration over the data anyway). The method is 'const' because the actual accumulation is done
/// by the builder and this class is unchanged except for various caches (like frequency mapper)
void UVWeightGridder::accumulate(const accessors::IConstDataAccessor& acc) const
{
   // it may be worth thinking about the mode where we don't bother figuring out field index but rather use 0
   // (although this index can always be ignored anyway by the builder class). This behaviour would mimic what we have in the
   // non-mosaicing gridders.
   indexField(acc);

   ASKAPCHECK(itsUVWeightBuilder, "weight builder is supposed to be set before calling to accumulate");

   const casacore::MVDirection imageCentre = getImageCentre();
   const casacore::MVDirection tangentPoint = getTangentPoint();
   ASKAPDEBUGASSERT(itsShape.nelements() >= 2);

   // the following code is borrowed from gridder, but OpenMP sections are removed (less benefits for weight gridding as the
   // effective "CF" is small and builder interface as it is doesn't support multithreaded weight addition)
   const casacore::Vector<casacore::RigidVector<double, 3> > &outUVW = acc.rotatedUVW(tangentPoint);

   const casacore::uInt nSamples = acc.nRow();
   const casacore::uInt nChan = acc.nChannel();
   const casacore::uInt nPol = acc.nPol();
   ASKAPCHECK(nPol > 0, "Accessor passed to the accumulate method should have at least one polarisation product");

   const casacore::Vector<casacore::Double>& frequencyList = acc.frequency();
   itsFreqMapper.setupMapping(frequencyList);

   // now loop over samples and add them to the grid
   ASKAPDEBUGASSERT(casa::uInt(nChan) <= frequencyList.nelements());
   const casa::Cube<casa::Bool>& flagCube = acc.flag();
   const casa::Cube<casa::Complex>& noiseCube = acc.noise();

   UVWeight uvWeightRW;
   for (casacore::uInt i=0; i<nSamples; ++i) {
        if (itsMaxPointingSeparation > 0.) {
           // need to reject samples, if too far from the image centre
           const casacore::MVDirection thisPointing  = acc.pointingDir1()(i);
           if (imageCentre.separation(thisPointing) > itsMaxPointingSeparation) {
               continue;
           }
        }
        if (itsFirstAccumulatedVis) {
            if (itsDoBeamAndFieldSelection) {
                itsSelectedBeam = acc.feed1()(i);
                itsSelectedPointing = acc.dishPointing1()(i);
                ASKAPLOG_DEBUG_STR(logger, "Using the data for beam "<<itsSelectedBeam<<
                  " and field at "<<printDirection(itsSelectedPointing)<<" for uv weight construction");
            } else {
               ASKAPLOG_DEBUG_STR(logger, "All data are used for uv weight construction");
            }
            itsFirstAccumulatedVis = false;
        }

        // the tolerance is hard coded the same way as in the proper gridders
        // although in principle we can make it configurable. Also note that itsPointingTolerance
        // is used for different purpose (and with offset beams). But there is no huge reason why
        // these tolerances should be different. We just have to match whatever the ordinary gridder 
        // is doing, at least for now
        if (itsDoBeamAndFieldSelection && ((itsSelectedBeam != acc.feed1()(i)) ||
            (itsSelectedPointing.separation(acc.dishPointing1()(i)) > 1e-6))) {
            continue;
        }

        uvWeightRW = itsUVWeightBuilder->addWeight(acc.feed1()(i), itsCurrentField, itsSourceIndex);
        ASKAPDEBUGASSERT(uvWeightRW.uSize() == itsShape[0]);
        ASKAPDEBUGASSERT(uvWeightRW.vSize() == itsShape[1]);

        for (casacore::uInt chan=0; chan<nChan; ++chan) {
             const double reciprocalToWavelength = frequencyList[chan]/casacore::C::c;
             if (chan == 0) {
                // check for ridiculous frequency to pick up a possible error with input file,
                // not essential for processing as such
                ASKAPCHECK((reciprocalToWavelength>0.1) && (reciprocalToWavelength<30000),
                    "Check frequencies in the input file as the order of magnitude is likely to be wrong, "
                    "comment this statement in the code if you're trying something non-standard. Frequency = "<<
                    frequencyList[chan]/1e9<<" GHz");
             }
             /*
             // commented out until AXA-2485 is sorted out
             // Scale U,V to integer pixels, ignore fractional terms and dependence on oversampling factor in the current code (see AXA-2485)
             const double uScaled=reciprocalToWavelength * outUVW(i)(0) / itsUCellSize;
             const int iu = askap::nint(uScaled) + itsShape(0) / 2;
             const double vScaled=reciprocalToWavelength * outUVW(i)(1) / itsVCellSize;
             const int iv = askap::nint(vScaled) + itsShape(1) / 2;
             */
             // same approach as used with gridder which takes fractional term into account and assumes oversampling factor of 1.
             const double uScaled=reciprocalToWavelength * outUVW(i)(0) / itsUCellSize;
             int iu = askap::nint(uScaled);
             const int fracu = askap::nint(itsOversample*(double(iu) - uScaled));
             if (fracu < 0) {
                 iu += 1;
             } else if (fracu >= itsOversample) {
                 iu -= 1;
             }
             const double vScaled=reciprocalToWavelength * outUVW(i)(1) / itsVCellSize;
             int iv = askap::nint(vScaled);
             const int fracv = askap::nint(itsOversample*(double(iv) - vScaled));
             if (fracv < 0) {
                 iv += 1;
             } if (fracv >= itsOversample) {
                 iv -= 1;
             }
             iu += itsShape(0) / 2;
             iv += itsShape(1) / 2;

             // mimic the behaviour of the orginary gridder w.r.t. partial polarisation, i.e. ignore the whole sample
             bool allPolGood=true;
             for (casacore::uInt pol=0; pol<nPol; ++pol) {
                  //if (flagCube(i, chan, pol)) {
                  if (flagCube(pol, chan, i)) {
                      allPolGood=false;
                      break;
                  }
             }
             if (allPolGood && itsFreqMapper.isMapped(chan) &&
                 (iu >= 0) && (iv >= 0) && (iu < itsShape(0)) && (iv < itsShape(1))) {
                 // obtain which channel of the weight grid this accessor channel is mapped to
                 const int imageChan = itsFreqMapper(chan);
                 ASKAPDEBUGASSERT(imageChan < uvWeightRW.nPlane());
                 // there is a short-cut here as we ignore the actual polarisation product produced for weight calculation. If X and Y have 
                 // notably different noise this could get us into trouble.

                 // we can run it through the converter, say for stokes I, but the estimate probably won't be better anyway, so keep it simple
                 // the check for more than one polarisation could in principle be optimised (e.g. scaling factor could be defined outside the loop) or even removed
                 // completely as absolute scaling is not important for uv-weights (but it is handy to get the matching behaviour to the case when the weights are 
                 // generated during the gridding)
                 //const float visNoise = casacore::square(casacore::real(noiseCube(i, chan, 0))) + (nPol > 1u ? casacore::square(casacore::real(noiseCube(i, chan, nPol-1))) : 0.f);
                 const float visNoise = casacore::square(casacore::real(noiseCube(0, chan, i))) + (nPol > 1u ? casacore::square(casacore::real(noiseCube(nPol-1, chan, i))) : 0.f);
                 const float visNoiseWt = (visNoise > 0.) ? 1./visNoise : 0.;
                 ASKAPCHECK(visNoiseWt>0., "Weight is supposed to be a positive number; visNoiseWt="<<
                            visNoiseWt<<" visNoise="<<visNoise);

                 ASKAPDEBUGASSERT(iu < uvWeightRW.uSize());
                 ASKAPDEBUGASSERT(iv < uvWeightRW.vSize());

                 uvWeightRW(iu,iv,imageChan) += visNoiseWt;
             }
        }
   }
}

/// @brief checks whether the current field has been updated
/// @details See currentField for the description of limitations. This method detects field changes in the field pointing (and numbers them in the 
/// order they are encountered). If at a later stage we find that the fields need to be numbered in a particular way, this can be implemented.
/// @note To match implementation of the gridder classes, we detect changes in the pointing of the first encountered beam. It has implications if
/// either 3rd axis is operated in a non-tracking way or accessor row structure is different from one iteration to another. I (MV) suspect it was done
/// this way because in early days we're trying to simulate equatorial vs. alt-az mounts and, technically, physical beam pointing matters.
/// @param[in] acc input const accessor to analyse
void UVWeightGridder::indexField(const accessors::IConstDataAccessor &acc) const
{
  // the code below uses a different (and more simple) approach to that of AProjectGridderBase, but may cause problems if we ever process data
  // with off-axis beams and non-standard operation of 3rd axis (or from an alt-az telescope).

  ASKAPDEBUGASSERT(acc.nRow()>0);
  const casacore::MVDirection firstPointing = acc.pointingDir1()(0);

  // some speed-up can probably be achieved if we check itsCurrentField first (most likely the field won't change from accessor to accessor), but 
  // on the other hand, it would only matter if we have many fields in the same measurement set/iteration

  for (size_t field = 0; field < itsKnownPointings.size(); ++field) {
       if (firstPointing.separation(itsKnownPointings[field])<itsPointingTolerance) {
           itsCurrentField = field;
           return;
       }
  }

  itsCurrentField = itsKnownPointings.size();
  itsKnownPointings.push_back(firstPointing);
  ASKAPLOG_DEBUG_STR(logger, "Found new field " << itsCurrentField << " with beam "<<acc.feed1()(0)<<" at "<<
            printDirection(firstPointing));
}

/// @brief obtain the tangent point
/// @details  This method extracts the tangent point (reference position) from the
/// coordinate system.
/// @return direction measure corresponding to the tangent point
casacore::MVDirection UVWeightGridder::getTangentPoint() const
{
   ASKAPCHECK(itsAxes.hasDirection(),"Direction axis is missing. axes="<<itsAxes);
   const casacore::Vector<casacore::Double> refVal(itsAxes.directionAxis().referenceValue());
   ASKAPDEBUGASSERT(refVal.nelements() == 2);
   const casacore::Quantum<double> refLon(refVal[0], "rad");
   const casacore::Quantum<double> refLat(refVal[1], "rad");
   const casacore::MVDirection out(refLon, refLat);
   return out;
}

// (MV): code duplication with the gridder classes tells me that we probably need to move it to Axes (and getTangentPoint too)
// there is a complication with shape and padding, though

/// @brief obtain the centre of the image
/// @details This method extracts RA and DEC axes from itsAxes and
/// forms a direction measure corresponding to the middle of each axis.
/// @return direction measure corresponding to the image centre
casacore::MVDirection UVWeightGridder::getImageCentre() const
{
   ASKAPCHECK(itsAxes.hasDirection(),"Direction axis is missing. axes="<<itsAxes);
   casacore::MDirection out;
   casacore::Vector<casacore::Double> centrePixel(2);
   ASKAPDEBUGASSERT(itsShape.nelements()>=2);
   ASKAPDEBUGASSERT(paddingFactor()>0);
   for (size_t dim=0; dim<2; ++dim) {
        centrePixel[dim] = double(itsShape[dim])/2./double(paddingFactor());
   }
   // MV: note, there were experiments running gridder code under OpenMP which would hit the issue with the lack of
   // thread safety for casacore routines. To abstract this out, syncHelper was written but it is a bit ugly to
   // include it outside of TableVisGridder. If thread safety is required here, one would need to change this code
   // to the way similar to TableVisGridder and probably made the access to syncHelper via a proper singleton pattern.
   // At this stage, I don't think we ever run this part from multiple threads, so I use casacore pixel to world conversion
   // directly.
   ASKAPCHECK(itsAxes.directionAxis().toWorld(out, centrePixel),
        "Unable to obtain world coordinates for the centre of the image. Something is wrong with the coordinate system");
   return out.getValue();
}

} // namespace synthesis

} // namespace askap

