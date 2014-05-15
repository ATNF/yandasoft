/// @file Simulator.h
///
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>

#ifndef ASKAP_SYNTHESIS_SIMULATOR_H
#define ASKAP_SYNTHESIS_SIMULATOR_H

//# Includes
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>
#include <casa/Quanta/Quantum.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MDirection.h>
#include <ms/MeasurementSets/MeasurementSet.h>

namespace askap
{
    namespace synthesis
    {
        /// @brief Simulator for synthesis observing 
        ///
        /// @details Cloned from casa::NewMSSimulator to allow substantial changes.
        /// @ingroup simulation

        class Simulator
        {
            public:

                /// Constructor from name and tile shape
                /// specification. This is used to define
                /// the tiles holdig the DATA column.
                /// Note that if the actual nr of corr or chan
                /// is less than given here, the actual tile
                /// size will be smaller.
                /// @param msname Name of measurement set to be constructed
                /// @param tileSize  Tile size in bytes
                /// @param ncorrTile Nr of correlations per tile
                /// @param nchanTile Nr of channels per tile
                explicit Simulator(const casa::String& msname,
                        int tileSize=32768,
                        int ncorrTile=4, int nchanTile=32);

                /// Constructor from existing MS
                /// @param ms Existing MeasurementSet object
                explicit Simulator(casa::MeasurementSet& ms);

                // Destructor
                ~Simulator();

                /// @brief Set the antenna and array data. 
                /// @details These are written immediately to the
                /// existing MS. The same model is used for the other init information
                /// @param telname Telescope name e.g. "ASKAP"
                /// @param x X coordinate for antennas (m)
                /// @param y Y coordinate for antennas (m)
                /// @param z Z coordinate for antennas (m)
                /// @param dishdiameter Dish diameters (m)
                /// @param offset Offset (m)
                /// @param mount Mount type (alt-az or equatorial)
                /// @param name Antenna name e.g. ASKAP23
                /// @param coordsystem Coordinate system local or global
                /// @param location MPosition of local coordinate system
                void initAnt(const casa::String& telname,
                        const casa::Vector<double>& x, const casa::Vector<double>& y,
                        const casa::Vector<double>& z,
                        const casa::Vector<double>& dishdiameter,
                        const casa::Vector<double>& offset,
                        const casa::Vector<casa::String>& mount,
                        const casa::Vector<casa::String>& name,
                        const casa::String& coordsystem, const casa::MPosition& location);

                /// @brief Set the observed fields
                /// @param sourcename Source name
                /// @param sourcedirection MDirection for source
                /// @param calcode Calibrator code
                void initFields(const casa::String& sourcename,
                            const casa::MDirection& sourcedirection,
                            const casa::String& calcode);

                /// @brief Set the Feeds
                /// @param mode "Perfect X Y" or "Perfect R L"
                /// @param x Offset of feed (radians)
                /// @param y Offset of feed (radians)
                /// @param pol Polarisation of feed e.g. "R L"
                void initFeeds(const casa::String& mode, const casa::Vector<double>& x,
                            const casa::Vector<double>& y,
                            const casa::Vector<casa::String>& pol);

                /// @brief Set the spectral windows information
                /// @param spwindowname Name of spectral window e.g. "Continuum6"
                /// @param nchan Number of channels e.g. 16
                /// @param startfreq Starting frequency e.g. "1.420GHz"
                /// @param freqinc Frequency increment e.g. "1MHz"
                /// @param freqres Frequency resolution e.g. "1MHz"
                /// @param stokes Stokes to be measured e.g. "XX XY YX YY"
                void initSpWindows(const casa::String& spwindowname, const int& nchan,
                        const casa::Quantity& startfreq, const casa::Quantity& freqinc,
                        const casa::Quantity& freqres, const casa::String& stokes);

                /// @brief Set maximum tolerable blockage before flagging
                /// @param fraclimit Maximim tolerable blockage (fraction) e.g. 0.1
                void setFractionBlockageLimit(const double fraclimit)
                {
                    fractionBlockageLimit_p = fraclimit;
                }

                /// @brief Set minimum allowed elevation before flagging
                /// @param ellimit Minimum allowed elevation e.g. "8deg"
                void setElevationLimit(const casa::Quantity& ellimit)
                {
                    elevationLimit_p = ellimit;
                }

                /// @brief Set autocorrelation weight (zero for no auto's)
                /// @param autocorrwt Weight of autocorrelation data
                void setAutoCorrelationWt(const float autocorrwt)
                {
                    autoCorrelationWt_p = autocorrwt;
                }

                /// @brief set noise rms (used to scale SIGMA column)
                /// @param[in] rms single visibility rms in Jy
                void setNoiseRMS(const float rms)
                {
                    itsNoiseRMS = rms;
                }

                /// @brief set relative Tsys/eff per antenna
                /// @details This method allows to simulate different Tsys 
                /// and efficiencies per antenna. The noise will be scaled
                /// up with the given factors. Empty array switches the
                /// scaling off and has the same effect as an array of 1
                /// with the number of elements equal to the number of antennas
                void setRelAntennaWeight(const casa::Vector<double> &wt)
                {
                   itsRelAntennaWeight.assign(wt.copy());
                }
                
                /// @brief access vector of relative antenna weights
                /// @details At the moment, we only use this to distinguish
                /// between the constant and variable noise cases
                /// (zero and non-zero length of the vector).
                /// @return const reference to the vector with relative antenna weights 
                const casa::Vector<double>& relAntennaWeight() const { return itsRelAntennaWeight; }

                /// @brief Set meaning of times for observing
                /// @param integrationtime Integration time e.g. "10s"
                /// @param usehourangles Use hour angles for observing times?
                /// @param reftime Reference time to which all times are referred
                void settimes(const casa::Quantity& integrationtime,
                        const bool usehourangles, const casa::MEpoch& reftime);

                /// @brief Observe a given source 
                /// @param sourcename Source to be observed (as defined in initFields)
                /// @param spwindowname Spectral window (as defined in InitSpw)
                /// @param starttime Starting time e.g. "1h"
                /// @param stoptime Stopping time e.g. "1h30m"
                /// @details Generate the empty data rows for the observing. All the
                /// relevant information must have been filled in by previous init
                /// calls. 
                void observe(const casa::String& sourcename,
                        const casa::String& spwindowname, const casa::Quantity& starttime,
                        const casa::Quantity& stoptime);

                /// @brief return area times sqrt(bandwidth*int_time)
                /// @details This quantity is used for automatic noise estimates. It is 
                /// composed from itsChanBandwidthForNoise and itsDishDiamForNoise.
                /// An exception is thrown if either array is inhomogeneous or 
                /// multiple spectral resolutions are simulated. This method is supposed
                /// to be called when the simulator is fully defined.
                /// @return antenna area (m^2) multiplied by square root of the product of bandwidth(Hz) 
                /// and integration time (s)
                double areaTimesSqrtBT() const;

                /// Convert local coordinates to global
                /// @param xreturned Converted x coordinate
                /// @param yreturned Converted y coordinate
                /// @param zreturned Converted z coordinate
                /// @param location Reference location
                /// @param xin Input x coordinate
                /// @param yin Input y coordinate
                /// @param zin Input z coordinate
                static void local2global(casa::Vector<double>& xreturned,
                        casa::Vector<double>& yreturned, casa::Vector<double>& zreturned,
                        const casa::MPosition& location, const casa::Vector<double>& xin,
                        const casa::Vector<double>& yin, const casa::Vector<double>& zin);

                /// Convert global coordinates to local
                /// @param xreturned Converted x coordinate
                /// @param yreturned Converted y coordinate
                /// @param zreturned Converted z coordinate
                /// @param location Reference location
                /// @param xin Input x coordinate
                /// @param yin Input y coordinate
                /// @param zin Input z coordinate
                static void longlat2global(casa::Vector<double>& xreturned,
                        casa::Vector<double>& yreturned, casa::Vector<double>& zreturned,
                        const casa::MPosition& location, const casa::Vector<double>& xin,
                        const casa::Vector<double>& yin, const casa::Vector<double>& zin);
            protected:
            
                /// Returns the fractional blockage of one antenna by another
                /// @param fraction1 Fraction of 1 blocked by 2
                /// @param fraction2 Fraction of 2 blocked by 1
                /// @param uvw UVW coordinates (same units as diameter)
                /// @param diam1 Diameter of antenna 1
                /// @param diam2 Diameter of antenna 2
                static void blockage(double &fraction1, double &fraction2,
                        const casa::Vector<double>& uvw, // uvw in same units as diam!
                        const double diam1, const double diam2);

                /// Provide nicely formatted direction
                /// @param direction Direction to be formatted
                static casa::String formatDirection(const casa::MDirection& direction);

                /// Provide nicely formatted time
                /// @param time Time to be formatted
                static casa::String formatTime(const double time);

            private:

                /// @brief Copy constructor
                /// @details made private to avoid it being called inadvertently
                /// @param mss Simulator to be copied
                explicit Simulator(const Simulator & mss);

                /// @brief assignment operator
                /// @details made private to avoid it being called inadvertently
                /// @param mss Simulator to be assigned from 
                /// @return reference to itself
                Simulator & operator=(const Simulator & mss);

                /// Fractional blockage limit
                double fractionBlockageLimit_p;
                /// Elevation limit
                casa::Quantity elevationLimit_p;
                /// Autocorrelation weight
                float autoCorrelationWt_p;
                /// Telescope name
                casa::String telescope_p;
                /// Integration time as a Quantity
                casa::Quantity qIntegrationTime_p;
                /// Use hour angles?
                bool useHourAngle_p;
                /// Is the hour angle defined?
                bool hourAngleDefined_p;
                /// Reference time
                casa::MEpoch mRefTime_p;
                /// Offset time as a double
                double t_offset_p;

                /// Measurement set points
                casa::MeasurementSet* ms_p;


                /// Restore default values
                void defaults();
                
                /// @brief dish diameter (to get noise)
                /// @details We have a simplified noise model at the moment with the constant noise
                /// for all visibilities. This data field contains the dish diameter in metres used
                /// to estimate noise. Non-positive number means that the array is inhomogeneous
                /// (this case is not supported yet and will cause an exception if an automatic noise
                /// estimate is requested). 
                double itsDishDiamForNoise;
                
                /// @brief channel bandwidth (to get noise)
                /// @details We have a simplified noise model at the moment with the constant noise
                /// for all visibilities. This data field contains the bandwidth of a single channel in Hz
                /// used to estimate noise. Non-positive number means that the observations contain a
                /// number of setups with different frequency increments. This case is not yet
                /// supported and will cause an exception if an automatic noise estimate is
                /// requested.
                /// @note The value is initialised with -100 in the constructor to be able to distinguish 
                /// uninitialised value when no spectral windows are processed and mismatching resolutions
                /// (the value is -1). These are internal flags, which are completely hidden from the user of
                /// this class.
                double itsChanBandwidthForNoise;

                /// @brief noise rms to scale SIGMA column
                /// @details For a proper theoretical noise estimate SIGMA column must be 
                /// initialised with the data in real units (only relative values, which do not change,
                /// matter for the final image). This field is 1. by default (to mimic the old behavior),
                /// but can be set to the proper noise figure via the appropriate setter method.
                double itsNoiseRMS;

                /// @brief relative Tsys/eff per antenna
                /// @details It is used to apply different scaling factor (to itsNoiseRMS) per antenna
                /// Empty array means that the numbers are identical for all antennas.
                casa::Vector<double> itsRelAntennaWeight;
        };

    }

}
#endif
