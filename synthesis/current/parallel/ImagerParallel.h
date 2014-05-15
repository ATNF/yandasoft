/// @file
///
/// ImagerParallel: Support for parallel applications using the measurement equation
/// classes.
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
///
#ifndef ASKAP_SYNTHESIS_IMAGERPARALLEL_H_
#define ASKAP_SYNTHESIS_IMAGERPARALLEL_H_

// ASKAPsoft includes
#include <askapparallel/AskapParallel.h>
#include <fitting/Solver.h>
#include <Common/ParameterSet.h>

// Local package includes
#include <parallel/MEParallelApp.h>
#include <measurementequation/IMeasurementEquation.h>
#include <calibaccess/ICalSolutionConstSource.h>

namespace askap
{
  namespace synthesis
  {
    /// @brief Support for parallel algorithms implementing imaging
    ///
    /// @details Provides imaging using the measurement equation framework.
    ///
    /// The control parameters are specified in a parset file. For example:
    /// @code
    ///  	Cimager.datacolumnset           = DATACOL     # default is DATA
    ///     Cimager.memorybuffers           = true      # default is false, i.e. write scratch data to the MS
    ///  	Cimager.dataset                 = [data/spw_1/sim.ms]
    ///  	#Feed                           = 5
    ///
    ///  	Cimager.Images.Names                    = [image.i.sim]
    ///  	Cimager.Images.image.i.sim.shape                = [4096,4096]
    ///  	Cimager.Images.image.i.sim.cellsize     = [6arcsec, 6arcsec]
    ///  	Cimager.Images.image.i.sim.frequency    = [1.164e9, 1.420e9]
    ///  	Cimager.Images.image.i.sim.nchan                = 1
    ///  	Cimager.Images.image.i.sim.direction    = [12h30m00.00, -45.00.00.00, J2000]
    ///
    ///  	#Cimager.gridder                          = WProject
    ///  	Cimager.gridder.WProject.wmax            = 8000
    ///  	Cimager.gridder.WProject.nwplanes        = 64
    ///  	Cimager.gridder.WProject.oversample     = 1
    ///  	Cimager.gridder.WProject.cutoff         = 0.001
    ///
    ///  	Cimager.gridder                                 = AntennaIllum
    ///  	Cimager.gridder.AntennaIllum.wmax               = 8000
    ///  	Cimager.gridder.AntennaIllum.nwplanes           = 64
    ///  	Cimager.gridder.AntennaIllum.oversample         = 1
    ///  	Cimager.gridder.AntennaIllum.diameter           = 12m
    ///  	Cimager.gridder.AntennaIllum.blockage           = 1m
    ///  	Cimager.gridder.AntennaIllum.maxfeeds           = 18
    ///  	Cimager.gridder.AntennaIllum.maxsupport         = 256
    ///
    ///  	Cimager.solver                                  = Dirty
    ///  	Cimager.solver.Dirty.cycles                     = 1
    ///  	Cimager.solver.Dirty.threshold                  = 10%
    ///  	Cimager.solver.Dirty.verbose                    = True
    ///
    ///  	Cimager.solver.Clean.algorithm                  = MultiScale
    ///  	Cimager.solver.Clean.niter                      = 1000
    ///  	Cimager.solver.Clean.gain                       = 0.7
    ///  	Cimager.solver.Clean.cycles                     = 5
    ///  	Cimager.solver.Clean.verbose                    = True
    ///
    ///  	Cimager.restore                                 = True
    ///  	Cimager.restore.beam                            = [30arcsec, 30arcsec, 0deg]
    /// @endcode
    /// @ingroup parallel
    class ImagerParallel : public MEParallelApp
    {
  public:

      /// @brief Constructor from ParameterSet
      /// @details The parset is used to construct the internal state. We could
      /// also support construction from a python dictionary (for example).
      /// The command line inputs are needed solely for MPI - currently no
      /// application specific information is passed on the command line.
      /// @param comms communication object 
      /// @param parset ParameterSet for inputs
      ImagerParallel(askap::askapparallel::AskapParallel& comms,
          const LOFAR::ParameterSet& parset);

      /// @brief Calculate the normalequations (runs in the prediffers)
      /// @details ImageFFTEquation and the specified gridder (set in the parset
      /// file) are used to calculate the normal equations. The image parameters
      /// are defined in the parset file.
      virtual void calcNE();

      /// @brief Solve the normal equations (runs in the solver)
      /// @details Either a dirty image can be constructed or the 
      /// multiscale clean can be used, as specified in the parset file.
      virtual void solveNE();

      /// @brief Write the results (runs in the solver)
      /// @details The model images are written as AIPS++ images. In addition,
      /// the images may be restored using the specified beam.
      /// @param[in] postfix this string is added to the end of each name
      /// (used to separate images at different iterations)
      virtual void writeModel(const std::string &postfix = std::string());

      /// @brief make sensitivity image
      /// @details This is a helper method intended to be called from writeModel. It
      /// converts the given weights image into a sensitivity image and exports it.
      /// This method is intended to be called if itsExportSensitivityImage is true.
      /// @param[in] wtImage weight image parameter name
      void makeSensitivityImage(const std::string &wtImage) const;

      /// @brief Helper method to zero all model images
      /// @details We need this for dirty solver only, as otherwise restored image 
      /// (which is crucial for faceting) will be wrong.
      void zeroAllModelImages() const;

  protected:
      
      /// @brief a helper method to extract peak residual
      /// @details This object actually manipulates with the normal equations. We need
      /// to be able to stop iterations on the basis of maximum residual, which is a 
      /// data vector of the normal equations. This helper method is designed to extract
      /// peak residual. It is then added to a model as a parameter (the model is 
      /// shipped around).
      /// @return absolute value peak of the residuals corresponding to the current normal
      /// equations
      double getPeakResidual() const;
      
      /// @brief helper method to indentify model parameters to broadcast
      /// @details We use itsModel to buffer some derived images like psf, weights, etc
      /// which are not required for prediffers. It just wastes memory and CPU time if
      /// we broadcast them. At the same time, some auxilliary parameters like peak
      /// residual value need to be broadcast (so the major cycle can terminate in workers).
      /// This method returns the vector with all parameters to be broadcast. By default
      /// it returns all parameter names, so it is overridden here to broadcast only
      /// model images and the peak_residual metadata.
      /// @return a vector with parameters to broadcast
      virtual std::vector<std::string> parametersToBroadcast() const;
      
      
  private:

      /// Calculate normal equations for one data set
      /// @param ms Name of data set
    /// @param discard Discard old equation?
      void calcOne(const string& dataset, bool discard=false);

      /// Do we want a restored image?
      bool itsRestore;
      
      /// @brief Do we want to keep scratch buffers in memory instead of writing them in a subtable?
      /// @details Turining this flag to true allows to work with a read-only dataset
      bool itsUseMemoryBuffers;
            
      /// @brief solution source to get calibration data from
      /// @details This object is initialised by workers. It knows how to
      /// retrieve calibration solutions (from a parset file, casa table or a database).
      /// Uninitialised shared pointer means that no calibration is done and the data are
      /// imaged as they are.
      boost::shared_ptr<accessors::ICalSolutionConstSource> itsSolutionSource;
            
      /// @brief void measurement equation
      /// @details Does nothing, just returns calls to predict and 
      /// calcNormalEquations. This shared pointer is initialized at the
      /// first use and then passed as the perfect MeasurementEquation, if
      /// required. We could have created a brand new object each time.
      boost::shared_ptr<IMeasurementEquation> itsVoidME;
      
      /// @brief export sensitivity image?
      /// @details If true, the weight image is converted to sensitivity image and 
      /// exported in addition to the original weight image.
      bool itsExportSensitivityImage;
      
      /// @brief fraction of peak sensitivity to crop on
      /// @details We normally don't want to divide by small numbers when calculating
      /// sensitivity images. This field gives the fraction of the maximum weight
      /// below which the sensitivity image will be set to 0.
      double itsExpSensitivityCutoff;
    };

  }
}
#endif
