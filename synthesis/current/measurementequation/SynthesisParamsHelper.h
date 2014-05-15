/// @file
/// @brief Helper functions for dealing with Params for synthesis
///
/// Adds some useful functions specific to synthesis
/// @todo Function to output nicely formatted axes
/// @todo Functions to read/write images
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
#ifndef SYNSYNTHESISPARAMSHELPER_H_
#define SYNSYNTHESISPARAMSHELPER_H_

#include <casa/Arrays/Array.h>
#include <casa/Arrays/Slicer.h>

#include <fitting/Params.h>
#include <fitting/Axes.h>

#include <Common/ParameterSet.h>
#include <images/Images/TempImage.h>
#include <images/Images/ImageInterface.h>
#include <coordinates/Coordinates/Projection.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>


#include <imageaccess/IImageAccess.h>
#include <boost/shared_ptr.hpp>

namespace askap
{
  namespace synthesis
  {
  
    /// @brief Helper functions for synthesis processing using Params
    /// @ingroup measurementequation
    class SynthesisParamsHelper 
    {
      public:
        /// @brief setup image handler
        /// @details This method uses the factory to setup a helper class handling the
        /// operations with images (default is casa). It is necessary to call this method
        /// at least once before any read or write operation can happen.
        /// @param[in] parset a parset file containing parameters describing which image handler to use
        /// @note The key parameter describing the image handler is "imagetype". By default, the
        /// casa image handler is created (however, a call to this method is still required)
        static void setUpImageHandler(const LOFAR::ParameterSet &parset);

        /// @brief configure default frequency frame
        /// @details All code workes in a single frequency frame (convertions are done, if
        /// necessary when the data are read (using conversion mechanism provided by the accessor).
        /// A call to this method sets up new default.
        /// @param[in] frame reference frame to use for all created images
        static void setDefaultFreqFrame(const casa::MFrequency::Ref &frame);
       
        /// @brief obtain image handler
        /// @details For some operations it may be necessary to access the (global) instance of the
        /// image handler. This method allows that. An exception is thrown if no image handler has
        /// been previously set up.
        /// @return a reference to image handler
        static accessors::IImageAccess& imageHandler();
        
        /// @brief zero all free model images
        /// @details I (MV) hope that this method is temporary. In the current design of the code we need to 
        /// discard the solution of the model updates in the case of a dirty image. Otherwise, the restored image
        /// is wrong. This method iterates over all free model image parameters and sets them to 0.
        /// @param[in] params collection of parameters
        static void zeroAllModelImages(const askap::scimath::Params::ShPtr& params);        
        
        /// @brief set up images according to the parset file
		/// @param[in] params Images to be created here
		/// @param[in] parset a parset object to read the parameters from
		/// @note (MV)This method is probably a duplication of the one of 
		/// add methods - needs to be cleared
		/// (MV, dec 2008) not any more. With faceting it is handy to have a separate method
		static void setUpImages(const askap::scimath::Params::ShPtr& params, const LOFAR::ParameterSet &parset);
		
		/// @brief load images according to the parset file
		/// @details This method is somewhat analogous to setUpImages, but it loads the images
		/// from the disk instead of setting them up from the scratch. Encapsulation of all loading
		/// of multiple images in a single method is required to provide a seamless handling of
		/// the faceted image.
		/// @param[in] params Images to be created here
		/// @param[in] parset a parset object to read the parameters from		
		static void loadImages(const askap::scimath::Params::ShPtr& params, const LOFAR::ParameterSet &parset);
        
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
        static void copyComponent(const askap::scimath::Params::ShPtr &params,
                 const LOFAR::ParameterSet &parset, 
                 const std::string &srcName, const std::string &baseKey = "sources.");
        
        /// @brief check whether parameter list defines at least one component
        /// @details Parameter lists can have a mixture of components and
        /// images defined. This method checks whether the given parameter
        /// list defines at least one component.
        /// @param[in] params a shared pointer to the parameter container
        /// @return true, if at least one component is defined
        static bool hasComponent(const askap::scimath::Params::ShPtr &params);
       
        /// @brief check whether parameter list defines at least one image
        /// @details Parameter lists can have a mixture of components and
        /// images defined. This method checks whether the given parameter
        /// list defines at least one image.
        /// @param[in] params a shared pointer to the parameter container
        /// @return true, if at least one image is defined
        static bool hasImage(const askap::scimath::Params::ShPtr &params);
        
      
        /// @brief Add a parameter as an image
        /// @param[in] ip Parameters
        /// @param[in] name Name of parameter
        /// @param[in] direction Strings containing [ra, dec, frame] defining tangent point
        /// @param[in] cellsize Cellsize as a string e.g. [12arcsec, 12arcsec]
        /// @param[in] shape Number of pixels in RA and DEC e.g. [256, 256]
        /// @param[in] ewprojection If true, SCP or NCP variant of SIN projection will be used
        /// @param[in] freqmin Minimum frequency (Hz)
        /// @param[in] freqmax Maximum frequency (Hz)
        /// @param[in] nchan Number of spectral channels
        /// @param[in] stokes Polarisation frame (vector of stokes enums)
        /// @param[in] centreDir strings containing [ra,dec,frame] defining direction of the image centre
        ///            empty vector defaults to the tangent point
        static void add(askap::scimath::Params& ip, const string& name, 
          const vector<string>& direction, 
          const vector<string>& cellsize, 
          const vector<int>& shape,
          const bool ewprojection,
          const double freqmin, const double freqmax, const int nchan, 
          const casa::Vector<casa::Stokes::StokesTypes> &stokes,
          const vector<string>& centreDir = vector<string>());

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
        static void add(askap::scimath::Params& ip, const string& name, 
          const vector<string>& direction, 
          const vector<string>& cellsize, 
          const vector<int>& shape,
          const bool ewprojection,
          const double freqmin, const double freqmax, const int nchan,
          const casa::Vector<casa::Stokes::StokesTypes> &stokes,
          const int nfacets, const int facetstep);
        
        /// @brief add a parameter as a merged faceted image
        /// @details Each facet is represented by a number of independent parameters with
        /// the appropriate names. This method looks at the coordinate systems of all
        /// subimages and forms a parameter representing merged image. It can then be
        /// populated with the data from the appropriate slices.
        /// @param[in] ip parameters
        /// @param[in] name Base name of the parameter (i.e. without .facet.0.0)
        /// @param[in] nfacets number of facets defined
        static void add(askap::scimath::Params& ip, const string &name,
              const int nfacets); 
                            
        /// @brief helper method to clip the outer edges of the image
        /// @details For experiments with faceting we want to be able to clip the outer
        /// edges of each model image (beyond the facet step) to zero. This is one way to
        /// reduce cross-talk problem (when facets overlap). This method encapsulates all
        /// the required operations. It takes facet step from the fake image axis FACETSTEP
        /// and does nothing if such a parameter doesn't exist or is larger than the shape
        /// along the directional axes.
        /// @param[in] ip parameters
        /// @param[in] name full name of the image (i.e. with .facet.x.y for facets)
        static void clipImage(const askap::scimath::Params &ip, const string &name);
        
        
        /// @brief helper method to store restoring beam for an image
        /// @details We have to carry restore beam parameters together with the image.
        /// This is done by creating 2 fake axes MAJMIN (with start = maj and end = min)
        /// and PA with position angle. All angles are given in radians. The presence of
        /// this fake axes distinguishes a restored image from model image. Restored image
        /// will have units Jy/beam instead of Jy/pixel and beam info will be added to the
        /// image.
        /// @param[in] ip parameters
        /// @param[in] name full name of the parameter representing this image
        /// @param[in] beam major, minor axes and position anlge as quantities
        static void setBeam(askap::scimath::Params &ip, const string &name,
                            const casa::Vector<casa::Quantum<double> > &beam);
        
        /// @brief find a parameter representing a PSF
        /// @details If multiple PSF parameters are present, the first encountered is returned
        /// @param[in] ip parameters
        /// @return full name of some PSF parameter or an empty string if it is not found
        static std::string findPSF(const askap::scimath::Params &ip);
        
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
        static casa::Vector<casa::Quantum<double> > fitBeam(const askap::scimath::Params &ip, const double cutoff = 0.05,
                                                            const std::string &name = "");
              
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
        static casa::Array<double> getFacet(askap::scimath::Params &ip, const string &name);                       
                  
        /// @brief Add a set of parameters from a parset
        /// @param ip Parameters
        /// @param parset ParameterSet
        /// @param baseKey basekey for parameters e.g. "Images."
        static void add(askap::scimath::Params& ip,
          const LOFAR::ParameterSet& parset,
          const std::string& baseKey);
          
        /// @brief Get a parameter from an image
        /// @param ip Parameters
        /// @param name Name of parameter
        /// @param imagename Name of image file
        static void loadImageParameter(askap::scimath::Params& ip, const string& name,
          const string& imagename);
          
        /// @brief Get parameters corresponding to all facets from a CASA image
        /// @param[in] ip Parameters
        /// @param[in] name Base name of the parameter (.facet.x.y will be added)
        /// @param[in] fileName Base name of the image file (.facet.x.y will be added)
        /// @param[in] nfacets Number of facets on each axis (assumed the same for both axes)
        static void getMultiFacetImage(askap::scimath::Params &ip, const string &name,
           const string &fileName, const int nfacets);
        
        /// @brief Save a parameter as a CASA image
        /// @param ip Parameters
        /// @param name Name of parameter
        /// @param imagename Name of image file
        static void saveImageParameter(const askap::scimath::Params& ip, const string& name,
          const string& imagename);
                  
        /// @brief Copy a parameter to a CASA TempImage
        /// Note that this will be a reference if possible
        /// @param ip Parameters
        /// @param name Name of parameter
        static boost::shared_ptr<casa::TempImage<float> > 
          tempImage(const askap::scimath::Params& ip, 
          const string& name);
       
        /// @brief Create a coordinate system for a parameter
        /// @param ip Parameters
        /// @param name Name of parameter
        static casa::CoordinateSystem 
          coordinateSystem(const askap::scimath::Params& ip, 
          const string& name);
                      
        /// @brief Create a direction coordinate for a parameter
        /// @param ip Parameters
        /// @param name Name of parameter
        static casa::DirectionCoordinate 
          directionCoordinate(const askap::scimath::Params& ip, 
          const string& name);
       
        /// @brief Update a parameter from an image
        /// @param ip Parameters
        /// @param name Name of parameter
        /// @param image Image to be drawn from 
        static void update(askap::scimath::Params& ip, const string& name, 
          const casa::ImageInterface<float>& image);
        
        /// @brief a helper template method to check whether the element is
        /// present in a container.
        /// @details This method is used to make the code more readable. It is very
        /// generic and can be moved to Base if needed elsewhere.
        /// @param[in] cont a container (stl)
        /// @param[in] val value to check
        /// @return true if an element equal to val is present in the container
        template<typename C, typename V>
        static bool hasValue(const C &cont, const V &val)
           { return cont.find(val) != cont.end(); }        
        
        /// @brief A helper method to parse strings of quantities
        /// @details Many parameters in parset file are given as quantities or
        /// vectors of quantities, e.g. [8.0arcsec,8.0arcsec]. This method allows
        /// to parse vector of strings corresponding to such parameter and return
        /// a vector of double values in the required units.
        /// @param[in] strval input vector of strings
        /// @param[in] unit required units (given as a string)
        /// @return vector of doubles with converted values
        static std::vector<double> convertQuantity(const std::vector<std::string> &strval,
                       const std::string &unit);

        /// @brief A helper method to parse string of quantities
        /// @details Many parameters in parset file are given as quantities or
        /// vectors of quantities, e.g. 8.0arcsec. This method allows
        /// to parse a single string corresponding to such a parameter and return
        /// a double value converted to the requested units.
        /// @param[in] strval input string
        /// @param[in] unit required units (given as a string)
        /// @return converted value
        static double convertQuantity(const std::string &strval,
                       const std::string &unit);                       
                               
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
        static void listFacets(const std::vector<std::string> &names,
                               std::map<std::string, int> &facetmap);

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
        static void listTaylor(const std::vector<std::string> &names,
                               std::map<std::string, int> &taylormap);

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
        static casa::Slicer facetSlicer(const askap::scimath::Params& ip, const std::string &name,
                              const casa::DirectionCoordinate &dc);

        /// @brief make a merged image parameter covering all given facets
        /// @details This method is very similar to another version of the add method which creates
        /// an image paramters covering all facets named in the appropriate fashion. Although doing 
        /// essentially the same job, this method works with any images, i.e. they are not necessarily
        /// regularly spaced and appropriately named facets. With time we can probably change how we
        /// do faceting and retire the old methods.
        /// @param[in] ip parameters
        /// @param[in] names names of all images to merge
        /// @param[in] mergedName name of the image to create
        static void add(askap::scimath::Params& ip, const std::vector<std::string> &names,
              const std::string &mergedName); 
        
    protected:
        /// @brief helper method to get projection
        /// @details We support both standard SIN projection and SCP/NCP variants (for East-West arrays).
        /// This method encapsulates the logic and returns a projection class
        /// @param[in] ewprojection true for SCP/NCP variant, false otherwise
        /// @param[in] dec declination in radians (unused for standard SIN projection)
        /// @return casa::Projection class
        static casa::Projection getProjection(const bool ewprojection, const double dec = 0.);
        
    private:    
        /// @brief image accessor
        static boost::shared_ptr<accessors::IImageAccess> theirImageAccessor;              

        /// @brief default frequency frame
        static casa::MFrequency::Ref theirFreqFrame;
    };

  }
}
#endif
