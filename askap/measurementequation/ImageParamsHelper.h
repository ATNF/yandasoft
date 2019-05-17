/// @copyright (c) 2007 CSIRO
/// @brief Helper class for dealing with Params representing images
/// @details Working on the faceting, it was found that a parser for
/// image parameter name was required. It should return a number of values, so a 
/// separate class seems to be a better alternative than a static member of the
/// existing SynthesisParamsHelper class. Some methods from the latter will probably
/// migrate eventually into this class.
///  
///
///
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
///
/// @author Max Voronkov <maxim.voronkov@csiro.au>
///

#ifndef IMAGE_PARAMS_HELPER_H
#define IMAGE_PARAMS_HELPER_H

#include <string>

namespace askap {

namespace synthesis {

/// @brief Helper class for dealing with Params representing images
/// @details Working on the faceting, it was found that a parser for
/// image parameter name was required. It should return a number of values, so a 
/// separate class seems to be a better alternative than a static member of the
/// existing SynthesisParamsHelper class. Some methods from the latter will probably
/// migrate eventually into this class.
/// @todo improve parsing to add polarisation/Urvashi's decomposition into Taylor series
/// @ingroup measurementequation
class ImageParamsHelper {
public:
   /// @brief empty constructor
   /// @details Full name to be specified later. This method of construction doesn't produce
   /// a valid object until parse method is called.
   ImageParamsHelper();
   
   /// @brief constructor with immediate parsing of a full name
   /// @details This version constructs the object and populates all fields with the parse
   /// results.
   /// @param[in] name full name to parse
   ImageParamsHelper(const std::string &name);
   
   /// @brief direct constructor of a facet name from constituents
   /// @details This method constructs the object directly from the actual name
   /// of the image and facet indices.
   /// @param[in] name actual name of the image (without suffixes)
   /// @param[in] xFacet facet index along the first axis
   /// @param[in] yFacet facet index along the second axis
   ImageParamsHelper(const std::string &name, int xFacet, int yFacet);
   
   /// @brief direct constructor of a taylor term from constituents
   /// @details This method constructs the object directly from the actual name
   /// of the image and given order.
   /// @param[in] name actual name of the image (without suffixes)
   /// @param[in] order order in the Taylor series
   ImageParamsHelper(const std::string &name, int order);

   /// @brief direct constructor of a faceted taylor term from constituents
   /// @details This method constructs the object directly from the actual name
   /// of the image, given order and facet indices.
   /// @param[in] name actual name of the image (without suffixes)
   /// @param[in] order order in the Taylor series
   /// @param[in] xFacet facet index along the first axis
   /// @param[in] yFacet facet index along the second axis
   ImageParamsHelper(const std::string &name, int order, int xFacet, int yFacet);
      
   /// @brief parse the given string
   /// @param[in] name full name to parse
   void parse(const std::string &name);
  
   /// @brief obtain the actual name of the parameter without all suffixes
   /// @return the name without suffixes
   inline const std::string& name() const { return itsName;}
   
   /// @brief obtain the full name of the image parameter
   /// @details This method composes the full name of the parameter from 
   /// the data stored internally. This returned full name should be the same 
   /// as one passed in the parse method or in the constructor. This method can
   /// be useful if this object is constructed directly without parsing a 
   /// string and effectively represents a reverse operation.
   /// @return full name of the image parameter
   std::string paramName() const;
   
   /// @brief obtain the full suffix
   /// @details This method composes the suffix from the facet and taylor 
   /// term. The result returned by paramName is just name()+suffix().
   /// We need the suffix separately from the full name to be able to
   /// index model in MSF simulations (the model may have a different name
   /// than the field).
   std::string suffix() const;
   
   /// @brief obtain the name of the image with just a facet suffix
   /// @details To have MSMFS algorithm working with facets one needs
   /// to be able to extract the name without suffix corresonding to
   /// the Taylor decomposition, but with the faceting suffixes preserved,
   /// if present. This method forms such a name.
   /// @return the name of the parameter without Taylor suffixes
   /// @note Facet suffix will not be added if the image is not a facet
   std::string facetName() const;

   /// @brief obtain the name of the image with just a taylor suffix
   /// @details To have MSMFS algorithm working with facets one needs
   /// to be able to extract the name without suffixes corresonding to
   /// facet, but with the taylor term suffixes preserved,
   /// if present. This method forms such a name.
   /// @return the name of the parameter without facet suffixes
   /// @note Taylor suffix will not be added if the image is not a Taylor term
   std::string taylorName() const;

   /// @brief make this object a facet
   /// @details It is sometimes necessary to merge faceting suffixes and
   /// Taylor term suffix. This method makes the current image a facet with
   /// given indices. 
   /// @param[in] xFacet facet index along the first axis
   /// @param[in] yFacet facet index along the second axis
   void makeFacet(int xFacet, int yFacet);

   /// @brief make this object a facet
   /// @details It is sometimes necessary to merge faceting suffixes and
   /// Taylor term suffix. This method makes the current image a facet with
   /// given indices. 
   /// @param[in] order order in the Taylor series
   void makeTaylorTerm(int order);
   
   /// @brief obtain the facet number along the first axis
   /// @return the facet number along the first axis
   int facetX() const;  

   /// @brief obtain the facet number along the second axis
   /// @return the facet number along the second axis
   int facetY() const;  
   
   /// @brief check whether this parameter corresponds to a facet
   /// @details
   /// @return true if this is a facet
   inline bool isFacet() const { return itsFacetX >= 0;}
   
   /// @brief check whether this object is valid
   /// @details Constructed with a default constructor an instance of this
   /// class can not be used. This method would return false for such object.
   /// @return true, if object was set up with some parameter name
   inline bool isValid() const { return itsFacetX != -2; } 
   
   /// @brief check whether this parameter corresponds to a Taylor term
   /// @details
   /// @return true if this is a term of a Taylor series
   inline bool isTaylorTerm() const { return itsOrder >= 0; }
   
   /// @brief obtain the order of the Taylor term
   /// @return the order of the Taylor term represented by this parameter
   int order() const;   
   
private:
   /// @brief name of the current parameter (cut before all suffixes)
   std::string itsName;

   /// @brief facet number along the first axis
   /// @note The value is negative for non-faceted image
   int itsFacetX;

   /// @brief facet number along the second axis
   /// @note The value is negative for non-faceted image
   int itsFacetY;
   
   /// @brief order in the Taylor series
   /// @note This value is negative for images which are
   /// not part of Taylor series
   int itsOrder;   
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef IMAGE_PARAMS_HELPER_H
