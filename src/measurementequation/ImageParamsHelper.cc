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

#include <measurementequation/ImageParamsHelper.h>

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

namespace askap {

namespace synthesis {

/// @brief empty constructor
/// @details Full name to be specified later. This method of construction doesn't produce
/// a valid object until parse method is called.
ImageParamsHelper::ImageParamsHelper() : itsFacetX(-2), itsFacetY(-2), itsOrder(-2) {}
   
/// @brief constructor with immediate parsing of a full name
/// @details This version constructs the object and populates all fields with the parse
/// results.
/// @param[in] name full name to parse
ImageParamsHelper::ImageParamsHelper(const std::string &name) : itsFacetX(-2), itsFacetY(-2), 
                                                                itsOrder(-2)
{
  parse(name);
}

/// @brief direct constructor of a facet name from constituents
/// @details This method constructs the object directly from the actual name
/// of the image and facet indices.
/// @param[in] name actual name of the image (without suffixes)
/// @param[in] xFacet facet index along the first axis
/// @param[in] yFacet facet index along the second axis
ImageParamsHelper::ImageParamsHelper(const std::string &name, int xFacet, int yFacet) :
              itsName(name), itsFacetX(xFacet), itsFacetY(yFacet), itsOrder(-2)
{  
}

/// @brief direct constructor of a taylor term from constituents
/// @details This method constructs the object directly from the actual name
/// of the image and given order.
/// @param[in] name actual name of the image (without suffixes)
/// @param[in] order order in the Taylor series
ImageParamsHelper::ImageParamsHelper(const std::string &name, int order) : itsName(name),
              itsFacetX(-1), itsFacetY(-1), itsOrder(order)
{
}

/// @brief direct constructor of a faceted taylor term from constituents
/// @details This method constructs the object directly from the actual name
/// of the image, given order and facet indices.
/// @param[in] name actual name of the image (without suffixes)
/// @param[in] order order in the Taylor series
/// @param[in] xFacet facet index along the first axis
/// @param[in] yFacet facet index along the second axis
ImageParamsHelper::ImageParamsHelper(const std::string &name, int order, int xFacet, int yFacet) :
              itsName(name), itsFacetX(xFacet), itsFacetY(yFacet), itsOrder(order)
{
}

   
/// @brief parse the given string
/// @param[in] name full name to parse
void ImageParamsHelper::parse(const std::string &name) 
{
  // suffixes should go in the order: taylor, facet
  size_t pos = name.rfind(".facet.");
  if (pos == std::string::npos) {
      // this is not a faceted image, just set flags and copy full name
      itsFacetX = -1;
      itsFacetY = -1;
      itsName = name;                      
  } else {
      // this is a single facet, we have to extract indices 
      itsName = name.substr(0,pos);
      pos+=7; // to move it to the start of numbers
      ASKAPCHECK(pos < name.size(), 
         "Name of the faceted parameter should contain facet indices at the end, you have "<<name);
      const size_t pos2 = name.find(".",pos);
      ASKAPCHECK((pos2 != std::string::npos) && (pos2+1<name.size()) && (pos2!=pos), 
          "Two numbers are expected in the parameter name for the faceted image, you have "<<name);
      itsFacetX = utility::fromString<int>(name.substr(pos,pos2-pos));
      itsFacetY = utility::fromString<int>(name.substr(pos2+1));
  }
  // facet-related part shoud now be stripped off, if presented
  pos = itsName.rfind(".taylor.");
  if (pos == std::string::npos) {
      // this is not a term of the Taylor series, just leave the name as is and set the flags
      itsOrder = -1;
  } else {
      // this is a term of a Taylor series
      const std::string newName = itsName.substr(0,pos);
      pos += 8; // to move it to the start of the number
      ASKAPCHECK(pos < itsName.size(),
          "Name of the parameter representing a taylor term should contain order at the end, you have "<<name); 
      itsOrder = utility::fromString<int>(itsName.substr(pos));
      itsName = newName;
  }
  // further parsing of the parameter name, if necessary, should go here.  
}

/// @brief obtain the facet number along the first axis
/// @return the facet number along the first axis
int ImageParamsHelper::facetX() const
{
  ASKAPDEBUGASSERT(itsFacetX>=0);
  return itsFacetX;
}

/// @brief obtain the facet number along the second axis
/// @return the facet number along the second axis
int ImageParamsHelper::facetY() const
{
  ASKAPDEBUGASSERT(itsFacetY>=0);
  return itsFacetY;
}

/// @brief obtain the order of the Taylor term
/// @return the order of the Taylor term represented by this parameter
int ImageParamsHelper::order() const
{
  ASKAPDEBUGASSERT(isTaylorTerm());
  return itsOrder;
}


/// @brief obtain the full name of the image parameter
/// @details This method composes the full name of the parameter from 
/// the data stored internally. This returned full name should be the same 
/// as one passed in the parse method or in the constructor. This method can
/// be useful if this object is constructed directly without parsing a 
/// string and effectively represents a reverse operation.
std::string ImageParamsHelper::paramName() const
{
  ASKAPDEBUGASSERT(isValid());
  return itsName + suffix();                            
}

/// @brief obtain the full suffix
/// @details This method composes the suffix from the facet and taylor 
/// term. The result returned by paramName is just name()+suffix().
/// We need the suffix separately from the full name to be able to
/// index model in MSF simulations (the model may have a different name
/// than the field).
std::string ImageParamsHelper::suffix() const
{
  ASKAPDEBUGASSERT(isValid());
  std::string suffix;
  if (isTaylorTerm()) {
      suffix += ".taylor."+utility::toString<int>(itsOrder);
  }
  
  if (isFacet()) {
      suffix += ".facet."+utility::toString<int>(itsFacetX)+"."+
                               utility::toString<int>(itsFacetY);
  }
  return suffix;
}


/// @brief obtain the name of the image with just a facet suffix
/// @details To have MSMFS algorithm working with facets one needs
/// to be able to extract the name without suffix corresonding to
/// the Taylor decomposition, but with the faceting suffixes preserved,
/// if present. This method forms such a name.
/// @return the name of the parameter without Taylor suffixes
/// @note Facet suffix will not be added if the image is not a facet
std::string ImageParamsHelper::facetName() const
{
  ASKAPDEBUGASSERT(isValid());
  ImageParamsHelper temp(*this);
  temp.itsOrder = -1;
  return temp.paramName();
}

/// @brief obtain the name of the image with just a taylor suffix
/// @details To have MSMFS algorithm working with facets one needs
/// to be able to extract the name without suffixes corresonding to
/// facet, but with the taylor term suffixes preserved,
/// if present. This method forms such a name.
/// @return the name of the parameter without facet suffixes
/// @note Taylor suffix will not be added if the image is not a Taylor term
std::string ImageParamsHelper::taylorName() const
{
  ASKAPDEBUGASSERT(isValid());
  ImageParamsHelper temp(*this);
  temp.itsFacetX = -1;
  temp.itsFacetY = -1;
  return temp.paramName();  
}

/// @brief make this object a facet
/// @details It is sometimes necessary to merge faceting suffixes and
/// Taylor term suffix. This method makes the current image a facet with
/// given indices. 
/// @param[in] xFacet facet index along the first axis
/// @param[in] yFacet facet index along the second axis
void ImageParamsHelper::makeFacet(int xFacet, int yFacet)
{
  ASKAPDEBUGASSERT(isValid());
  ASKAPDEBUGASSERT(xFacet>=0);
  ASKAPDEBUGASSERT(yFacet>=0);
  itsFacetX = xFacet;
  itsFacetY = yFacet;
}

/// @brief make this object a facet
/// @details It is sometimes necessary to merge faceting suffixes and
/// Taylor term suffix. This method makes the current image a facet with
/// given indices. 
/// @param[in] order order in the Taylor series
void ImageParamsHelper::makeTaylorTerm(int order)
{
  ASKAPDEBUGASSERT(isValid());
  ASKAPDEBUGASSERT(order>=0);
  itsOrder = order;
}



} // namespace synthesis

} // namespace askap

