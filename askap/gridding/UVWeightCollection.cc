/// @file
/// @brief Collection of weights in the uv-domain
/// @details This class is related to traditional weighting. It represents a
/// set of grids with weight in the uv-domain. The extra dimension with arbitrary 
/// index is required for complex cases like joint imaging where (at some point) we
/// may want to build or apply weights depending on, say, beam index. 
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
#include <askap/gridding/UVWeightCollection.h>

// casa includes

// other 3rd party
#include <Blob/BlobArray.h>
#include <Blob/BlobSTL.h>


namespace askap {

namespace synthesis {

/// @brief add a new weight into collection
/// @details This method adds a new weight into collection setting up
/// an array with the given dimensions and populating it with zeros. 
/// If the weight already exists for the particular index, only cross-check of
/// dimensions is performed (so the method is safe to call multiple times, if desired)
/// @param[in] index integer index of the weight in the collection 
/// @param[in] uSize size in the direction of u-coordinate
/// @param[in] vSize size in the direction of v-coordinate
/// @param[in] nPlanes number of planes for the given weight (3rd dimension)
void UVWeightCollection::add(casacore::uInt index, casacore::uInt uSize, casacore::uInt vSize, casacore::uInt nPlanes)
{
   const std::map<casacore::uInt, casacore::Cube<float> >::const_iterator ci = itsData.find(index);
   if (ci != itsData.end()) {
       // just check the shape
       ASKAPASSERT(ci->second.nrow() == uSize);
       ASKAPASSERT(ci->second.ncolumn() == vSize);
       ASKAPASSERT(ci->second.nplane() == nPlanes);
   } else {
       // use explicit typing to trigger the right overload (and avoid implict copy or move)
       itsData.emplace(index, casacore::Cube<float>(uSize, vSize, nPlanes, 0.f));
   }
}

/// @brief add a new weight into collection or update existing weight
/// @details This version of add method works with a cube which is expected to contain some
/// weight information to be added. The behaviour depends on whether a weight exist for the given
/// index or not. If it does, the cube given in the parameters is simply added to the existing weight
/// cube. In this case, the passed cube is unchanged but must conform in shape. If the index doesn't 
/// exist, then a new element in the collection is setup from the cube passed as a parameter. And in 
/// this case, due to reference semantics of casacore arrays, both the parameter and the element 
/// in the collection would point to the same storage.
/// @param[in] index integer index of the weight in the collection 
/// @param[in] wt input weight cube
void UVWeightCollection::add(casacore::uInt index, const casacore::Cube<float> &wt)
{
   const std::map<casacore::uInt, casacore::Cube<float> >::iterator ci = itsData.find(index);
   if (ci != itsData.end()) {
       ci->second += wt;
   } else {
       itsData.emplace(index, wt);
   }
}

/// @brief obtain weight for the given index
/// @details This method returns the given item from the collection as a UVWeight object (restricting
/// the interface to the essential methods of the casacore Cube class). Note, an exception is thrown
/// if the given index doesn't exist.
/// @param[in] index integer index of the weight in the collection 
/// @return weight cube for the given index wrapped into a UVWeight object (reference semantics is used)
/// @note we may implement caching in the future, but at this stage it seems the access pattern via this 
/// method will be in the outer loop and then iterating within the same UVWeight object, so perhaps
/// no need to cache.
UVWeight UVWeightCollection::get(casacore::uInt index) const
{
   const std::map<casacore::uInt, casacore::Cube<float> >::const_iterator ci = itsData.find(index);
   ASKAPCHECK(ci != itsData.end(), "UVWeight with index "<<index<<" doesn't exist in the collection");
   // the following relies on C++17 if the object is non-copyable
   return UVWeight(ci->second);
}

/// @brief obtain weight for writing
/// @details This method returns a non-const reference to the actual cube object. If called in the 
/// read-only setting, it will be wrapped by UVWeight implicitly (as that object has the approproate 
/// constructor set up) making it an equivalent of the const method. This method is expected to be used
/// in the weight builder where we can benefit from casacore cube interface. It is hidden behind the interface
/// class to hinder breaking encapsulation in gridders.
/// @param[in] index integer index of the weight in the collection . Note, an exception is raised if the index
///                  doesn't exist in the collection
/// @return non-const reference to the weight cube for the given index
casacore::Cube<float>& UVWeightCollection::get(casacore::uInt index)
{
   const std::map<casacore::uInt, casacore::Cube<float> >::iterator ci = itsData.find(index);
   ASKAPCHECK(ci != itsData.end(), "UVWeight with index "<<index<<" doesn't exist in the collection");
   return ci->second;
}

/// @brief check that the index exists in the collection
/// @details We probably don't need this method, but add just in case
/// @param[in] index integer index of the weight in the collection 
/// @return true, if the weight with the given index is present (and, therefore, get method can be used for it)
bool UVWeightCollection::exists(casacore::uInt index) const
{
   const std::map<casacore::uInt, casacore::Cube<float> >::const_iterator ci = itsData.find(index);
   return ci != itsData.end();
}

/// @brief merge content from another collection
/// @details This method is equivalent to a set of get and add calls done for all indices present in the input collection.
/// Weights corresponding to new indices (i.e. not present in this class at the time of running this method) are added by
/// reference. In the typical use case of this method this would be equivalent to the ownership transfer of the particular 
/// weight cube.
/// @param[in] src input collection to merge from
void UVWeightCollection::merge(const UVWeightCollection &src)
{
   for (std::map<casacore::uInt, casacore::Cube<float> >::const_iterator ci = src.itsData.begin(); ci != src.itsData.end(); ++ci) {
        std::map<casacore::uInt, casacore::Cube<float> >::iterator match = itsData.find(ci->first);
        if (match != itsData.end()) {
            // check matching shape, although we could've left this to the array but the error message will be less readable
            ASKAPCHECK(match->second.shape() == ci->second.shape(), "UVWeight shape mismatch in an attempt to merge two collections, input shape = "
                       <<ci->second.shape()<<" existing element with index = "<<ci->first<<" has shape "<<match->second.shape());
            match->second += ci->second;
        } else {
            itsData.emplace(*ci);
        }
   }
}

/// @brief obtain a set of indices in this collection
/// @details As the indices are sparse we need a way to obtain a list of indices in the given collection. The C++ way of 
/// doing it is to have begin/end iterators. However, to do this without exposing internal structure (i.e. the map class) 
/// would require a specialised iterator type. For now just build and return a set of indices for simplicity although it
/// would imply more overhead. In practice, the number of indices we use is expected to be small, so it is hard to justify
/// extra complexity at this stage.
std::set<casacore::uInt> UVWeightCollection::indices() const
{
   std::set<casacore::uInt> result;
   for (std::map<casacore::uInt, casacore::Cube<float> >::const_iterator ci = itsData.begin(); ci != itsData.end(); ++ci) {
        result.insert(ci->first);
   }
   return result;
}

/// @brief write the object to a blob stream
/// @param[in] os the output stream
void UVWeightCollection::writeToBlob(LOFAR::BlobOStream& os) const
{
   os.putStart("UVWeightCollection", theirPayloadVersion); 
   os << itsData;
   os.putEnd(); 
}

/// @brief read the object from a blob stream
/// @param[in] is the input stream
/// @note Not sure whether the parameter should be made const or not
void UVWeightCollection::readFromBlob(LOFAR::BlobIStream& is)
{
   const int version = is.getStart("UVWeightCollection");
   ASKAPCHECK(version == theirPayloadVersion, "Payload version mismatch reading UVWeightCollection from blob");
   is >> itsData;
   is.getEnd();
}

} // namespace synthesis

} // namespace askap

