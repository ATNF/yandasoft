/// @file genoffsets.cc
///
/// test code to help with off-axis direction calculation
///
/// @copyright (c) 2014 CSIRO
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


#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".genoffsets");

#include <askap/AskapError.h>
#include <dataaccess/ParsetInterface.h>
#include <askap/AskapUtil.h>


// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/ArrayMath.h>

// LOFAR
#include <Common/ParameterSet.h>

// std
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

// boost
#include <boost/algorithm/string.hpp>

using namespace askap;
using namespace askap::accessors;

// the following code was copied from SynthesisParamsHelper. It probably needs to go to Base

/// @brief A helper method to parse string of quantities
/// @details Many parameters in parset file are given as quantities or
/// vectors of quantities, e.g. 8.0arcsec. This method allows
/// to parse a single string corresponding to such a parameter and return
/// a double value converted to the requested units.
/// @param[in] strval input string
/// @param[in] unit required units (given as a string)
/// @return converted value
double convertQuantity(const std::string &strval,
                      const std::string &unit)
{
    casa::Quantity q;
      
    casa::Quantity::read(q, strval);
    return q.getValue(casa::Unit(unit));
}
        
/// @brief A helper method to parse string of quantities
/// @details Many parameters in parset file are given as quantities or
/// vectors of quantities, e.g. [8.0arcsec,8.0arcsec]. This method allows
/// to parse vector of strings corresponding to such parameter and return
/// a vector of double values in the required units.
/// @param[in] strval input vector of strings
/// @param[in] unit required units (given as a string)
/// @return vector of doubles with converted values
std::vector<double> convertQuantity(const std::vector<std::string> &strval,
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
//    

class OffsetManager {
public:
    /// @return number of offsets
    size_t size() const;
    
    /// @brief add explicit offset
    /// @param[in] x true offset in degrees
    /// @param[in] y true offset in degrees
    void add(double x, double y);
    
    /// @brief add offset corresponding to two sources
    /// @details This method adds the offset of the source w.r.t. given reference direction 
    /// @param[in] dir source direction
    /// @param[in] ref reference direction
    /// @param[in] optional factor to scale that offset 
    void add(const casa::MVDirection &dir, const casa::MVDirection &ref, const double factor = 1.);

    /// @brief add an offset by forming a triangle with two existing ones 
    /// @details helper method to add extra offset to the list of offsets to form an equilateral triangle with two 
    /// specified points. There are two possible such points selected by the 3rd parameter (raAdd true and false)
    /// @param[in] pt1 index of the first vertex
    /// @param[in] pt2 index of the second vertex
    /// @param[in] raAdd roughly speaking left or right as there are two possible triangles
    void complementTriangle(size_t pt1, size_t pt2, bool raAdd);

    /// @brief access to the stored offsets
    /// @param[in] index index of the offset
    /// @return pair of offsets in degrees
    std::pair<double, double> operator[](size_t index) const; 
    
private:
    /// ra offsets in degrees
    std::vector<double> itsXOffsets;
    /// dec offsets in degrees
    std::vector<double> itsYOffsets;
};

/// @brief access to the stored offsets
/// @param[in] index index of the offset
/// @return pair of offsets in degrees
std::pair<double, double> OffsetManager::operator[](size_t index) const
{
   ASKAPDEBUGASSERT(itsXOffsets.size() == itsYOffsets.size());
   ASKAPASSERT(index < itsXOffsets.size());
   return std::pair<double,double>(itsXOffsets[index], itsYOffsets[index]);
} 


/// @return number of offsets
size_t OffsetManager::size() const {
   ASKAPDEBUGASSERT(itsXOffsets.size() == itsYOffsets.size());
   return itsXOffsets.size();
};

/// @brief add explicit offset
/// @param[in] x true offset in ra in degrees
/// @param[in] y true offset in dec in degrees
void OffsetManager::add(double x, double y)
{
   itsXOffsets.push_back(x);
   itsYOffsets.push_back(y);
}


/// @brief add an offset by forming a triangle with two existing ones 
/// @details helper method to add extra offset to the list of offsets to form an equilateral triangle with two 
/// specified points. There are two possible such points selected by the 3rd parameter (raAdd true and false)
/// @param[in] pt1 index of the first vertex
/// @param[in] pt2 index of the second vertex
/// @param[in] raAdd roughly speaking left or right as there are two possible triangles
void OffsetManager::complementTriangle(size_t pt1, size_t pt2, bool raAdd)
{
   ASKAPDEBUGASSERT(pt1 < itsXOffsets.size());
   ASKAPDEBUGASSERT(pt2 < itsXOffsets.size());
   ASKAPDEBUGASSERT(itsXOffsets.size() == itsYOffsets.size());
   const double xNew = itsXOffsets[pt1] + 0.5 * (itsXOffsets[pt2] - itsXOffsets[pt1]) + (raAdd ? +1. : -1.) * 
                       (itsYOffsets[pt2] - itsYOffsets[pt1]);
   const double yNew = itsYOffsets[pt1] + 0.5 * (itsYOffsets[pt2] - itsYOffsets[pt1]) + (raAdd ? -1. : +1.) * 
                       (itsXOffsets[pt2] - itsXOffsets[pt1]);
   itsXOffsets.push_back(xNew);
   itsYOffsets.push_back(yNew);
}

/// @brief add offset corresponding to two sources
/// @details This method adds the offset of the source w.r.t. given reference direction 
/// @param[in] dir source direction
/// @param[in] ref reference direction
/// @param[in] optional factor to scale that offset 
void OffsetManager::add(const casa::MVDirection &dir, const casa::MVDirection &ref, const double factor)
{
  const double offset1 = sin(dir.getLong() - ref.getLong()) * cos(dir.getLat());
  const double offset2 = sin(dir.getLat()) * cos(ref.getLat()) - cos(dir.getLat()) * sin(ref.getLat())
                                                  * cos(dir.getLong() - ref.getLong());
  add(factor * offset1 / casa::C::pi * 180., factor * offset2 / casa::C::pi * 180.);
}

/// @brief helper method to extract direction string
/// @param[in] dirStr vector of strings representing direction
/// @return direction
casa::MVDirection  extractDir(const std::vector<std::string> &dirStr) {
   ASKAPCHECK(dirStr.size() >= 2, "Expect at least two elements in the direction string; you have="<<dirStr);
   const casa::MVDirection dir(convertQuantity(dirStr[0],"rad"), convertQuantity(dirStr[1],"rad"));
   return dir;
}
 
/// @brief extract index
/// @param[in] vec vector of strings
/// @param[in] val value
/// @return index
size_t extractIndex(const std::vector<std::string> &vec, const std::string &val) 
{
  const std::vector<std::string>::const_iterator ci = std::find(vec.begin(), vec.end(),val);
  ASKAPCHECK(ci != vec.end(), "Source "<<val<<" is not defined in the list");
  return ci - vec.begin();
} 

void makeOffsets(const LOFAR::ParameterSet &parset) {
  std::vector<casa::MVDirection> directions;
  const std::vector<std::string> srcNames = parset.getStringVector("sources",std::vector<std::string>());
  for (std::vector<std::string>::const_iterator ci = srcNames.begin(); ci != srcNames.end(); ++ci) {
       directions.push_back(extractDir(parset.getStringVector("sources." + *ci)));                      
  }
  const size_t nBeams = parset.getInt32("nbeams");
  OffsetManager ofm;
  for (size_t beam=0; beam<nBeams; ++beam) {
       const std::string beamKey = "beam"+utility::toString<size_t>(beam+1);
       ASKAPCHECK(parset.isDefined(beamKey), "Missing definition for beam "<<beam + 1);
       const std::string descr = boost::trim_copy(parset.getString(beamKey));
       //std::cout<<descr<<std::endl;
       
       const size_t pos1 = descr.find("(");
       ASKAPCHECK(pos1 != std::string::npos, "Unable to parse: "<<descr);
       const size_t pos2 = descr.find(")",pos1);
       ASKAPCHECK(pos2 != std::string::npos, "Unable to parse: "<<descr);
       const std::string operation = descr.substr(0,pos1);       
       std::vector<std::string> parts;
       std::string operands = descr.substr(pos1+1,pos2-pos1-1);
       boost::split(parts, operands, boost::is_any_of(","));
       if (operation.find("triangle") != std::string::npos) {
           ASKAPCHECK(parts.size() == 3, "Expect exactly 3 elements in the triangle description");
           const int index1 = utility::fromString<int>(boost::trim_copy(parts[0]));
           const int index2 = utility::fromString<int>(boost::trim_copy(parts[1]));
           ASKAPCHECK((index1 > 0) && (index1 <= static_cast<int>(ofm.size())),  
                       "First beam index of the triangle description should be " 
                       "positive and should relate only to beams which are already defined");
           ASKAPCHECK((index2 > 0) && (index2 <= static_cast<int>(ofm.size())),  
                      "Second beam index of the triangle description should be "
                      "positive and should relate only to beams which are already defined");
           const std::string mode = boost::trim_copy(parts[2]);
           ASKAPCHECK((mode == "left") || (mode == "right"), 
                      "Triangle should be either left or right (3rd parameter), you have "<<mode);
                                              
           ofm.complementTriangle(static_cast<size_t>(index1) - 1, static_cast<size_t>(index2) - 1, mode == "left");
       } else {
           ASKAPCHECK(parts.size() == 2, "Expect exactly 2 elements in the simple offset description");
           const std::string param1 = boost::trim_copy(parts[0]);
           const std::string param2 = boost::trim_copy(parts[1]);
           if (operation.find("offset") != std::string::npos) {
               ofm.add(utility::fromString<double>(param1), utility::fromString<double>(param2));
           } else {
               const double factor = (pos1 == 0) ? 1. : utility::fromString<double>(operation);
               const size_t index1  = extractIndex(srcNames, param1);
               const size_t index2  = extractIndex(srcNames, param2);
               ASKAPDEBUGASSERT(index1 < directions.size());
               ASKAPDEBUGASSERT(index2 < directions.size());
                   
               ofm.add(directions[index1],directions[index2],factor);
           }
       }
  }
  // reference position
  // default - Sun at MRO
  casa::MPosition mroPos(casa::MVPosition(casa::Quantity(370.81, "m"), casa::Quantity(116.6310372795, "deg"), 
                          casa::Quantity(-26.6991531922, "deg")), casa::MPosition::Ref(casa::MPosition::WGS84));
  casa::Quantity q;
  ASKAPCHECK(casa::MVTime::read(q, "today"), "MVTime::read failed");
  std::cout<<"Current UT MJD: "<<q<<std::endl;
  const casa::MEpoch::Ref utcRef(casa::MEpoch::UTC);
  casa::MEpoch when(casa::MVEpoch(q), utcRef);
  casa::MeasFrame frame(mroPos, when);
  casa::MVDirection refDir = casa::MDirection::Convert(casa::MDirection(casa::MDirection::SUN), casa::MDirection::Ref(casa::MDirection::J2000,frame))().getValue();
  if (parset.isDefined("refdir")) {
      refDir = extractDir(parset.getStringVector("refdir"));
  }
  // finally doing the output
  const double factor = -1;
  
  for (size_t i=0; i<ofm.size(); ++i) {
       std::pair<double,double> offsetsInDeg = ofm[i];
       const double offset1 = offsetsInDeg.first / 180. * casa::C::pi;
       const double offset2 = offsetsInDeg.second / 180. * casa::C::pi;
       casa::MVDirection testDir = refDir;
       testDir.shift(offset1*factor,offset2*factor, casa::True);
       std::cout<<"offset ("<<offsetsInDeg.first<<","<<offsetsInDeg.second<<") applied to test direction: "<<printDirection(testDir)<<std::endl;
  }
  
};

int main(int argc, char **argv) {
  try {
     if (argc!=2) {
         std::cerr<<"Usage: "<<argv[0]<<" parset"<<std::endl;
	 return -2;
     }

     casa::Timer timer;

     timer.mark();
     std::cerr<<"Initialization: "<<timer.real()<<std::endl;
     timer.mark();
     makeOffsets(LOFAR::ParameterSet(argv[1]));
     std::cerr<<"Job: "<<timer.real()<<std::endl;
     
  }
  catch(const AskapError &ce) {
     std::cerr<<"AskapError has been caught. "<<ce.what()<<std::endl;
     return -1;
  }
  catch(const std::exception &ex) {
     std::cerr<<"std::exception has been caught. "<<ex.what()<<std::endl;
     return -1;
  }
  catch(...) {
     std::cerr<<"An unexpected exception has been caught"<<std::endl;
     return -1;
  }
  return 0;
}
