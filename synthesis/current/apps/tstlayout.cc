///
/// @file
///
/// This program analyses antenna layout (written to investigate snap-shot imaging limitations)
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

// Std includes
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

// Askapsoft includes
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, "");

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <CommandLineParser.h>
#include <simulation/Simulator.h>
#include <askapparallel/MPIComms.h>
#include <casa/Quanta/MVDirection.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Matrix.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <measures/Measures/UVWMachine.h>
#include <measures/Measures/MPosition.h>
#include <images/Images/PagedImage.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <utils/EigenDecompose.h>
#include <measurementequation/SynthesisParamsHelper.h>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::askapparallel;
using namespace LOFAR;

/// @brief load layout
/// @details
/// @param[in] fname input file name
/// @param[out] x coordinate X
/// @param[out] y coordinate Y
/// @param[out] z coordinate Z
void loadLayout(const std::string &fname, casa::Vector<double> &x, casa::Vector<double> &y, casa::Vector<double> &z)
{
  ParameterSet parset(fname);
  
  const std::string telName = parset.getString("antennas.telescope");
  ASKAPLOG_INFO_STR(logger, "Loading " << telName);
  std::ostringstream oos;
  oos << "antennas." << telName << ".";
  ParameterSet antParset(parset.makeSubset(oos.str()));
  
  ASKAPCHECK(antParset.isDefined("names"), "Subset (antennas."<<telName<<") of the antenna definition parset does not have 'names' keyword.");
  std::vector<string> antNames(antParset.getStringVector("names"));
  const int nAnt = antNames.size();
  ASKAPCHECK(nAnt > 0, "No antennas defined in parset file");
  
  std::string coordinates = antParset.getString("coordinates", "local");
  ASKAPCHECK((coordinates == "local") || (coordinates == "global"), "Coordinates type unknown");

  /// Csimulator.ASKAP.scale=0.333
  const float scale = antParset.getFloat("scale", 1.0);

  /// Now we get the coordinates for each antenna in turn
  casa::Vector<double> xBuf(nAnt);
  casa::Vector<double> yBuf(nAnt);
  casa::Vector<double> zBuf(nAnt);
  
  /// Antenna information in the form:
  /// antennas.ASKAP.antenna0=[x,y,z]
  /// ...
  for (int iant = 0; iant < nAnt; ++iant) {
       std::vector<float> xyz = antParset.getFloatVector(antNames[iant]);
       ASKAPCHECK(xyz.size()>=3, "Error loading ant="<<iant+1<<", xyz.size()="<<xyz.size());
       xBuf[iant] = xyz[0] * scale;
       yBuf[iant] = xyz[1] * scale;
       zBuf[iant] = xyz[2] * scale;
  }

  /// Csimulator.ASKAP.location=[+115deg, -26deg, 192km, WGS84]
  if((coordinates == "local") || (coordinates == "longlat") ) {
     const casa::MPosition location = asMPosition(antParset.getStringVector("location"));
      
     casa::MVAngle mvLong = location.getAngle().getValue()(0);
     casa::MVAngle mvLat = location.getAngle().getValue()(1);

     ASKAPLOG_INFO_STR(logger, "Using local coordinates for the antennas: Reference position = "
                              << mvLong.string(casa::MVAngle::ANGLE, 7) << " "
                              << mvLat.string(casa::MVAngle::DIG2, 7));
     x.resize(xBuf.nelements());                              
     y.resize(yBuf.nelements());                              
     z.resize(zBuf.nelements());                            
     if (coordinates == "local") {  
         Simulator::local2global(x, y, z, location, xBuf, yBuf, zBuf);      
     } else {
      Simulator::longlat2global(x, y, z, location, xBuf, yBuf, zBuf);
     }
  } else if (coordinates == "global") {
      x.assign(xBuf);
      y.assign(yBuf);
      z.assign(zBuf);
      ASKAPLOG_INFO_STR(logger, "Using global coordinates for the antennas");
  } else {
        ASKAPLOG_INFO_STR(logger, "Unknown coordinate system type: " << coordinates);
  }

  ASKAPLOG_INFO_STR(logger, "Successfully defined " << nAnt
          << " antennas of " << telName);
  
}

/// @brief read layout, form baselines
/// @details
/// @param[in] fname input file name
/// @param[out] x coordinate X
/// @param[out] y coordinate Y
/// @param[out] z coordinate Z
void getBaselines(const std::string &fname, casa::Vector<double> &x, casa::Vector<double> &y, casa::Vector<double> &z)
{
  casa::Vector<double> xAnt,yAnt,zAnt;
  loadLayout(fname,xAnt,yAnt,zAnt);
  ASKAPCHECK((xAnt.nelements() == yAnt.nelements()) && (xAnt.nelements() == zAnt.nelements()), 
             "Expect the same number of elements in xAnt, yAnt, zAnt");
  /*
  // to export the layout in ASCII format
  for (casa::uInt ant=0; ant<xAnt.nelements(); ++ant) {
       std::cout<<std::scientific<<std::setprecision(15)<<xAnt[ant]<<" "<<yAnt[ant]<<" "<<zAnt[ant]<<std::endl;
  }
  throw 1;
  //
  */         
  const casa::uInt nBaselines = (xAnt.nelements()*(xAnt.nelements()-1)) / 2;
  x.resize(nBaselines);
  y.resize(nBaselines);
  z.resize(nBaselines);
  for (casa::uInt ant1=0,baseline=0;ant1<xAnt.nelements(); ++ant1) {
       for (casa::uInt ant2=0;ant2<ant1; ++ant2,++baseline) {
            ASKAPASSERT(baseline < nBaselines);
            x[baseline] = xAnt[ant1] - xAnt[ant2];
            y[baseline] = yAnt[ant1] - yAnt[ant2];
            z[baseline] = zAnt[ant1] - zAnt[ant2];            
       }
  }
  ASKAPLOG_INFO_STR(logger, "Formed "<<nBaselines<<" baselines");
}

/// @brief obtain uvws
/// @details for given declination and hour angle
/// @param[in] x baseline coordinates X
/// @param[in] y baseline coordinates Y
/// @param[in] z baseline coordinates Z
/// @param[in] dec declination of the tangent point (in radians)
/// @param[in] H0 hour angle of the tanget point (at longitude of 0, in radians)
/// @param[out] u baseline coordinate u
/// @param[out] v baseline coordinate v
/// @param[out] w baseline coordinate w
void calculateUVW(const casa::Vector<double> &x, const casa::Vector<double> &y, const casa::Vector<double> &z,
                  double dec, double H0, casa::Vector<double> &u, casa::Vector<double> &v,casa::Vector<double> &w)
{
   ASKAPDEBUGASSERT(x.nelements() == y.nelements());
   ASKAPDEBUGASSERT(x.nelements() == z.nelements());
   ASKAPDEBUGASSERT(y.nelements() == z.nelements());
   const casa::uInt size = x.nelements();
   ASKAPDEBUGASSERT(size != 0);
   u.resize(size);
   v.resize(size);
   w.resize(size);
   const double sDec = sin(dec);
   const double cDec = cos(dec);
   const double sH0 = sin(H0);
   const double cH0 = cos(H0);
   for (casa::uInt row = 0; row<size; ++row) {
        u[row] = sH0 * x[row] + cH0 * y[row];
        v[row] = -sDec * cH0 * x[row] + sDec * sH0 * y[row] + cDec * z[row];
        w[row] = cDec * cH0 * x[row] - cDec * sH0 * y[row] + sDec * z[row]; 
   }
}                  

/// @brief analyse the layout
/// @details
/// @param[in] x baseline coordinate X
/// @param[in] y baseline coordinate Y
/// @param[in] z baseline coordinate Z
/// @return vector normal to the layout
casa::Vector<double> analyseBaselines(casa::Vector<double> &x, casa::Vector<double> &y, casa::Vector<double> &z)
{
  casa::Matrix<double> normalMatr(3,3,0.);
  const casa::uInt nBaselines = x.nelements();
  for (casa::uInt b = 0; b<nBaselines; ++b) {
       normalMatr(0,0) += casa::square(x[b]);
       normalMatr(1,1) += casa::square(y[b]);
       normalMatr(2,2) += casa::square(z[b]);
       normalMatr(0,1) += x[b]*y[b];
       normalMatr(0,2) += x[b]*z[b];
       normalMatr(1,2) += y[b]*z[b];
  }
  normalMatr(1,0) = normalMatr(0,1);
  normalMatr(2,0) = normalMatr(0,2);
  normalMatr(2,1) = normalMatr(1,2);
  
  casa::Vector<double> eVal;
  casa::Matrix<double> eVect;
  scimath::symEigenDecompose(normalMatr,eVal,eVect);
 
  ASKAPLOG_INFO_STR(logger, "Normal matrix: "<<normalMatr);
  ASKAPLOG_INFO_STR(logger, "eVal: "<<eVal);
  ASKAPLOG_INFO_STR(logger, "eVect: "<<eVect);
  casa::Vector<double> normalVector(eVect.column(2).copy());
  const double norm = casa::sum(casa::square(normalVector));
  normalVector /= norm;
  ASKAPLOG_INFO_STR(logger, "Normalised vector normal to the best fit plane: "<<normalVector);
  double maxDeviation = -1;
  for (casa::uInt b = 0; b<nBaselines; ++b) {
       // (baseline,normalVector)
       const double dotProduct = x[b]*normalVector[0] + y[b]*normalVector[1] + z[b]*normalVector[2];
       if (fabs(dotProduct) > maxDeviation) {
           maxDeviation = fabs(dotProduct);
       }
  }
  ASKAPLOG_INFO_STR(logger, "Largest deviation from the plane is "<<maxDeviation<<" metres");
  return normalVector;
}

/// @brief analyse the uvw
/// @details
/// @param[in] u baseline coordinates U
/// @param[in] v baseline coordinates V
/// @param[in] w baseline coordinates W
/// @param[out] nv normalised vector orthogonal to the fitted plane in uvw-coordinates
/// @param[in] verbose if true, debug messages will be written into log
/// @return largest residual w-term (negative value means the fit has failed)
double analyseUVW(casa::Vector<double> &u, casa::Vector<double> &v, casa::Vector<double> &w, 
                  casa::Vector<double> &nv, bool verbose = false)
{
  nv.resize(3);
  casa::Matrix<double> normalMatr(3,3,0.);
  const casa::uInt nBaselines = u.nelements();
  for (casa::uInt b = 0; b<nBaselines; ++b) {
       normalMatr(0,0) += casa::square(u[b]);
       normalMatr(1,1) += casa::square(v[b]);
       normalMatr(2,2) += casa::square(w[b]);
       normalMatr(0,1) += u[b]*v[b];
       normalMatr(0,2) += u[b]*w[b];
       normalMatr(1,2) += v[b]*w[b];
  }
  normalMatr(1,0) = normalMatr(0,1);
  normalMatr(2,0) = normalMatr(0,2);
  normalMatr(2,1) = normalMatr(1,2);
  
  casa::Vector<double> eVal;
  casa::Matrix<double> eVect;
  scimath::symEigenDecompose(normalMatr,eVal,eVect);
  
  if (verbose) {
      ASKAPLOG_DEBUG_STR(logger, "(uvw) eVal: "<<eVal);
  }
  casa::Vector<double> normalVector(eVect.column(2).copy());
  const double norm = casa::sum(casa::square(normalVector));
  normalVector /= norm;
  if (verbose) {
      ASKAPLOG_DEBUG_STR(logger, "Normalised vector normal to the best fit uvw plane: "<<normalVector);
  }
  nv = normalVector;
  if (fabs(normalVector[2]) > 1e-6) {
      normalVector[0] /= normalVector[2];
      normalVector[1] /= normalVector[2];
      normalVector[0] *= -1.;
      normalVector[1] *= -1.;
      if (verbose) {      
          ASKAPLOG_DEBUG_STR(logger, "Best fit plane w = u * "<<normalVector[0]<<" + v * "<<normalVector[1]);
      }
      
      double maxDeviation = -1;
      for (casa::uInt b = 0; b<nBaselines; ++b) {
          // residual w-term
          const double resW = w[b] - u[b]*normalVector[0] - v[b] * normalVector[1];
          if (fabs(resW) > maxDeviation) {
              maxDeviation = fabs(resW);
          }
      }
      if (verbose) {
          ASKAPLOG_DEBUG_STR(logger, "Largest residual w-term  is "<<maxDeviation<<" metres");
      }
      return maxDeviation;      
  } else {
     if (verbose) {
         ASKAPLOG_DEBUG_STR(logger, "w is independent on u and v in this layout. Fitting failed");
     }
  }
  return -1.;
}

/// @brief map residual w-term
/// @details This method calls calculateUVW/analyseUVW for different dec and H0 and makes an image
/// showing the behavior of the residual w-term.
/// @param[in] name name of the output image
/// @param[in] longitude reference longiude (i.e. where hour angle is zero; in radians)
/// @param[in] x baseline coordinate X
/// @param[in] y baseline coordinate Y
/// @param[in] z baseline coordinate Z
void mapResidualWTerm(const std::string &name, const double longitude, 
            const casa::Vector<double> &x, const casa::Vector<double> &y, const casa::Vector<double> &z)
{
  // sizes
  const int xSize = 512;
  const int ySize = 512;
  // start/stop values
  const double startDec = -casa::C::pi/2.;
  const double stopDec = 0.;
  const double incH = 2.*casa::C::pi/double(xSize)*2./3.;
  const double incDec = -(stopDec - startDec)/double(ySize);
  
  // buffer
  casa::Matrix<float> buf(xSize,ySize,0.);

  // coordinate system
  casa::CoordinateSystem coords;
  casa::Matrix<double> xform(2,2,0.);
  xform.diagonal() = 1.;
  const casa::DirectionCoordinate dc(casa::MDirection::HADEC, 
                  casa::Projection(casa::Projection::STG), casa::C::pi+incH , startDec,
                  incH,incDec,xform,xSize/2,1.);  
  coords.addCoordinate(dc);
  // filling the buffer
  casa::Vector<casa::Double> pixel(2,0.);
  casa::MVDirection world;
  casa::Vector<double> u,v,w;
  casa::Vector<double> unused;
  for (int xx = 0; xx<xSize; ++xx) {
       for (int yy=0; yy<ySize; ++yy) {
            pixel[0] = double(xx);
            pixel[1] = double(yy);
            if (dc.toWorld(world,pixel)) {
                calculateUVW(x,y,z,world.getLat(),world.getLong()-longitude,u,v,w);
                const double resW = analyseUVW(u,v,w,unused);                
                if (resW>=0) {
                    buf(xx,yy) = resW;
                }
            }
       }
  }
  
  // writing the actual image
  casa::PagedImage<casa::Float> result(casa::TiledShape(buf.shape()), coords, name);
  casa::ArrayLattice<casa::Float> lattice(buf);
  result.copyData(lattice);  
}            

/// @brief test uvw-rotation
/// @details this method computes uvw's for two directions on the sky using the first
/// principles and then does uvw-rotation from one direction to another and compares two
/// sets of uvws
/// @param[in] x baseline coordinate X
/// @param[in] y baseline coordinate Y
/// @param[in] z baseline coordinate Z
void testUVWRotation(const casa::Vector<double> &x, const casa::Vector<double> &y, const casa::Vector<double> &z)
{
  ASKAPLOG_INFO_STR(logger, "Test of the uvw-rotation");
  const casa::MVDirection tangent(SynthesisParamsHelper::convertQuantity("12h30m00.000","rad"),
                                  SynthesisParamsHelper::convertQuantity("-45.00.00.000","rad"));
  const casa::MDirection tangentDir(tangent, casa::MDirection::J2000);
  const double raOffset = 0./180.*casa::C::pi;
  const double decOffset = 2./180.*casa::C::pi;
  casa::MDirection offsetDir(tangentDir);
  offsetDir.shift(-raOffset,decOffset,casa::True);
  casa::Vector<double> uOffset,vOffset,wOffset;
  // uvw at longitude where the region is close to transit
  calculateUVW(x,y,z,offsetDir.getValue().getLat(),casa::C::pi-offsetDir.getValue().getLong(),uOffset,vOffset,wOffset);
  casa::Vector<double> uTangent,vTangent,wTangent;
  calculateUVW(x,y,z,tangentDir.getValue().getLat(),casa::C::pi-tangentDir.getValue().getLong(),uTangent,vTangent,wTangent);
  casa::UVWMachine machine(offsetDir,tangentDir,false, true);
  //casa::UVWMachine machine(tangentDir,offsetDir,false, true);
  for (casa::uInt row=0; row<x.nelements(); ++row) {
       casa::Vector<double> buf(3);
       buf[0]=uOffset[row];
       buf[1]=vOffset[row];
       buf[2]=wOffset[row];
       double delay;
       machine.convertUVW(delay, buf);
       if (row%50 == 0) {
           ASKAPLOG_INFO_STR(logger,"row="<<row<<" uvw_rotated-uvw_original: "<<buf[0]-uTangent[row]<<" "<<buf[1]-vTangent[row]<<" "<<
                                 buf[2]-wTangent[row]);
       }
       uOffset[row] = buf[0];
       vOffset[row] = buf[1];
       wOffset[row] = buf[2];       
  }
  casa::Vector<double> nvTangent, nvOffset;
  const double wMaxTangent = analyseUVW(uTangent,vTangent,wTangent,nvTangent);
  const double wMaxOffset = analyseUVW(uOffset,vOffset,wOffset,nvOffset);
  ASKAPLOG_INFO_STR(logger,"Largest w-term deviation: "<<wMaxTangent<<" metres for uvw's calculated for original direction");
  ASKAPLOG_INFO_STR(logger,"  and: "<<wMaxOffset<<" metres for uvw's calculated for offset direction and rotated back");
  const double dotProduct = nvTangent[0]*nvOffset[0] + nvTangent[1]*nvOffset[1] + nvTangent[2]*nvOffset[2];
  ASKAPLOG_INFO_STR(logger,"The angle between best fit planes for original uvw distribution and the distribution calculated for the direction, which");
  ASKAPLOG_INFO_STR(logger,"is offset "<< raOffset/casa::C::pi*180.<<" deg in ra and "<<decOffset/casa::C::pi*180.<<
            " deg in dec and rotated back to the original direction, is "<< acos(dotProduct)/casa::C::pi*180.<<" deg");
}

/// @brief main
/// @param[in] argc number of arguments
/// @param[in] argv vector of arguments
/// @retirm status
int main(int argc, char **argv) {
  try {
     cmdlineparser::Parser parser; // a command line parser
     
     // command line parameter
     cmdlineparser::GenericParameter<std::string> cfgName;
     parser.add(cfgName, cmdlineparser::Parser::throw_exception);
     cmdlineparser::GenericParameter<std::string> imgName;
     parser.add(imgName,cmdlineparser::Parser::return_default);
     parser.process(argc, argv);
     
     // Initialize MPI (also succeeds if no MPI available).
     askap::askapparallel::MPIComms comms(argc, argv);
     ASKAPLOG_INIT("askap.log_cfg");
     
     casa::Vector<double> x,y,z;
     getBaselines(cfgName,x,y,z);
     casa::Vector<double> normalVector = analyseBaselines(x,y,z);
     ASKAPDEBUGASSERT(normalVector.nelements() == 3);
     const double centreLong = atan2(normalVector[1],normalVector[0]);
     const double centreLat = asin(normalVector[2]);
     ASKAPLOG_INFO_STR(logger, "Centre of the layout is at longitude "<<centreLong*180/casa::C::pi<<
                       " and latitude "<<centreLat*180/casa::C::pi);

     testUVWRotation(x,y,z);

     if (imgName.defined()) { 
         ASKAPLOG_INFO_STR(logger, "Image showing largest residual w-term vs. local hour angle and declination will be stored in "<<
                                    imgName.getValue());
         mapResidualWTerm(imgName,centreLong,x,y,z);                           
     }
  }
  catch(const cmdlineparser::XParser &) {
     std::cerr<<"Usage "<<argv[0]<<" cfg_name [output_image]"<<std::endl;
	 return -2;    
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

  
