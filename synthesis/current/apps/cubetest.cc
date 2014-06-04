/// @file cubetest.cc
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

#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/Arrays/IPosition.h>
#include <casa/BasicSL/Complex.h>


#include <iostream>

using namespace casa;
using namespace std;

int main() 
{
  Cube<Float> cube(1,1,4,-1.);
  cout<<"initial cube shape: "<<cube.shape()<<endl;
  size_t index = 0;
  for (size_t row = 0; row<cube.nrow(); ++row) {
       for (size_t column = 0; column<cube.ncolumn();++column) {
            for (size_t plane = 0; plane<cube.nplane();++plane,++index) {
                 cube(row,column,plane) = Float(index);
            }
       }
  }
  cout<<cube.yzPlane(0).shape()<<endl;
  Matrix<Float> m = cube(IPosition(3,0,0,0),IPosition(3,0,0,3),IPosition(3,1,1,1)).nonDegenerate(IPosition(2,1,2)); //cube.yzPlane(0);
  cout<<m.shape()<<endl;
  cout<<m<<endl;
}

