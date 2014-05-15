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

