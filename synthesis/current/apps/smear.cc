#include <iostream>

#include <images/Images/PagedImage.h>
#include <images/Images/ImageInterface.h>
#include <casa/Arrays/Array.h>
#include <utils/PaddingUtils.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/Lattices/LatticeFFT.h>
#include <lattices/Lattices/LatticeExpr.h>
#include <askap/AskapError.h>

using namespace casa;
using namespace askap;

int main() {
try {
  PagedImage<Float> img("total.model");
  PagedImage<Float> out(img.shape(), img.coordinates(),"out.img");
  IPosition doubleShape(img.shape());
  ASKAPASSERT(doubleShape.nelements()>=2);
  doubleShape(0)*=2;
  doubleShape(1)*=2;
  ArrayLattice<casa::Complex> scratch(doubleShape);
  scratch.set(0.);
  scimath::PaddingUtils::inject(scratch,img);
  //scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(img)));
  LatticeFFT::cfft2d(scratch, True);

  ArrayLattice<casa::Complex> scratch2(scratch.shape());
  scratch2.set(0.);
  IPosition pos(scratch2.shape().nelements(),0);
  IPosition inpos(scratch.shape().nelements(),0);
  int nsteps = 10;
  for (int nx=0;nx<scratch.shape()[0];++nx) {
       for (int ny=0; ny<scratch.shape()[1];++ny) {
            pos(0)=nx; pos(1)=ny;
            Complex val(0.,0.);
            for (int steps=-nsteps;steps<=nsteps;++steps) {
                 double scale = 1.+ double(steps)/double(2*nsteps+1)*0.005;
                 inpos(0)=int(double(scratch.shape()[0]/2)+
                           (double(nx-scratch.shape()[0]/2)+0.5)*scale);
                 inpos(1)=int(double(scratch.shape()[1]/2)+
                           (double(ny-scratch.shape()[1]/2)+0.5)*scale);
                 if ((inpos(0)<0) || (inpos(1)<0) || (inpos(0)>=scratch.shape()(0)) || (inpos(1)>=scratch.shape()(1))) {
                    continue;
                 }
                 val+=scratch.getAt(inpos)/float(2*nsteps+1);
            }
            scratch2.putAt(val,pos);
       }
  }
 
  LatticeFFT::cfft2d(scratch2,False);
  //out.copyData(casa::LatticeExpr<casa::Float>(real(scratch2)));
  scimath::PaddingUtils::extract(out,scratch2);
}
catch(AskapError &ae) {
  std::cerr<<ae.what()<<std::endl;
}
catch(AipsError &ae) {
  std::cerr<<ae.what()<<std::endl;
}
catch(...) {
  std::cerr<<"Unexpected exception"<<std::endl;
}

  return 0;
}
