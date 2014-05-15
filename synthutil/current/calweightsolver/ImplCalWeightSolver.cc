/// @file
/// 
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>



// ImplCalWeightSolver -  implementation of the algorithm which
// solves for the best FPA weights for an optimum calibration on a
// given sky model
//
//

#include <calweightsolver/ImplCalWeightSolver.h>
#include <components/ComponentModels/Flux.h>
#include <components/ComponentModels/SpectralModel.h>
#include <components/ComponentModels/ComponentShape.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MeasConvert.h>
#include <casa/OS/Path.h>
#include <casa/BasicMath/Math.h>
#include <casa/Utilities/Assert.h>
#include <calweightsolver/EigenSolve.h>
#include <calweightsolver/SkyCatalogTabWriter.h>
#include <iomanip>

#include <askap/AskapError.h>

using namespace casa;
using namespace std;

using namespace askap;

// ImplCalWeightSolver

ImplCalWeightSolver::ImplCalWeightSolver() throw() :
	 vp_real(NULL),vp_imag(NULL) {}

ImplCalWeightSolver::~ImplCalWeightSolver() throw(casa::AipsError) 
{
  if (vp_real!=NULL) delete vp_real;
  if (vp_imag!=NULL) delete vp_imag;
}

// set up calculation for a given pointing centre and sky model
void ImplCalWeightSolver::setSky(const casa::MDirection &ipc,
       const casa::String &clname) throw(casa::AipsError)
{
  cl=ComponentList(Path(clname),True);
  pc=ipc;
}

// set up the voltage pattern from a disk-based image
void ImplCalWeightSolver::setVP(const casa::String &namer,
		                const casa::String &namei) 
	 throw(casa::AipsError)
{
  if (vp_real!=NULL) delete vp_real;
  if (vp_imag!=NULL) delete vp_imag;
  vp_real=new PagedImage<Float>(namer);
  vp_imag=new PagedImage<Float>(namei);
  if (vp_real->shape()!=vp_imag->shape())
      throw AipsError("The shape of real and imagingary parts should be identical");
  if (vp_real->shape().nelements()<2) 
      throw AipsError("VP image should have at least 2 directional axes");
}

/// @brief make synthetic beam
/// @details This method constructs synthetic primary beam for the given weights.
/// @param[in] name output image name
/// @param[in] weights vector of weights
void ImplCalWeightSolver::makeSyntheticPB(const std::string &name, 
	                     const casa::Vector<casa::Complex> &weights)
{
  ASKAPASSERT(vp_real!=NULL);
  ASKAPASSERT(vp_imag!=NULL);
  const casa::CoordinateSystem cs = vp_real->coordinates();
  PagedImage<float> result(name);
  result.setCoordinateInfo(cs);
  
}

// calculate visibility matrix for given feed_offsets
// uvw - a vector with the uvw coordinates (in the units of wavelength)
void ImplCalWeightSolver::formVisMatrix(const Matrix<Double> &feed_offsets,
             const Vector<Double> &uvw) const throw(casa::AipsError)
{
  if (uvw.nelements()!=3)
     throw AipsError("uvw should be a vector with 3 elements");

  if (vp_real==NULL) 
     throw AipsError("A vp image should be set before calling formVisMatrix");

  if (vp_imag==NULL) 
     throw AipsError("A vp image should be set before calling formVisMatrix");

  if (!cl.nelements()) 
      throw AipsError("At least one sky-component should be set");

  // assume that the number of measurements is equivalent to the number
  // of feeds. More sets of weights will not help as the data will not
  // be linearly independent to the first order
  vismatrix.resize(feed_offsets.nrow(),feed_offsets.nrow());
  vismatrix=0.;

  for (uInt comp=0;comp<cl.nelements();++comp) {
       const SkyComponent &skc=cl.component(comp);
       if (skc.spectrum().type()!=ComponentType::CONSTANT_SPECTRUM)
	   throw AipsError("Non-constant spectrum sources are not supported");
       // we have to convert reference system of the sky component to
       // that used for the phase/pointing centre
       const MDirection &srcdir_ownref=skc.shape().refDirection();
       // srcdir -> direction in the same reference system as used
       // for the phase/pointing centre
       MVDirection srcdir=MDirection::Convert(srcdir_ownref.getRef(),
		       pc.getRef())(srcdir_ownref.getValue()).getValue();
       // true offsets w.r.t the dish centre
       Double l0=sin(srcdir.getLong()-
                   pc.getValue().getLong())*cos(srcdir.getLat());
       Double m0=sin(srcdir.getLat())*cos(pc.getValue().getLat())-
		 cos(srcdir.getLat())*sin(pc.getValue().getLat())*
	          cos(srcdir.getLong()-pc.getValue().getLong());
       // component flux
       Flux<Double> flux=skc.flux();
       flux.convertPol(ComponentType::STOKES);
       flux.convertUnit("Jy");
	   
       for (uInt feed=0;feed<vismatrix.ncolumn();++feed) {
	    Complex visbuf(1.,0.);
	    // following formulae are approximate. We need to replace 
	    // them with exact formulae for the angular offset of 
	    // the source w.r.t. the appropriate feed centre
	    Double l=l0-feed_offsets(feed,0);
	    Double m=m0-feed_offsets(feed,1);
	    if (!getVPValue(visbuf,l,m)) continue;

            
	    /*  // Full interferometric visibility
	    Double phasor=2*M_PI*(uvw[0]*l+uvw[1]*m+
		   uvw[2]*(sqrt(1.-square(l)-square(m))-1.));
	    visbuf*=Complex(cos(phasor),sin(phasor))*
	           casa::Double(real(flux.value()[0])/
		        sqrt(1.-square(l)-square(m)));
	    */
	    // assume that phase terms will be cancelled (with the help of
	    // the dish + A^hA multiplies to the conjugate. We still need
	    // to prove it theoretically).
	    visbuf*=casa::Double(real(flux.value()[0]));

	    for (uInt measurement=0;measurement<vismatrix.nrow();
	         ++measurement) 
		    vismatrix(measurement,feed)+=visbuf;
       }
  }
}

// this version fills vismatrix with sum(F_i E_k*E_l^H), where 
// E is the voltage pattern value at the source position and
// F_i is the flux of the ith source
// pa - parallactic angle to rotate all source offsets (in radians)
// if skycat!="", a table with this name will be filled with the offsets
// w.r.t. the dish pointing centre
void ImplCalWeightSolver::formVPMatrix(const casa::Matrix<casa::Double> &feed_offsets,
                    casa::Double pa, const casa::String &skycat)
                    const throw(casa::AipsError)
{
  if (vp_real==NULL) 
     throw AipsError("A vp image should be set before calling formVPMatrix");

  if (vp_imag==NULL) 
     throw AipsError("A vp image should be set before calling formVPMatrix");

  if (!cl.nelements()) 
      throw AipsError("At least one sky-component should be set");

  AlwaysAssert(feed_offsets.ncolumn()==2,AipsError);
  
  auto_ptr<SkyCatalogTabWriter> psctwr;
  if (skycat!="")
      psctwr.reset(new SkyCatalogTabWriter(skycat));
  

  // resulting matrix is Nfeeds x Nfeeds
  vismatrix.resize(feed_offsets.nrow(),feed_offsets.nrow());
  vismatrix=0.;

  ofstream os("comps.dat");
  for (uInt comp=0;comp<cl.nelements();++comp) {
       const SkyComponent &skc=cl.component(comp);
       if (skc.spectrum().type()!=ComponentType::CONSTANT_SPECTRUM)
	   throw AipsError("Non-constant spectrum sources are not supported");
       // we have to convert reference system of the sky component to
       // that used for the phase/pointing centre
       const MDirection &srcdir_ownref=skc.shape().refDirection();
       // srcdir -> direction in the same reference system as used
       // for the phase/pointing centre
       MVDirection srcdir=MDirection::Convert(srcdir_ownref.getRef(),
		       pc.getRef())(srcdir_ownref.getValue()).getValue();
       // true offsets w.r.t the dish centre
       Double l0_buf=sin(srcdir.getLong()-
                   pc.getValue().getLong())*cos(srcdir.getLat());
       Double m0_buf=sin(srcdir.getLat())*cos(pc.getValue().getLat())-
		 cos(srcdir.getLat())*sin(pc.getValue().getLat())*
	          cos(srcdir.getLong()-pc.getValue().getLong());
       // parallactic angle rotation		  
       Double l0=l0_buf*cos(pa)+m0_buf*sin(pa);
       Double m0=-l0_buf*sin(pa)+m0_buf*cos(pa);
       
       /*
       // to use the code with some dummy model instead of the one
       // loaded from component list (mainly for debugging)
       if (comp==0) {
           l0=M_PI/180.; m0=-M_PI/180.;
       } else if (comp==1) {
           l0=-M_PI/180.; m0=M_PI/180.;       
       } else break;
       */
       
       cout<<"Source "<<comp+1<<" is at "<<setw(10)<<l0/M_PI*180.<<
             " "<<m0/M_PI*180.<<endl;
       
       // component flux
       Flux<Double> flux=skc.flux();
       flux.convertPol(ComponentType::STOKES);
       flux.convertUnit("Jy");
       os<<setw(10)<<l0/M_PI*180*3600<<" "<<m0/M_PI*180*3600<<" "<<real(flux.value()[0])<<endl;

       if (psctwr.get()!=NULL && real(flux.value()[0])>0.3)
           psctwr->addComponent(l0/M_PI*180.,m0/M_PI*180.,
	                      real(flux.value()[0]),"AzEl");
			      
       for (uInt feed1=0;feed1<vismatrix.nrow();++feed1) {
	    Complex visbuf1(1.,0.);
	    // following formulae are approximate. We need to replace 
	    // them with exact formulae for the angular offset of 
	    // the source w.r.t. the appropriate feed centre
	    Double l1=l0-feed_offsets(feed1,0);
	    Double m1=m0-feed_offsets(feed1,1);
	    if (!getVPValue(visbuf1,l1,m1)) continue;
	    for (uInt feed2=0;feed2<vismatrix.ncolumn();++feed2) {
	         Complex visbuf2(1.,0.);
                 Double l2=l0-feed_offsets(feed2,0);
	         Double m2=m0-feed_offsets(feed2,1);
	         if (!getVPValue(visbuf2,l2,m2)) continue;
		 vismatrix(feed1,feed2)+=visbuf1*conj(visbuf2)*
		                 casa::Double(real(flux.value()[0]));
	    }
       }
  }  
}

// solve for a basis in the space of weights, which is best for calibration
// in the sense that the flux from known sources is maximized
// The result matrix contains the basis vectors in its columns
// ndim is a number of basis vectors required (should be less than or
// equal to the number of feeds.
// pa - parallactic angle to rotate all source offsets (in radians)
// if skycat!="", a table with this name will be filled with the offsets
// w.r.t. the dish pointing centre
casa::Matrix<casa::Complex>
ImplCalWeightSolver::calBasis(const casa::Matrix<casa::Double> &feed_offsets,
               casa::uInt ndim, casa::Double pa, const casa::String &skycat)
     	                      const throw(casa::AipsError)
{
  if (ndim>feed_offsets.nrow())
      throw AipsError("The dimension of the basis can not be greater than the number of feeds");
  formVPMatrix(feed_offsets,pa,skycat);
  //QFSumOptimizer qfso;
  //qfso.solve(vismatrix,ndim);
  //return qfso.getBasis();
  return casa::Matrix<casa::Complex>();
}


// solve for eigenvectors for the VP matrix. The first vector
// (column=0) corresponds to the maximum attainable flux,
// the last one (column=Nfeeds-1) corresponds to the
// weight set for an optimal rejection of all known sources.
// pa - parallactic angle to rotate all source offsets (in radians)
// if skycat!="", a table with this name will be filled with the offsets
// w.r.t. the dish pointing centre
casa::Matrix<casa::Complex>
ImplCalWeightSolver::eigenWeights(const casa::Matrix<casa::Double>
                    &feed_offsets, casa::Double pa,
		    const casa::String &skycat)  const throw(casa::AipsError)
{
   formVPMatrix(feed_offsets,pa,skycat);
   /*
    for (uInt feed1=0;feed1<vismatrix.nrow();++feed1)         
        for (uInt feed2=0;feed2<vismatrix.ncolumn();++feed2)
            cout<<"A("<<feed1<<","<<feed2<<"): matrix "<<abs(vismatrix(feed1,feed2))<<endl;
    
    throw AipsError("Debug exception");
   */

   EigenSolver es;
   es.solveEigen(vismatrix);
   
   Vector<Double> eval = es.getEigenValues();
   ofstream os("eval.dat");
   for (uInt k=0;k<eval.nelements();++k) {
        os<<k+1<<" "<<eval[k]<<std::endl;
   }
   /*
   Vector<Complex> vec=es.getEigenVectors().column(0);
   Vector<Complex> tst(vec.nelements(),0.);
   for (uInt k=0;k<tst.nelements();++k)
        for (uInt i=0;i<vec.nelements();++i)
	     tst[k]+=vec[i]*vismatrix(k,i)/es.getEigenValues()[0];
   for (uInt i=0;i<tst.nelements();++i)
        cout<<i<<" "<<tst[i]<<" "<<vec[i]<<endl;
   throw AipsError("Debug exception");
   */

   return es.getEigenVectors();
}


// an auxiliary function to extract a Value of the Voltage Pattern
// at the given offset (in radians). Return True if successful, and
// False if the requested offset lies outside the model
// The parameter val is multiplied by the VP value
Bool ImplCalWeightSolver::getVPValue(Complex &val, Double l, Double m) const
	               throw(AipsError)
{
  Vector<Double> world(vp_real->shape().nelements()); // world coordinate
  // centre of the vp image in pixels
  IPosition cntr=vp_real->shape();
  for (uInt i=0;i<cntr.nelements();++i)
       if (i<2) cntr[i]/=2;
       else cntr[i]=0;

  // just to have something useful for degenerate axes
  vp_real->coordinates().toWorld(world,cntr);
  world[0]+=l; world[1]+=m;
  Vector<Double> pixel; // offset in pixels
  if (!vp_real->coordinates().toPixel(pixel,world))
      throw AipsError(String("access to the real image: ")+
		      vp_real->coordinates().errorMessage());
  cntr[0]=Int(pixel[0]);
  cntr[1]=Int(pixel[1]);
  if (cntr[0]>=vp_real->shape()[0] || cntr[0]<0 ||
      cntr[1]>=vp_real->shape()[1] || cntr[1]<0) return False;
  val*=Complex((*vp_real)(cntr),(*vp_imag)(cntr));    
  return True;
}

casa::Matrix<casa::Complex>
    ImplCalWeightSolver::solveWeights(const casa::Matrix<casa::Double> &feed_offsets,
	     const casa::Vector<casa::Double> &uvw) const
                                  throw(casa::AipsError)
{
  if (vp_real==NULL || vp_imag==NULL)
      throw AipsError("VP image must be set before calling solveWeights");
  formVisMatrix(feed_offsets,uvw);
  /*
  for (uInt meas=0;meas<vismatrix.ncolumn();++meas) {
       cout<<"Measurement "<<meas<<endl;
       for (uInt feed=0;feed<vismatrix.nrow();++feed)
            cout<<"Feed "<<feed<<": visibility "<<abs(vismatrix(meas,feed))<<endl;
  }
  throw AipsError("Debug exception");
  */  
  //EqSolver solver(vismatrix);
  //return solver.solveWeights();
  return casa::Matrix<casa::Complex>();
}

