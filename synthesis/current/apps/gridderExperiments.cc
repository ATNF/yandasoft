
// casa
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>

// System includes
#include <stdexcept>
#include <iostream>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <CommandLineParser.h>

#include <gridding/SphFuncVisGridder.h>
#include <dataaccess/DataAccessorStub.h>
#include <fitting/Axes.h>

ASKAP_LOGGER(logger, ".test");

using namespace askap;
using namespace askap::synthesis;
//using namespace askap::scimath;

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>


inline casa::DComplex getComplex(const gsl_complex gc) {
     return casa::DComplex(GSL_REAL(gc), GSL_IMAG(gc));
}


class testGridder : public SphFuncVisGridder
{
public:
      testGridder() {
          std::cout<<"Test gridder, used for debugging"<<std::endl;
          casa::Matrix<casa::DComplex> B(5,5);
          const double c = casa::C::pi*6./2.;
          /*
          fillMatrixB(B,c,15);
          std::cout<<B<<std::endl;
          
          casa::Vector<casa::DComplex> V(B.nrow());
          casa::DComplex eVal = optimumEigenVector(B,V);
          std::cout<<"eigen value "<<eVal<<" vector: "<<V<<std::endl;
          
          casa::Vector<casa::DComplex> P(V.nelements(),casa::DComplex(0.,0.));
          for (size_t i = 0; i<P.nelements(); ++i) {
               P[i] = -eVal*V[i];
               for (size_t k = 0; k<V.nelements(); ++k) {
                    P[i] += B(i,k)*V(k);
               }
          }
          
          std::cout<<P<<std::endl;
          
          
          casa::Vector<casa::DComplex> vals(6.);
          calcValsAtRegularGrid(vals, V, eVal, false);
                    
          std::cout<<vals<<std::endl;
          */
          //sphFunc(c,1.,0.0001,16); return;
          
          std::ofstream os ("cf.dat");
          const size_t nPoints = 100;
          for (size_t i=0; i<nPoints; ++i) {
               const double x = 2./double(nPoints)*double(i);
               const double sfval = sphFunc(c,1.,x,16);
               //const double sfval = (abs(x)<1 ? real(vals[i])/sqrt(1.-x*x) : 0.);
               os<<x<<" "<<sfval<<" "<<grdsf(x)<<std::endl;
          }
                
      }
      
      /// @brief lth derivative of kth Legendre polynomial at 1.0
      /// @details Calculate the value of the lth derivative of the Legendre
      /// polynomial at 1.0 using recursive formula. It might be possible to
      /// join several loops together and speed the algorithm up a bit, but
      /// we will worry about the optimisation later (if we see that it is 
      /// useful). 
      /// @param[in] l order of the derivative
      /// @param[in] k order of the polynomial
      /// @return value of the derivative at 1.0
      static double derivativeOfLegendrePolynomial(casa::uInt l, casa::uInt k) {
         if (l > k) {
             return 0.;
         }
         // initialise with 0-order derivative
         double res = sqrt(double(2*k+1)/2.);
         for (size_t order = 0; order < l; ++order) {
              res *= (double(k*(k+1)) - double(order*(order+1)))/ double(2*(order+1));
         }
         return res;
      }
      
      /// @brief calculate values at the regular grid
      /// @details The spheroidal function is approximated as a series with coefficients
      /// which are the values at regular grid points pi*N/c. This method fills a vector
      /// with such values spanning N from 0 to size-1. The formulas are slightly different
      /// for odd and even order of spherodial functions.
      /// @param[in] vals values vector to fill, should already be sized to the required size
      /// @param[in] eVec eigen vector in Legendre space
      /// @param[in] eVal eigen value corresponding to eVec
      /// @param[in] isOdd if true, the calculated function is assumed to be of the odd order
      static void calcValsAtRegularGrid(casa::Vector<casa::DComplex> &vals, 
             const casa::Vector<casa::DComplex> &eVec, const casa::DComplex &eVal, const bool isOdd) 
      {
         ASKAPASSERT(vals.nelements()>1);
         ASKAPASSERT(eVec.nelements()>1);
         ASKAPASSERT(casa::abs(eVal) != 0.);
         // the value at 0 is calculated through direct series expansion
         vals.set(0.);
         // all odd Legendre polynomials are anti-symmetric, so will be the spheroidal function
         if (!isOdd) {
             double p0 = 1; // The value of Legendre polynomial at x=0
             for (size_t order=0; 2*order<eVec.nelements(); ++order) {
                  vals[0] += eVec[order] * p0;
                  p0 *= -double(order+1)/double(order+2);
             }
         }
         // now filling the values
         for (size_t N=1; N<vals.nelements(); ++N) {
              for (size_t k = (isOdd ? 1 : 0); k<eVec.nelements(); k+=2) {
                   // Ink is the coefficient in the eigenvector space
                   // see formula (49) in Karoui & Moumni
                   // for the function of an odd order, the value is pure imaginary
                   // so we store just the imaginary part 
                   double Ink = 0;
                   for (size_t l=1; l < k/2; ++l) {
                        if (isOdd) {
                            Ink += negateForOdd(l+1)/pow(casa::C::pi*double(N),2*l+1)*
                                   derivativeOfLegendrePolynomial(2*l,k);
                        } else {
                            Ink += negateForOdd(l+1)/pow(casa::C::pi*double(N),2*l)*
                                   derivativeOfLegendrePolynomial(2*l-1,k);
                        }
                   }
                   Ink *= 2*negateForOdd(N);
                   vals[N] += eVec[k] * Ink * (isOdd ? casa::DComplex(0.,1.) : casa::DComplex(1.,0.));
              }
              // all function values for N>0 should be divided by the eigenvalue
              vals[N] /= eVal;
         }             
      }
      
      /// @brief sum of the Legendre series
      /// @details This helper method sums Legendre series for the given coefficients and the origin
      /// @param[in] coeffs vector with coefficients (element index r is incremented by two)
      /// @param[in] x abcissa 
      /// @param[in] m parameter m of the Legendre function (corresponding to resulting Smn(c,eta))
      /// @param[in] rEven true, if series starts from r=0, false if from r=1 (n-m of Smn is even or odd)
      static double sumLegendreSeries(const casa::Vector<double> &coeffs, double x, int m, bool rEven) {
           ASKAPASSERT(m>=0);
           const int nOrders = int(coeffs.nelements())*2 + (rEven ? 0 : 1);
           double *vals = new double[nOrders+1];
           
           int status = gsl_sf_legendre_sphPlm_array(nOrders + m, m, x, vals);
           double result = 0.;
           for (casa::uInt elem = 0; elem<coeffs.nelements(); ++elem) {
                const int r = 2*elem + (rEven ? 0 : 1);
                //const int l = r + m;
                ASKAPASSERT(r < nOrders + 1);
                result += coeffs[elem]*vals[r];
           }
           
           delete[](vals);
           ASKAPCHECK(status == GSL_SUCCESS, "Error calculating associated Legendre functions, status="<<status);           
           return -result;
      }
      
      /// @brief calculate spheroidal function via Legendre decomposition
      /// @details This algorithm decomposes the spheroidal function into series with associated Legendre functions
      /// with some special normalisation
      /// @param[in] c parameter c of the spheroidal function (bandwidth or a measure of the support size in our case)
      /// @param[in] alpha parameter alpha of the spheroidal function (weighting exponent in our case)
      /// @param[in] eta argument of the function
      /// @param[in] nterms number of terms in the decomposition
      static double sphFunc(const double c, const double alpha, const double eta, const casa::uInt nterms)
      {
        if (fabs(double(eta))>=1) {
            return 0.;
        }
        casa::Matrix<double> hlp(nterms,nterms,0.);
        const bool rEven = true;      
        fillHelperMatrix(hlp,c,int(alpha),rEven);
        
        casa::Vector<double> coeffs;
        legendreCoeffs(hlp,coeffs);
        
        // force normalisation to 1. at eta=0., functions corresponding to n=0 are even, so such normalisation
        // should not cause any problems
        const double res = sumLegendreSeries(coeffs,eta,int(alpha),rEven) / sumLegendreSeries(coeffs,0.,int(alpha),rEven);
        return res * pow(1.-eta*eta, -alpha/2.);
      }
      
      /// @brief calculate spheroidal function via Bessel decomposition
      /// @details This algorithm decomposes the spheroidal function into series with Bessel functions
      /// @param[in] c parameter c of the spheroidal function (bandwidth or a measure of the support size in our case)
      /// @param[in] alpha parameter alpha of the spheroidal function (weighting exponent in our case)
      /// @param[in] eta argument of the function
      /// @param[in] nterms number of terms in the decomposition
      /// @param[in] mSize optional matrix size for the dependent eigenproblem, 0 means the miminal size 
      ///                  sufficient to produce nterms in the decomposition
      static double sphFunc1(const double c, const double alpha, const double eta, const casa::uInt nterms, 
                     const casa::uInt mSize = 0)
      {
        ASKAPCHECK(alpha>-0.5, "The case of alpha<=-0.5 has not been tested (although might work), you have alpha="<<alpha);
        ASKAPASSERT(nterms>1);
        if (eta == 0.) {
            return 1.;
        }
        casa::Vector<double> coeffs(nterms, 0.);
        calcBesselCoeffs(c, alpha, coeffs, mSize);
        ASKAPDEBUGASSERT(coeffs.nelements() == nterms);
        
        // value at (0,0) used for normalisation    
        //const double sfAt0_0 = sphFuncAt0_0(alpha,coeffs);
        
        // first order of Bessel function in the series
        const double startOrder = alpha >= -0.5 ? alpha + 0.5 - int(alpha+0.5) : alpha + 0.5;
        const int nBesselVals = alpha >= -0.5 ? int(alpha+1.5+2*nterms) : 2*int(nterms);
        ASKAPASSERT(nBesselVals > 0);
        ASKAPASSERT(int(nBesselVals) >= 2*int(nterms)+int(alpha+0.5));

        casa::Vector<double> besselVals(nBesselVals,0.);
        // calculate series of Bessel function values, orders go from startOrder to startOrder+nBesselVals-1
        bessel(startOrder,c*std::abs(eta),besselVals);
        
        double sum = 0.;
        for (casa::uInt i = 0; i<nterms; ++i) {
             const int index = 2*int(i) + int(alpha+0.5);
             ASKAPCHECK((index>=0) && (index<nBesselVals), "Invalid index = "<<index<<", alpha="<<alpha<<" term="<<i);
             sum += coeffs[i] * besselVals[index];
        }
        ASKAPASSERT(coeffs[0]!=0.);
        sum *= pow(2./(c*std::abs(eta)),alpha+0.5)*gamma(alpha+1.5)/coeffs[0];
        
        return sum;
      }
      
      /// @brief calculate spheroidal function at (0,0)
      /// @details This helper method calculates the value of spheroidal function for c=0 and eta=0
      /// @param[in] alpha parameter alpha of the spheroidal function (weighting exponent in our case)
      /// @param[in] coeffs vector with the coefficients
      /// @return value of the function
      static double sphFuncAt0_0(const double alpha, const casa::Vector<double> &coeffs) {
          double res = 0.;
          for (casa::uInt i=0; i<coeffs.nelements(); ++i) {
               res += gamma(0.5+double(i))/gamma(alpha+double(i+1))*coeffs[i];
          }
          res *= casa::C::_1_sqrtpi / pow(2.,alpha);
          return res;
      }
      
      /// @brief gamma function 
      /// @details this is a helper wrapper over GSL's gamma function implementation
      /// @param[in] x argument of the gamma function
      /// @return value of the gamma function
      /// @note an exception is thrown if there is an error
      static double gamma(double x) {
         gsl_sf_result res;
         const int status = gsl_sf_gamma_e(x, &res);
         ASKAPCHECK(status == GSL_SUCCESS, "Error in calculation of gamma function for x="<<x<<", status="<<status);
         return res.val;
      }
      
      /// @brief regular cylindrical Bessel function
      /// @details This is a wrapper on top of the GSL routine to calculate Bessel function. Ideally we
      /// want to be able to generate a sequence of functions of different orders, but at the same 
      /// argument (series expansion). Libraries seem to provide calculation of a sequence at different 
      /// arguments, but for the same order. So this code may be sub-optimal. But it is fine for now.
      /// Another optimisation which might be possible is to take into account that the only order of 
      /// Bessel function we use is half plus integer.
      /// @param[in] nu order of the Bessel function (can be fractional)
      /// @param[in] x argument
      /// @return value of the Bessel function
      /// @note an exception is thrown if there is an error
      static double bessel(const double nu, const double x) {
         gsl_sf_result res;
         const int status = gsl_sf_bessel_Jnu_e(nu,x, &res);
         ASKAPCHECK(status == GSL_SUCCESS, "Error in calculation of Bessel function for nu="<<nu<<" x="<<x<<", status="<<status);
         return res.val;       
      }
      
      /// @brief set of values for Bessel function
      /// @details This version calculates a set of values for the regular cylindrical Bessel function
      /// for a sequence of orders.
      /// @param[in] nu starting order
      /// @param[in] x argument
      /// @param[in] vals vector to be filled with values corresponding to orders nu,nu+1,...nu+vals.nelements()-1
      /// @note an exception is thrown if there is an error. Vector vals must have non-zero size to indicate the number
      /// of orders required
      static void bessel(const double nu, const double x, casa::Vector<double> &vals) {
         ASKAPASSERT(vals.nelements()>=1);
         for (casa::uInt k = 0; k<vals.nelements(); ++k) {
              vals[k] = bessel(nu+double(k),x);
         } 
      }
      
      /// @brief Bessel series expansion coefficients
      /// @details This is a helper method to compute series coefficients for decomposition
      /// of a given spheroidal functions via Bessel functions
      /// @param[in] c parameter c of the spheroidal function (bandwidth or a measure of the support size in our case)
      /// @param[in] alpha parameter alpha of the spheroidal function (weighting exponent in our case)
      /// @param[in] coeffs vector to fill with the coefficients (must already be resized to a 
      ///                   required number of coefficients)
      /// @param[in] mSize optional matrix size for the dependent eigenproblem, 0 means the miminal size sufficient to
      ///            produce coeffs.nelements() coefficient. Positive number should not be below coeffs.nelements().
      static void calcBesselCoeffs(const double c, const double alpha, casa::Vector<double> &coeffs,
                                   const casa::uInt mSize = 0)
      {
        ASKAPASSERT(coeffs.nelements()>1);
        const casa::uInt matrSize = (mSize == 0 ? coeffs.nelements() : mSize);
        ASKAPCHECK(matrSize >= coeffs.nelements(), "Requested matrix size of "<<matrSize<<
                   " should not be less than the number of requested coefficients ("<<coeffs.nelements()<<")");
        ASKAPCHECK(2*alpha != -3.,"Implemented formulas don't work for alpha = -1.5");
        const double cSquared = c*c;
        // buffers
        casa::Vector<double> bufA(2*matrSize+1,0.);                   
        casa::Vector<double> bufB(2*matrSize+1,0.);
        casa::Vector<double> bufC(2*matrSize+1,0.);
        casa::Vector<double> diag(matrSize,0.);
        casa::Vector<double> sdiag2(matrSize-1,0.);
        
        // fill the buffers
        bufB[0] = cSquared / (2*alpha+3);
        bufC[0] = cSquared * (2*alpha+2) / (2*alpha+3);
        for (casa::uInt k = 2; k <= 2*matrSize; k+=2) {
             // k is a 1-based index into buffers
             bufA[k] = cSquared * double(k*(k-1))/(2*alpha+2*k-1)/(2*alpha+2*k+1);
             bufB[k] = cSquared * (double(k)*(2*alpha+k+1)+(2*alpha-1+2*k*k+2*k*(2*alpha+1)))/
                                  (2*alpha+2*k-1) / (2*alpha+2*k+3);
             bufC[k] = cSquared * (2*alpha+k+1) * (2*alpha+k+2) / (2*alpha+2*k+1) / (2*alpha+2*k+3);
             // recursion relation for matrix coefficients
             ASKAPDEBUGASSERT(k>=2);
             diag[k/2-1] = bufB[k-2];
             // sub-diagonal has one less element, exclude the last one
             if (k<2*matrSize) {
                 sdiag2[k/2] = bufA[k] * bufC[k-2];
             }
        }
        const double eVal = smallestEigenValue(diag,sdiag2);
        coeffs[0] = 1.; // arbitrary scaling factor which will be normalised later when we sum the series
        // now calculate the ratios of adjacent coefficients
        for (casa::uInt elem = coeffs.nelements() - 1; elem>0; --elem) {
             double factor = eVal - bufB[2*elem];
             if (elem + 1 < coeffs.nelements()) {
                 factor += bufA[2*elem+2]*coeffs[elem+1];
             }
             if (factor != 0) {
                 coeffs[elem] = -bufC[2*elem-2]/factor;
             } else {
                 coeffs[elem] = 0.;
             }
        }  
        std::cout<<coeffs<<std::endl;
        // and bootstrap all coefficients using the ratios and the first arbitrary defined element
        for (casa::uInt elem = 1; elem<coeffs.nelements(); ++elem) {
             coeffs[elem] *= coeffs[elem-1];
        }
      }
      
      /// @brief fill matrix which has the same eigenvalues/vectors as the original problem
      /// @details See equation (20) in Aquino and Casta\~no (2002)
      /// @param[in] B matrix to fill (should already be sized to required number of terms)
      /// @param[in] c bandwidth of the prolate spheroidal function
      /// @param[in] m parameter m of the prolate spheroidal function Smn(c,eta)
      /// @param[in] rEven true, if series starts from r=0, false if from r=1 (n-m of Smn is even or odd)
      static void fillHelperMatrix(casa::Matrix<double> &B, const double c,int m, bool rEven) {
          ASKAPASSERT(B.nrow() == B.ncolumn());
          ASKAPASSERT(B.nrow() > 1);
          const double cSquared = c*c;
          for (casa::uInt row = 0; row<B.nrow(); ++row) {
               const int r = int(row)*2 + (rEven ? 0 : 1);
               const int l = r + m; // order of Legendre function P_l^m
               B(row,row) = double(l*(l+1)) + cSquared*(double(2*l+3)*(l+m)*(l-m)+double(2*l-1)*(l+m+1)*(l-m+1)) /
                               (double(2*l+1)*(2*l-1)*(2*l+3));
               if (row>=1) {
                   B(row,row-1) = cSquared/double(2*l-1)*sqrt(double(l+m)*(l+m-1)*(l-m)*(l-m-1)/(double(2*l+1)*(2*l-3)));
               }
               if (row+1<B.nrow()) {
                   B(row,row+1) = cSquared/double(2*l+3)*sqrt(double(l+m+1)*(l+m+2)*(l-m+1)*(l-m+2)/
                                  (double(2*l+1)*(2*l+5)));
               }
               
          }
      }      
      
      /// @brief coefficients in Legendre series
      /// @details This method solves eigenvalue problem and obtains eigenvector corresponding to
      /// the smallest eigenvalue (for function Smn(c,eta) this means n=0). Coefficients are in the same order
      /// as elements of matrix B, i.e. in steps of 2 starting from even or odd depending whether n-m is even or odd
      /// @param[in] B matrix to solve
      /// @param[in] coeffs output coefficients for Legendre series (to be resized to match the size of B)
      /// @return eigenvalue
      /// @note an exception is thrown if there is an error solving eigensystem
      static double legendreCoeffs(const casa::Matrix<double> &B, casa::Vector<double> &coeffs) 
      {
         ASKAPASSERT(B.nrow() == B.ncolumn());
         coeffs.resize(B.nrow());
         
         gsl_matrix *A = gsl_matrix_alloc(B.nrow(),B.nrow());
         gsl_matrix *eVec = gsl_matrix_alloc(B.nrow(),B.nrow());
         //gsl_matrix_set_zero(A);         
         gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(B.nrow());
         gsl_vector *eVal = gsl_vector_alloc(coeffs.nelements());

         // fill the matrix (a bit of an overkill, but it is faster to reuse existing code
         // than to write something for tridiagonal matrix)
         for (casa::uInt row = 0; row<B.nrow(); ++row) {
              for (casa::uInt col = 0; col<B.ncolumn(); ++col) {
                   gsl_matrix_set(A, row, col, B(row,col));
              }
         }
         
         const int status = gsl_eigen_symmv(A,eVal,eVec,work);
         double result = -1.;
         casa::uInt optIndex = 0;
         if (status == GSL_SUCCESS) {
             for (casa::uInt elem = 0; elem<coeffs.nelements(); ++elem) {
                  const double val = gsl_vector_get(eVal,elem);
                  if ((elem == 0) || (val<result)) {
                      result = val;
                      optIndex = elem;
                  }
             }
         }
         
         // extract the appropriate eigenvector
         for (size_t i=0; i<B.nrow(); ++i) {
              coeffs[i] = gsl_matrix_get(eVec,i,optIndex);                             
         }         

         gsl_matrix_free(A);         
         gsl_matrix_free(eVec);         
         gsl_eigen_symmv_free(work);
         gsl_vector_free(eVal);
         
         ASKAPCHECK(status == GSL_SUCCESS, "Error solving eigenproblem in legendreCoeffs, status="<<status);
         
         /*
         // consistency check
         casa::Vector<double> test(coeffs.nelements(),0);
         for (size_t i=0;i<coeffs.nelements(); ++i) {
              for (size_t k=0;k<coeffs.nelements(); ++k) {
                   test[i]+=B(i,k)*coeffs[k];
              }
              test[i]-=result*coeffs[i];
         }
         std::cout<<test<<std::endl;
         */
         
         return result;         
      }
    
      
      /// @brief smallest eigenvalue of a symmetric tridiagonal matrix
      /// @details This helper method finds the smallest eigenvalue of a symmetric tridiagonal
      /// matrix.
      /// @param[in] diag main diagonal of the matrix
      /// @param[in] sdiag2 squares of the subdiagonal of the matrix
      /// @return smallest eigenvalue
      static double smallestEigenValue(const casa::Vector<double> &diag, const casa::Vector<double> &sdiag2) 
      {
         ASKAPASSERT(diag.nelements() == sdiag2.nelements() + 1);
         ASKAPASSERT(diag.nelements() > 1);
         
         gsl_matrix *A = gsl_matrix_alloc(diag.nelements(),diag.nelements());
         gsl_matrix_set_zero(A);         
         gsl_eigen_symm_workspace *work = gsl_eigen_symm_alloc(diag.nelements());
         gsl_vector *eVal = gsl_vector_alloc(diag.nelements());
         
         // fill the matrix (a bit of an overkill, but it is faster to reuse existing code
         // than to write something for tridiagonal matrix)
         for (casa::uInt elem = 0; elem<diag.nelements(); ++elem) {
              gsl_matrix_set(A, elem, elem, diag[elem]);
              if ((elem + 1 < diag.nelements()) != 0) {
                  ASKAPASSERT(sdiag2[elem]>=0.);
                  gsl_matrix_set(A, elem, elem+1, sqrt(fabs(double(sdiag2[elem]))));              
                  gsl_matrix_set(A, elem+1, elem, sqrt(fabs(double(sdiag2[elem]))));              
              }
         }
         const int status = gsl_eigen_symm(A,eVal, work);
         double result = -1.;
         if (status == GSL_SUCCESS) {
             for (casa::uInt elem = 0; elem<diag.nelements(); ++elem) {
                  const double val = gsl_vector_get(eVal,elem);
                  if ((elem == 0) || (val<result)) {
                      result = val;
                  }
             }
         }

         gsl_matrix_free(A);         
         gsl_eigen_symm_free(work);
         gsl_vector_free(eVal);
         
         ASKAPCHECK(status == GSL_SUCCESS, "Error solving eigenproblem for symmetric tridiagonal matrix, status="<<status);
         return result;
      }
            
      /// @brief do eigen decomposition, get optimum eigen vector/value
      /// @details Solve for eigenvalues and eigen vectors of the helper matrix,
      /// Find the largest by absolute value and extract appropriate eigenvector
      /// @param[in] B matrix to decompose      
      /// @param[out] V optimum eigenvector (will be resized to B.nrow())
      /// @return largest eigenvalue (by absolute value)
      static casa::DComplex optimumEigenVector(const casa::Matrix<casa::DComplex> &B, casa::Vector<casa::DComplex> &V)
      {
         ASKAPASSERT(B.nrow() == B.ncolumn());
         ASKAPASSERT(B.nrow() > 1);
         gsl_matrix *A = gsl_matrix_alloc(B.nrow()*2,B.ncolumn()*2);
         gsl_vector_complex *eVal = gsl_vector_complex_alloc(B.nrow()*2);
         gsl_eigen_nonsymmv_workspace *work = gsl_eigen_nonsymmv_alloc(B.nrow()*2);
         gsl_matrix_complex *eVec = gsl_matrix_complex_alloc(B.nrow()*2,B.ncolumn()*2);
         
         for (size_t row=0; row<B.nrow(); ++row) {
              for (size_t col=0; col<B.ncolumn(); ++col) {
                   const double reB = casa::real(B(row,col));
                   const double imB = casa::imag(B(row,col));
                   gsl_matrix_set(A, row*2, col*2, reB);
                   gsl_matrix_set(A, row*2+1, col*2+1, reB);
                   gsl_matrix_set(A, row*2, col*2+1, -imB);
                   gsl_matrix_set(A, row*2+1, col*2, imB);                   
              }
         }
         const int status = gsl_eigen_nonsymmv(A,eVal,eVec, work);
         
         casa::Complex peakVal(0.,0.);
         if (status == 0) {
             // eigenproblem solved successfully
             size_t peakIndex = 0;
             
             // search for peak eigenvalue
             for (size_t el=0; el<B.nrow()*2; ++el) {
                  casa::DComplex val = getComplex(gsl_vector_complex_get(eVal,el));              
                  std::cout<<"el="<<el<<" "<<val<<std::endl;
                  if ((el == 0) || (casa::abs(val) > casa::abs(peakVal))) {
                      peakIndex = el;
                      peakVal = val;
                  }
             }
             
             /*
             peakIndex = 7;
             peakVal = getComplex(gsl_vector_complex_get(eVal,peakIndex));              
             */
             
             std::cout<<"peak Index="<<peakIndex<<" peakValue="<<peakVal<<std::endl;
             // extract the appropriate eigenvector
             V.resize(B.nrow());
             for (size_t i=0; i<B.nrow(); ++i) {
                  V[i] = getComplex(gsl_matrix_complex_get(eVec,  2*i,peakIndex))+casa::DComplex(0.,1.)*
                             getComplex(gsl_matrix_complex_get(eVec, 2*i+1,peakIndex));                             
             }
         }
         gsl_matrix_complex_free(eVec);
         gsl_eigen_nonsymmv_free(work);
         gsl_vector_complex_free(eVal);
         gsl_matrix_free(A);
         ASKAPCHECK(status == 0, "Eigen problem solution has failed in optimumEigenVector");
         return peakVal;
      }
      
      /// @brief helper method to evaluate (-1)^l
      /// @param[in] l integer power index
      /// @return 1 if l is even, -1 otherwise
      static inline int negateForOdd(const casa::uInt l) { return (l%2 == 0) ? 1 : -1; } 
      
      /// @brief fill matrix B which has the same eigenvalues as the original problem
      /// @details See equation (8) in Karoui & Moumni (2008)
      /// @param[in] B matrix to fill (should already be sized)
      /// @param[in] c bandwidth of the prolate spheroidal function
      /// @param[in] nterms number of terms in the series expansion
      static void fillMatrixB(casa::Matrix<casa::DComplex> &B, const double c, const int nterms = 15) {
          ASKAPASSERT(B.nrow() == B.ncolumn());
          ASKAPASSERT(B.nrow() > 1);
          ASKAPASSERT(nterms>1);
          // supplementary matrix is used to hold moments (row is the moment number, starting from 0),
          // of the normalised Legendre polynomials (column is the order of the polynomial, starting from 0)
          // Note, the matrix is rectangular. The last row is not used to fill
          // the matrix B, but is required to construct other elements of this supplementary matrix through
          // the recursion formula.
          casa::Matrix<double> moments(nterms+1,B.nrow(),0.);
          ASKAPASSERT(moments.ncolumn() >= 2);
          // fill the first two columns
          for (casa::uInt l=0; l<moments.nrow(); ++l) {
               moments(l,0) = (1. + negateForOdd(l))/(sqrt(2.)*(l + 1));
               moments(l,1) = sqrt(1.5)*(1. + negateForOdd(l+1))/(l + 2);               
          }
          // now fill other columns, if any
          for (casa::uInt k=1; k + 1 < moments.ncolumn(); ++k) {
               for (casa::uInt l=0; l + 1 < moments.nrow(); ++l) {
                    moments(l,k+1) = sqrt(double((2*k+1)*(2*k+3))/(k+1)/(k+1)) * moments(l+1,k) -
                          double(k)/(k+1)*sqrt(double(2*k+3)/(2*int(k)-1)) * moments(l, k-1);
               } 
          }
          
          // now fill the matrix B (approximation of the matrix for the Helmoltz equation operator in the 
          // Legendre basis)
          double coeff = 1.; // (c^l / l!)
          B.set(casa::DComplex(0.,0.));
          for (int l = 0; l<nterms; ++l) {
               if (l != 0) {
                   coeff *= c / l;
               }
               // i^l goes in sequence 1,i,-1,-i,1,...
               const casa::DComplex iPwrl = casa::DComplex((l % 4 > 1) ? -1. : 1.) * 
                             ((l % 2 == 1) ? casa::DComplex(0.,1.) : casa::DComplex(1.,0.));
               // fill the actual elements of the matrix              
               for (casa::uInt row = 0; row<B.nrow(); ++row) {
                    for (casa::uInt col = 0; col<B.ncolumn(); ++col) {
                         B(row,col) += iPwrl * coeff * moments(l,row) * moments(l, col); 
                    }
               }
          }
      }
};

// Main function
int main(int argc, const char** argv)
{

  try {
  
        // Put everything in scope to ensure that all destructors are called
        // before the final message
        {
            cmdlineparser::Parser parser; // a command line parser
            testGridder gridder;
            
            const double cellSize=10*casa::C::arcsec;

            casa::Matrix<double> xform(2,2,0.);
            xform.diagonal().set(1.);
               
            scimath::Axes axes;
            axes.addDirectionAxis(casa::DirectionCoordinate(casa::MDirection::J2000, 
                      casa::Projection(casa::Projection::SIN), 0.,0.,cellSize,cellSize,xform,256.,256.));
        
            accessors::DataAccessorStub acc(true);
            
            const casa::IPosition shape(4,256,256,1,1);
            gridder.initialiseGrid(axes,shape,false);
            gridder.grid(acc);
            casa::Array<double> grid;
            gridder.finaliseGrid(grid);        
        }
    } catch (const cmdlineparser::XParser &ex) {
        ASKAPLOG_FATAL_STR(logger, "Command line parser error, wrong arguments " << argv[0]);
        std::cerr << "Usage: " << argv[0] << " [-inputs parsetFile]"
                      << std::endl;
    } catch (const askap::AskapError& x) {
        ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
        std::cerr << "Askap error in " << argv[0] << ": " << x.what()
                      << std::endl;
        exit(1);
    } catch (const std::exception& x) {
        ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
        std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what()
                      << std::endl;
        exit(1);
    }

    return 0;
}

