/// @file
/// 
/// @brief Unit tests for SynthesisParamsHelper.
/// @details SynthesisParamsHelper class contains utilities to simplify
/// handling of parameters representing images. This unit test is intended to
/// test this functionality.
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

#ifndef SYNTHESIS_PARAMS_HELPER_TEST_H
#define SYNTHESIS_PARAMS_HELPER_TEST_H

#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <utils/ImageUtils.h>
#include <measurementequation/ImageParamsHelper.h>
#include <cppunit/extensions/HelperMacros.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>


#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

#include <vector>
#include <map>
#include <string>

namespace askap
{
  namespace synthesis
  {
    
    class SynthesisParamsHelperTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(SynthesisParamsHelperTest);
      CPPUNIT_TEST(testListFacet);
      CPPUNIT_TEST(testListTaylor);
      CPPUNIT_TEST(testFacetCreationAndMerging);
      CPPUNIT_TEST(testCoordinates);
      CPPUNIT_TEST(testClipImage);
      CPPUNIT_TEST(testGaussianPreconditioner);
      //CPPUNIT_TEST(test4Debugging);
      CPPUNIT_TEST_SUITE_END();
      
      private:

      public:
        void test4Debugging()
        {
           ifstream is("temp.dat");
           ASKAPASSERT(is);
           int npoints;
           is>>npoints;
           casa::Matrix<casa::Float> arr(npoints,npoints,0.);
           for (int x=0;x<npoints;++x) {
                for (int y=0; y<npoints;++y) {
                     ASKAPASSERT(is);
                     float buf;
                     is>>buf;
                     arr(x,y)=buf;
                }
           }
           scimath::saveAsCasaImage("temp.img",arr);
        }

        void testListFacet()
        {
           std::vector<std::string> names;
           std::map<std::string,int> facetmap;
           names.push_back("image.i.src.facet.0.0");
           names.push_back("image.i.src.facet.0.1");
           names.push_back("image.i.src.facet.1.0");
           names.push_back("image.i.src.facet.1.1");
           names.push_back("image.i.src2");
           SynthesisParamsHelper::listFacets(names,facetmap);
           CPPUNIT_ASSERT(facetmap.size()==2);
           CPPUNIT_ASSERT(facetmap.find("image.i.src")!=facetmap.end());
           CPPUNIT_ASSERT(facetmap["image.i.src"] == 2);
           CPPUNIT_ASSERT(facetmap.find("image.i.src2")!=facetmap.end());
           CPPUNIT_ASSERT(facetmap["image.i.src2"] == 1);           
        }
        
        void testListTaylor()
        {
           std::vector<std::string> names;
           std::map<std::string,int> taylormap;
           names.push_back("image.src.taylor.0");
           names.push_back("image.src.taylor.1");
           names.push_back("image.src.taylor.2");
           names.push_back("image.src2");
           SynthesisParamsHelper::listTaylor(names,taylormap);
           CPPUNIT_ASSERT(taylormap.size()==2);
           CPPUNIT_ASSERT(taylormap.find("image.src")!=taylormap.end());
           CPPUNIT_ASSERT(taylormap["image.src"] == 3);
           CPPUNIT_ASSERT(taylormap.find("image.src2")!=taylormap.end());
           CPPUNIT_ASSERT(taylormap["image.src2"] == 1);                     
        }
        
        void testGaussianPreconditioner()
        { 
           askap::scimath::Params params;
           makeParameter(params,"psf.testsrc",1);
           const askap::scimath::Axes axes = params.axes("psf.testsrc");
           CPPUNIT_ASSERT(axes.hasDirection());
           casa::Vector<casa::Double> increments = axes.directionAxis().increment();
           CPPUNIT_ASSERT(increments.nelements() == 2);
           CPPUNIT_ASSERT(fabs(fabs(increments[0])-fabs(increments[1]))<1e-6);
           casa::IPosition shape = params.value("psf.testsrc").shape();
           CPPUNIT_ASSERT(shape.nonDegenerate().nelements() == 2);
           CPPUNIT_ASSERT(shape[0] == shape[1]);
           casa::Array<float> dirty(shape,0.);
           casa::Array<float> psfArray(shape,0.);
           casa::Matrix<float> psf(psfArray.nonDegenerate());
           psf(shape[0]/2,shape[1]/2) = 1.;           
           const double factor = 4.*log(2.) * fabs(increments[0]) * shape[0] / casa::C::pi;
           GaussianTaperPreconditioner gp(factor/SynthesisParamsHelper::convertQuantity("20arcsec","rad"));
           gp.doPreconditioning(psfArray,dirty);  
           
           // update the parameter
           casa::Array<double> temp(psfArray.shape());
           casa::convertArray<double,float>(temp,psfArray);
           params.update("psf.testsrc",temp);
           //

           casa::Vector<casa::Quantum<double> > fit = SynthesisParamsHelper::fitBeam(params,0.05,"psf.testsrc");
           CPPUNIT_ASSERT(fit.nelements() == 3);
           // the cell size is 8 arcsec, so the tolerance of 0.5 arcsec seems good enough
           CPPUNIT_ASSERT(fabs(fit[0].getValue("arcsec")-20.)<0.5);
           CPPUNIT_ASSERT(fabs(fit[1].getValue("arcsec")-20.)<0.5);
           /*
           // writing the image for debugging
           LOFAR::ParameterSet parset;
           SynthesisParamsHelper::setUpImageHandler(parset);
           SynthesisParamsHelper::saveImageParameter(params,"psf.testsrc","test.img");
           */
        }
        
        void testFacetCreationAndMerging()
        {
           askap::scimath::Params params;
           makeParameter(params,"testsrc",2,128);
           // checking the content
           std::map<std::string,int> facetmap;
           SynthesisParamsHelper::listFacets(params.freeNames(),facetmap);
           CPPUNIT_ASSERT(facetmap.find("testsrc")!=facetmap.end());
           CPPUNIT_ASSERT(facetmap["testsrc"] == 2);
           // adding a merged image
           SynthesisParamsHelper::add(params,"testsrc",2);
           params.fix("testsrc");
           CPPUNIT_ASSERT(params.freeNames().size() == 4); 
           CPPUNIT_ASSERT(params.names().size() == 5); 
           const std::vector<std::string> &facets = params.freeNames();
           for (std::vector<std::string>::const_iterator ci = facets.begin(); ci!=facets.end(); ++ci) {
                 CPPUNIT_ASSERT(SynthesisParamsHelper::getFacet(params,*ci).shape() == 
                          casa::IPosition(4,128,128,1,1));
           }
        }
        
        void testClipImage()
        {
           askap::scimath::Params params;
           const int facetStep = 128;
           makeParameter(params,"testsrc",2,facetStep);
           params.value("testsrc.facet.0.0").set(1.);
           SynthesisParamsHelper::clipImage(params,"testsrc.facet.0.0");
           casa::Array<double> arr = params.value("testsrc.facet.0.0");
           const casa::IPosition shape = arr.shape(); 
           ASKAPDEBUGASSERT(shape.nelements()>=2);
           casa::IPosition index(shape.nelements(),0);                     
           for (index[0]=0; index[0]<shape[0]; ++index[0]) {
                for (index[1]=0; index[1]<shape[1]; ++index[1]) {
                     const bool isCentre = (index[0]>=(shape[0]-facetStep)/2) && 
                                     (index[0]<(shape[0]+facetStep)/2) &&
                                     (index[1]>=(shape[1]-facetStep)/2) && 
                                     (index[1]<(shape[1]+facetStep)/2);
                     CPPUNIT_ASSERT((arr(index)>0.5) == isCentre);                
                }
           }
        }
        
        void testCoordinates()
        {
           doCoordinateAlignmentTest(128,2);
           doCoordinateAlignmentTest(256,2);
           doCoordinateAlignmentTest(64,3);           
        }
        
      protected:
        /// @brief actual test of coordinate alignment
        /// @details
        /// @param[in] facetStep offset between two adjacent facets in pixels
        /// @param[in] nFacets number of facets along each axis
        void doCoordinateAlignmentTest(const int facetStep, const int nFacets)
        {
           askap::scimath::Params params;
           makeParameter(params,"testsrc",nFacets,facetStep);
           // adding a merged image
           SynthesisParamsHelper::add(params,"testsrc",nFacets);

           for (int facetX = 0; facetX<nFacets; ++facetX) {
                for (int facetY = 0; facetY<nFacets; ++facetY) {
                     ImageParamsHelper iph("testsrc",facetX,facetY);
                     
                     casa::IPosition blc(4,0),trc(4,0);
                     
                     //std::cout<<std::endl<<"facet "<<facetX<<" "<<facetY<<std::endl;
                     
                     blc[0] = iph.facetX()*facetStep;
                     trc[0] = blc[0]+facetStep-1;
                     blc[1] = iph.facetY()*facetStep;
                     trc[1] = blc[1]+facetStep-1;
                     //std::cout<<"blc="<<blc<<" trc="<<trc<<std::endl;                     

                     casa::IPosition blc2,trc2;
                     getCorners(params,iph.name(),iph.paramName(),facetStep, blc2,trc2);
                     CPPUNIT_ASSERT((blc2.nelements()>=2) && (trc2.nelements()>=2));
                     //std::cout<<"blc="<<blc2<<" trc="<<trc2<<std::endl;               

                     // there is no exact match between two images, although we're using
                     // the same projection. Probably it is the second order effect 
                     // resulted from approximation of the sphere by a plane.
                     CPPUNIT_ASSERT(casa::abs(blc[0]-blc2[0])<=5);
                     CPPUNIT_ASSERT(casa::abs(blc[1]-blc2[1])<=5);
                     CPPUNIT_ASSERT(casa::abs(trc[0]-trc2[0])<=5);
                     CPPUNIT_ASSERT(casa::abs(trc[1]-trc2[1])<=5);
                     
                }
           }
        }
        
        /// @brief a helper method to find corners of the patch in a bigger image
        /// @details
        /// If this method is proved to be useful, it can be moved to SynthesisParamsHelper.
        /// @param[in] params parameter container
        /// @param[in] fullName name of the full image
        /// @param[in] patchName name of the patch
        /// @param[in] patchSize the size of the patch to work with, the array stored under
        ///                       patchName in the parameter container can have a larger size.
        ///                      At this stage the same size is assumed for both directional axes
        ///                      and the array is centered.
        /// @param[out] blc bottom left corner of the patch inside the full image
        /// @param[out] trc top right corner of the patch inside the full image
        static void getCorners(askap::scimath::Params &params, const std::string &fullName,
                     const std::string &patchName, const int patchSize, 
                     casa::IPosition &blc, casa::IPosition &trc)
        {
           const casa::Array<double> fullImage = params.value(fullName);
      
           blc = fullImage.shape();
           trc = fullImage.shape();
           CPPUNIT_ASSERT(blc.nelements()>=2);
           // adjust extra dimensions
           for (size_t i=2;i<blc.nelements();++i) {
                blc[i] = 0;
                CPPUNIT_ASSERT(trc[i]!=0);
                trc[i] -= 1;           
           }
      
           const casa::IPosition patchShape = params.value(patchName).shape();
           ASKAPDEBUGASSERT(patchShape.nelements()>=2);
           ASKAPDEBUGASSERT((patchSize<=patchShape[0]) && (patchSize<=patchShape[1]));
      
           ASKAPDEBUGASSERT(patchSize>=1);

           
           const casa::DirectionCoordinate csPatch = 
                    SynthesisParamsHelper::directionCoordinate(params,patchName);
           const casa::DirectionCoordinate csFull = 
                    SynthesisParamsHelper::directionCoordinate(params,fullName);
           casa::Vector<double> world(2);
      
           // first get blc
           casa::Vector<double> blcPixel(2);
           blcPixel(0)=double((patchShape[0]-patchSize)/2);
           blcPixel(1)=double((patchShape[1]-patchSize)/2);
          
           //std::cout<<blcPixel<<endl;
          
           csPatch.toWorld(world,blcPixel);
           csFull.toPixel(blcPixel,world);
          
           //std::cout<<blcPixel<<endl;

           // now get trc
           casa::Vector<double> trcPixel(2);
           trcPixel[0]=double((patchShape[0]+patchSize)/2-1);
           trcPixel[1]=double((patchShape[1]+patchSize)/2-1);
           ASKAPDEBUGASSERT((trcPixel[0]>0) && (trcPixel[1]>0));
           
           //std::cout<<trcPixel<<endl;
           
           csPatch.toWorld(world,trcPixel);
           csFull.toPixel(trcPixel,world);
           
           //std::cout<<trcPixel<<endl;
           
           for (size_t dim=0;dim<2;++dim) {
                const int pix1 = int(blcPixel[dim]);
                const int pix2 = int(trcPixel[dim]);
                blc[dim] = pix1>pix2 ? pix2 : pix1;
                trc[dim] = pix1>pix2 ? pix1 : pix2;
           }                                                    
        }
      
        /// @brief a helper method to make a parameter representing a test faceted image
        /// @details
        /// @param[in] params parameter container
        /// @param[in] name name of the parameter
        /// @param[in] nfacets number of facets
        /// @param[in] facetstep step in pixels between facet centres          
        static void makeParameter(askap::scimath::Params &params, const std::string &name,
                           const int nfacets, const int facetstep = 256)  
        {
           std::vector<std::string> direction(3);
           direction[0]="12h30m00.0";
           direction[1]="-15.00.00.00";
           direction[2]="J2000";
           std::vector<int> shape(2,256);
           std::vector<std::string> cellsize(2,"8arcsec");
           casa::Vector<casa::Stokes::StokesTypes> stokes(1, casa::Stokes::I);
           if (nfacets != 1) {
               SynthesisParamsHelper::add(params,name,direction,cellsize,shape,false,1.4e9,
                                      1.4e9,1,stokes, nfacets,facetstep);
           } else {
               SynthesisParamsHelper::add(params,name,direction,cellsize,shape,false,1.4e9,
                                      1.4e9,1,stokes);
           }
        }
        
   };
    
  } // namespace synthesis
} // namespace askap

#endif // #ifndef SYNTHESIS_PARAMS_HELPER_TEST_H

