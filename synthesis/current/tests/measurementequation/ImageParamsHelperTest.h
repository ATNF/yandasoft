/// @file
/// 
/// @brief Unit tests for ImageParamsHelper.
/// @details ImageParamsHelper class simplifies parsing the image parameter name in the
/// complex cases such as faceting and multi-frequency decomposition.
/// @note This class is also used inside one of the unit tests for SynthesisParamsHelper,
/// but is not tested comprehensively there.
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

#ifndef IMAGE_PARAMS_HELPER_TEST_H
#define IMAGE_PARAMS_HELPER_TEST_H

#include <measurementequation/ImageParamsHelper.h>

namespace askap
{
  namespace synthesis
  {
    
    class ImageParamsHelperTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(ImageParamsHelperTest);
      CPPUNIT_TEST(testVoidParse);
      CPPUNIT_TEST(testParseFacet);
      CPPUNIT_TEST(testExplicitFacet);
      CPPUNIT_TEST(testExplicitTaylorTerm);
      CPPUNIT_TEST(testParseTaylorTerm);
      CPPUNIT_TEST(testParseFacetTaylorTerm);
      CPPUNIT_TEST(testExplicitFacetTaylorTerm);
      CPPUNIT_TEST(testPartialNames);
      CPPUNIT_TEST(testMakeFacet);
      CPPUNIT_TEST(testMakeTaylor);
      CPPUNIT_TEST_SUITE_END();
    protected:
    
       /// @brief check correct settings for faceted parameter 
       /// @param[in] iph class to test
       void checkFacet(const ImageParamsHelper &iph) {
          CPPUNIT_ASSERT(iph.isValid());
          CPPUNIT_ASSERT(iph.isFacet());
          CPPUNIT_ASSERT(iph.name() == "image.test");
          CPPUNIT_ASSERT(iph.facetX() == 1);
          CPPUNIT_ASSERT(iph.facetY() == 2);
       }  

       /// @brief check correct settings for faceted parameter 
       /// @param[in] iph class to test
       void checkTaylor(const ImageParamsHelper &iph) {
          CPPUNIT_ASSERT(iph.isValid());
          CPPUNIT_ASSERT(iph.isTaylorTerm());
          CPPUNIT_ASSERT(iph.name() == "image.test");
          CPPUNIT_ASSERT(iph.order() == 3);
       }       
       
    public:
       void testVoidParse() {
          ImageParamsHelper iph("image.cmp.test");
          CPPUNIT_ASSERT(iph.isValid());
          CPPUNIT_ASSERT(!iph.isFacet());
          CPPUNIT_ASSERT(iph.paramName() == iph.name());
          CPPUNIT_ASSERT(iph.paramName() == "image.cmp.test");
          CPPUNIT_ASSERT_EQUAL(size_t(0),iph.suffix().size());
       }
       
       void testParseFacet() {
          ImageParamsHelper iph("image.test.facet.1.2");
          CPPUNIT_ASSERT(!iph.isTaylorTerm());
          checkFacet(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.facet.1.2");
          CPPUNIT_ASSERT_EQUAL(std::string(".facet.1.2"), iph.suffix());
       }

       void testExplicitFacet() {                    
          ImageParamsHelper iph("image.test",1,2);
          CPPUNIT_ASSERT(!iph.isTaylorTerm());
          checkFacet(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.facet.1.2");
          CPPUNIT_ASSERT_EQUAL(std::string(".facet.1.2"), iph.suffix());
       }
       
       void testParseTaylorTerm() {
          ImageParamsHelper iph("image.test.taylor.3");
          CPPUNIT_ASSERT(!iph.isFacet());
          checkTaylor(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.taylor.3");
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.3"), iph.suffix());          
       }
       
       void testExplicitTaylorTerm() {
          ImageParamsHelper iph("image.test",3);
          CPPUNIT_ASSERT(!iph.isFacet());
          checkTaylor(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.taylor.3");
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.3"), iph.suffix());          
       }
       
       void testParseFacetTaylorTerm() {
          ImageParamsHelper iph("image.test.taylor.3.facet.1.2");
          checkFacet(iph);
          checkTaylor(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.taylor.3.facet.1.2");          
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.3.facet.1.2"), iph.suffix());          
       }

       void testExplicitFacetTaylorTerm() {
          ImageParamsHelper iph("image.test",3,1,2);
          checkFacet(iph);
          checkTaylor(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.taylor.3.facet.1.2");          
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.3.facet.1.2"), iph.suffix());          
       }
       
       void testPartialNames() {
          ImageParamsHelper iph("image.test",3,1,2);
          CPPUNIT_ASSERT(iph.facetName() == "image.test.facet.1.2");
          CPPUNIT_ASSERT(iph.taylorName() == "image.test.taylor.3");
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.3.facet.1.2"), iph.suffix());          
       }
       
       void testMakeFacet() {
          ImageParamsHelper iph("image.test",3);
          CPPUNIT_ASSERT(!iph.isFacet());
          iph.makeFacet(1,2);
          checkFacet(iph);
          checkTaylor(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.taylor.3.facet.1.2");
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.3.facet.1.2"), iph.suffix());          
          
          ImageParamsHelper iph2("planeimage");
          CPPUNIT_ASSERT(!iph2.isFacet());
          CPPUNIT_ASSERT(!iph2.isTaylorTerm());
          iph2.makeFacet(3,0);
          CPPUNIT_ASSERT(iph2.paramName() == "planeimage.facet.3.0");
          CPPUNIT_ASSERT_EQUAL(std::string(".facet.3.0"), iph2.suffix());                    
          CPPUNIT_ASSERT(iph2.isFacet());
          CPPUNIT_ASSERT(!iph2.isTaylorTerm());          
       }
       
       void testMakeTaylor() {
          ImageParamsHelper iph("image.test",1,2);
          CPPUNIT_ASSERT(!iph.isTaylorTerm());
          iph.makeTaylorTerm(3);
          checkFacet(iph);
          checkTaylor(iph);
          CPPUNIT_ASSERT(iph.paramName() == "image.test.taylor.3.facet.1.2");          
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.3.facet.1.2"), iph.suffix());                    

          ImageParamsHelper iph2("planeimage");
          CPPUNIT_ASSERT(!iph2.isFacet());
          CPPUNIT_ASSERT(!iph2.isTaylorTerm());
          iph2.makeTaylorTerm(0);
          CPPUNIT_ASSERT(iph2.paramName() == "planeimage.taylor.0");
          CPPUNIT_ASSERT_EQUAL(std::string(".taylor.0"), iph2.suffix());          
          CPPUNIT_ASSERT(!iph2.isFacet());
          CPPUNIT_ASSERT(iph2.isTaylorTerm());
       }
       

    }; // class ImageParamsHelperTest
    
  } // namespace synthesis

} // namespace askap

#endif // #ifndef IMAGE_PARAMS_HELPER_TEST_H

