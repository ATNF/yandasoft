/// @file
///
/// Unit test for UVWeightParamsHelper class 
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

#include <askap/gridding/UVWeightParamsHelper.h>
#include <askap/scimath/fitting/Params.h>
#include <cppunit/extensions/HelperMacros.h>

namespace askap {

namespace synthesis {

class UVWeightParamsHelperTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(UVWeightParamsHelperTest);
   CPPUNIT_TEST(testAdd);
   CPPUNIT_TEST(testAdd2);
   CPPUNIT_TEST_EXCEPTION(testDuplicateAdd, AskapError);
   CPPUNIT_TEST(testRemove);
   CPPUNIT_TEST_EXCEPTION(testAddRedundant, AskapError);
   CPPUNIT_TEST_EXCEPTION(testAddWrongIndex, AskapError);
   CPPUNIT_TEST(testGet);
   CPPUNIT_TEST(testGetIndexTranslator);
   CPPUNIT_TEST(testCopy);
   CPPUNIT_TEST_SUITE_END();
public:
   void setUp() {
      itsParams.reset(new scimath::Params());
   }
  
   void tearDown() {
      itsParams.reset();
   }

   void testAdd() {
      UVWeightParamsHelper hlp(itsParams);
      CPPUNIT_ASSERT(!hlp.exists("gc"));
      CPPUNIT_ASSERT(itsParams);

      UVWeightCollection collection;
      collection.add(3u, 5u,7u, 1u);

      casacore::Cube<float> wtCube(10,15,1, 1.);
      collection.add(1u, wtCube);

      boost::shared_ptr<GenericUVWeightIndexTranslator> ttor(new GenericUVWeightIndexTranslator(1u,2u,3u));

      hlp.addUVWeights("gc", collection, ttor);

      CPPUNIT_ASSERT(hlp.exists("gc"));
      // check results of addition directly in Params
      CPPUNIT_ASSERT(itsParams->has("uvweight.gc.3"));
      CPPUNIT_ASSERT(itsParams->has("uvweight.gc.1"));
      CPPUNIT_ASSERT(!itsParams->has("uvweight.gc.2"));
      CPPUNIT_ASSERT(itsParams->has("uvweight_indices.gc"));
      CPPUNIT_ASSERT(!itsParams->isFree("uvweight.gc.3"));
      CPPUNIT_ASSERT(!itsParams->isFree("uvweight.gc.1"));
      CPPUNIT_ASSERT(!itsParams->isFree("uvweight_indices.gc"));

      const casacore::Vector<float> translationIndices = itsParams->valueF("uvweight_indices.gc");
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3u), translationIndices.nelements());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.f, translationIndices[0], 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.f, translationIndices[1], 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3.f, translationIndices[2], 1e-5);

      const casacore::Matrix<float> wt1 = itsParams->valueF("uvweight.gc.3").nonDegenerate();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(5u), wt1.nrow());
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(7u), wt1.ncolumn());
      for (casacore::uInt row = 0; row < wt1.nrow(); ++row) {
           for (casacore::uInt col = 0; col < wt1.ncolumn(); ++col) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, wt1(row,col), 1e-5);
           }
      }

      const casacore::Matrix<float> wt2 = itsParams->valueF("uvweight.gc.1").nonDegenerate();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(10u), wt2.nrow());
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(15u), wt2.ncolumn());
      for (casacore::uInt row = 0; row < wt2.nrow(); ++row) {
           for (casacore::uInt col = 0; col < wt2.ncolumn(); ++col) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1.f, wt2(row,col), 1e-5);
           }
      }
   }

   void testAdd2() {
      // another addition test, but with trivial index mapping (i.e. all indices ignored)
      UVWeightParamsHelper hlp(itsParams);
      CPPUNIT_ASSERT(!hlp.exists("cena"));
      CPPUNIT_ASSERT(itsParams);

      UVWeightCollection collection;
      collection.add(0u, 5u,7u, 1u);
      hlp.addUVWeights("cena", collection);

      CPPUNIT_ASSERT(hlp.exists("cena"));
      CPPUNIT_ASSERT(itsParams->has("uvweight.cena.0"));
      CPPUNIT_ASSERT(!itsParams->has("uvweight_indices.cena"));
      CPPUNIT_ASSERT(!itsParams->isFree("uvweight.cena.0"));

      // check parameter value directly via Params
      const casacore::Matrix<float> wt = itsParams->valueF("uvweight.cena.0").nonDegenerate();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(5u), wt.nrow());
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(7u), wt.ncolumn());
      for (casacore::uInt row = 0; row < wt.nrow(); ++row) {
           for (casacore::uInt col = 0; col < wt.ncolumn(); ++col) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, wt(row,col), 1e-5);
           }
      }

      // check that the index translator (the trivial one) behaves as expected
      const boost::shared_ptr<IUVWeightIndexTranslator> ttor = hlp.getIndexTranslator("cena");
      CPPUNIT_ASSERT(ttor);
      CPPUNIT_ASSERT_EQUAL(0u, ttor->indexOf(10u,0u,0u));
      CPPUNIT_ASSERT_EQUAL(0u, ttor->indexOf(0u,10u,0u));
      CPPUNIT_ASSERT_EQUAL(0u, ttor->indexOf(0u,0u,10u));
   }

   void testDuplicateAdd() {
      testAdd();
      // now addition of a collection with the same name/indices should throw an exception
      // (itsParam is only reinitialised between tests). Note, we cannot just call testAdd
      // again here - it does some additional checks and they will trigger first. So do the
      // second addition from scratch
      UVWeightParamsHelper hlp(itsParams);
      UVWeightCollection collection;
      collection.add(3u, 5u,7u, 1u);

      boost::shared_ptr<GenericUVWeightIndexTranslator> ttor(new GenericUVWeightIndexTranslator(10u,20u,30u));
      // the next line should throw an exception
      hlp.addUVWeights("gc", collection, ttor);
   }

   void testRemove() {
      testAdd();
      testAdd2();
      // now remove some of the parameters just added (note, itsParam is reinitialised between the tests)
      UVWeightParamsHelper hlp(itsParams);
      hlp.remove("gc");
      CPPUNIT_ASSERT(hlp.exists("cena"));
      // now a new call to testAdd should succeed, the parameters added by testAdd2 should remain intact
      testAdd();
   }

   void testAddRedundant() {
      UVWeightParamsHelper hlp(itsParams);
      CPPUNIT_ASSERT(!hlp.exists("cena"));
      CPPUNIT_ASSERT(itsParams);

      UVWeightCollection collection;
      collection.add(0u, 5u,7u, 1u);
      // add inaccessible weight given the trivial index translator
      collection.add(1u, 5u,7u, 1u);
      // the following should throw an exception
      hlp.addUVWeights("cena", collection);
   }

   void testAddWrongIndex() {
      // similar test to testAddRedundant, but with single non-zero index
      UVWeightParamsHelper hlp(itsParams);
      CPPUNIT_ASSERT(!hlp.exists("cena"));
      CPPUNIT_ASSERT(itsParams);

      UVWeightCollection collection;
      collection.add(1u, 5u,7u, 1u);

      // the following should throw an exception
      hlp.addUVWeights("cena", collection);
   }
  
   void testGet() {
      testAdd();
      testAdd2();
      UVWeightParamsHelper hlp(itsParams);
      const boost::shared_ptr<UVWeightCollection> collectionPtr = hlp.getUVWeights("gc");
      checkCollectionFromTestAdd(collectionPtr);
   }

   void testGetIndexTranslator() {
      testAdd();
      UVWeightParamsHelper hlp(itsParams);

      // check that the index translator behaves as expected 
      // note, the trivial translator is checked in testAdd2
      const boost::shared_ptr<IUVWeightIndexTranslator> ttor = hlp.getIndexTranslator("gc");
      CPPUNIT_ASSERT(ttor);
      CPPUNIT_ASSERT_EQUAL(1u, ttor->indexOf(1u,0u,0u));
      CPPUNIT_ASSERT_EQUAL(2u, ttor->indexOf(0u,1u,0u));
      CPPUNIT_ASSERT_EQUAL(3u, ttor->indexOf(0u,0u,1u));
   }

   void testCopy() {
      testAdd();
      testAdd2();
      UVWeightParamsHelper hlp(itsParams);
      scimath::Params dest;
      hlp.copyTo(dest, "gc");
      // test an alternative constructor from reference as well
      UVWeightParamsHelper hlp2(dest);
      CPPUNIT_ASSERT(hlp2.exists("gc"));

      // check that the required stuff has been copied
      const boost::shared_ptr<IUVWeightIndexTranslator> ttor = hlp2.getIndexTranslator("gc");
      CPPUNIT_ASSERT(ttor);
      CPPUNIT_ASSERT_EQUAL(1u, ttor->indexOf(1u,0u,0u));
      CPPUNIT_ASSERT_EQUAL(2u, ttor->indexOf(0u,1u,0u));
      CPPUNIT_ASSERT_EQUAL(3u, ttor->indexOf(0u,0u,1u));

      const boost::shared_ptr<UVWeightCollection> collectionPtr = hlp2.getUVWeights("gc");
      checkCollectionFromTestAdd(collectionPtr);
   }

protected:
   // helper method to check collection to be what is expected after testAdd
   void checkCollectionFromTestAdd(const boost::shared_ptr<UVWeightCollection> &collectionPtr) {
      CPPUNIT_ASSERT(collectionPtr);
      // get non-const interface to get cube interface for access, although we use it read-only
      UVWeightCollection& collection = *collectionPtr;
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(2u), collection.indices().size());
      CPPUNIT_ASSERT(collection.exists(1u));
      CPPUNIT_ASSERT(collection.exists(3u));

      const casacore::Cube<float>& wt1 = collection.get(3u);
      CPPUNIT_ASSERT_EQUAL(casacore::IPosition(3,5,7,1), wt1.shape());
      for (casacore::uInt row = 0; row < wt1.nrow(); ++row) {
           for (casacore::uInt col = 0; col < wt1.ncolumn(); ++col) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, wt1(row, col, 0), 1e-5);
           }
      }

      const casacore::Cube<float>& wt2 = collection.get(1u);
      CPPUNIT_ASSERT_EQUAL(casacore::IPosition(3,10,15,1), wt2.shape());
      for (casacore::uInt row = 0; row < wt2.nrow(); ++row) {
           for (casacore::uInt col = 0; col < wt2.ncolumn(); ++col) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1.f, wt2(row, col, 0), 1e-5);
           }
      }
   }

private:
   boost::shared_ptr<scimath::Params> itsParams;
};
    
} // namespace synthesis

} // namespace askap

