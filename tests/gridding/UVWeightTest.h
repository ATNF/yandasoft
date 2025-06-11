/// @file
///
/// Unit test for UVWeight and related classes (used in traditional weighting)
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

#include <askap/gridding/UVWeightCollection.h>
#include <askap/gridding/GenericUVWeightIndexTranslator.h>
#include <askap/gridding/GenericUVWeightAccessor.h>
#include <askap/gridding/GenericUVWeightBuilder.h>
#include <askap/gridding/UVWeightIndexTranslationHelper.h>
#include <askap/scimath/utils/EstimatorAdapter.h>

#include <cppunit/extensions/HelperMacros.h>

#include <casacore/casa/Arrays/IPosition.h>
#include <askap/askap/AskapUtil.h>

#include <Blob/BlobString.h>
#include <Blob/BlobOBufString.h>
#include <Blob/BlobIBufString.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>

namespace askap {

namespace synthesis {

class UVWeightTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(UVWeightTest);
   CPPUNIT_TEST(testUVWeight);
   CPPUNIT_TEST(testUVWeightFromCube);
   CPPUNIT_TEST(testUVWeightCollection);
   CPPUNIT_TEST_EXCEPTION(testUVWeightCollectionBadIndex, AskapError);
   CPPUNIT_TEST_EXCEPTION(testUVWeightCollectionBadIndexConstAccess, AskapError);
   CPPUNIT_TEST(testUVWeightCollectionMerge);
   CPPUNIT_TEST(testUVWeightCollectionIndices);
   CPPUNIT_TEST(testUVWeightCollectionSerialisation);
   CPPUNIT_TEST(testIndexTranslation);
   CPPUNIT_TEST(testGenericUVWeightAccessor);
   CPPUNIT_TEST(testGenericUVWeightAccessorViaPtr);
   CPPUNIT_TEST(testGenericUVWeightBuilder);
   CPPUNIT_TEST(testGenericUVWeightBuilderMerge);
   CPPUNIT_TEST(testGenericUVWeightBuilderFinalise);
   CPPUNIT_TEST(testGenericUVWeightBuilderSerialisation);
   CPPUNIT_TEST_EXCEPTION(testGenericUVWeightBuilderUninitialised, AskapError);
   CPPUNIT_TEST(testGenericUVWeightBuilderMergeViaAdapter);
   CPPUNIT_TEST_SUITE_END();
public:
   
   void testUVWeight() {
      UVWeight wtGrid1(10u,15u,1u);
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid1.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid1.vSize());
      CPPUNIT_ASSERT_EQUAL(1u, wtGrid1.nPlane());
      CPPUNIT_ASSERT(!wtGrid1.empty());

      // reference semantics in copy constructor
      UVWeight wtGrid2(wtGrid1);
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid2.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid2.vSize());
      CPPUNIT_ASSERT_EQUAL(1u, wtGrid2.nPlane());
      CPPUNIT_ASSERT(!wtGrid2.empty());

      UVWeight wtGrid3;
      CPPUNIT_ASSERT(wtGrid3.empty());
      // reference semantics in assignment operator
      wtGrid3 = wtGrid1;
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid3.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid3.vSize());
      CPPUNIT_ASSERT_EQUAL(1u, wtGrid3.nPlane());
      CPPUNIT_ASSERT(!wtGrid3.empty());
      
      // the following assignment should propagate to both weight grids (reference semantics)
      wtGrid1(5,7,0) = 1.5;
      for (casacore::uInt u = 0; u < wtGrid1.uSize(); ++u) {
           for (casacore::uInt v = 0; v < wtGrid1.vSize(); ++v) {
                const float expected = ( (u == 5) && (v == 7) ) ? 1.5 : 0.;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid1(u,v,0u), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid2(u,v,0u), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid3(u,v,0u), 1e-6);
           }
      }
   }

   void testUVWeightFromCube() {
      casacore::Cube<float> wtCube(10,15,2, 1.);
      UVWeight wtGrid(wtCube);
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid.vSize());
      CPPUNIT_ASSERT_EQUAL(2u, wtGrid.nPlane());
      CPPUNIT_ASSERT(!wtGrid.empty());
 
      // check reference semantics
      wtCube(5,7,0) = 0.5;
      for (casacore::uInt u = 0; u < wtGrid.uSize(); ++u) {
           for (casacore::uInt v = 0; v < wtGrid.vSize(); ++v) {
                for (casacore::uInt plane = 0; plane < wtGrid.nPlane(); ++plane) {
                     const float expected = ( (u == 5) && (v == 7)  && (plane == 0) ) ? 0.5 : 1.;
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid(u,v,plane), 1e-6);
                }
           }
      }
   }
   
   
   void testUVWeightCollection() {
       UVWeightCollection collection;
       collection.add(3u, 5u,7u, 1u);

       casacore::Cube<float> wtCube(10,15,1, 1.);
       collection.add(1u, wtCube);

       // check reference semantics
       wtCube(5,7,0) = 0.5;
       // test sparse indices
       for (casa::uInt index = 0; index < 4u; ++index) {
            CPPUNIT_ASSERT(collection.exists(index) == (index % 2 == 1));
       }
       // test shapes (this relies on the non-const interface)
       CPPUNIT_ASSERT(collection.get(1u).shape() == casacore::IPosition(3, 10, 15, 1));
       CPPUNIT_ASSERT(collection.get(3u).shape() == casacore::IPosition(3, 5, 7, 1));
  
       // test values
       for (casacore::uInt u = 0; u < wtCube.nrow(); ++u) {
           for (casacore::uInt v = 0; v < wtCube.ncolumn(); ++v) {
                if ((u < 5u) && (v<7u)) {
                    const float value = collection.get(3u)(u, v, 0u);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., value, 1e-6);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., getWeightViaConstInterface(collection, 3u, u,v, 0u), 1e-6);
                }
                const float expected = (u == 5) && (v == 7) ? 0.5 : 1.;
                const float actual = collection.get(1u)(u, v, 0u);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, actual, 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, getWeightViaConstInterface(collection, 1u, u,v, 0u), 1e-6);
           }
       }
   }

   void testUVWeightCollectionBadIndex() {
       UVWeightCollection collection;
       collection.get(0u);
   }

   void testUVWeightCollectionBadIndexConstAccess() {
       UVWeightCollection collection;
       collection.add(3u, 5u,7u, 1u);
       getWeightViaConstInterface(collection, 1u, 0u, 0u, 0u);
   }
  
   void testUVWeightCollectionMerge() {
       UVWeightCollection collection;
       collection.add(3u, 5u,7u, 1u);

       casacore::Cube<float> wtCube(10,15,1, 1.);
       collection.add(1u, wtCube);
        
       // second collection
       UVWeightCollection collection2;
       casacore::Cube<float> wtCube2(10,15,1, 2.);
       collection2.add(1u, wtCube2);

       casacore::Cube<float> wtCube3(3,5,1, 3.);
       collection2.add(2u, wtCube3);

       // merge
       collection.merge(collection2);

       // set values to check the behaviour w.r.t. the reference semantics, where applicable
       wtCube(5,7,0) = 0.5;
       wtCube2(5,7,0) = 0.25;
       wtCube3(2,3,0) = 0.125;

       // test sparse indices
       for (casa::uInt index = 0; index < 4u; ++index) {
            CPPUNIT_ASSERT(collection.exists(index) == (index != 0));
            CPPUNIT_ASSERT(collection2.exists(index) == ((index != 0) && (index != 3)));
       }
       // test shapes (this relies on the non-const interface)
       CPPUNIT_ASSERT(collection.get(1u).shape() == casacore::IPosition(3, 10, 15, 1));
       CPPUNIT_ASSERT(collection.get(2u).shape() == casacore::IPosition(3, 3, 5, 1));
       CPPUNIT_ASSERT(collection.get(3u).shape() == casacore::IPosition(3, 5, 7, 1));

       // test values
       for (casacore::uInt u = 0; u < wtCube.nrow(); ++u) {
           for (casacore::uInt v = 0; v < wtCube.ncolumn(); ++v) {
                if ((u < 5u) && (v<7u)) {
                    const float value = collection.get(3u)(u, v, 0u);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., value, 1e-6);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., getWeightViaConstInterface(collection, 3u, u,v, 0u), 1e-6);
                }

                if ((u < 3u) && (v<5u)) {
                    const float value = collection.get(2u)(u, v, 0u);
                    const float expected = (u == 2) && (v == 3) ? 0.125 : 3.;
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, value, 1e-6);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, getWeightViaConstInterface(collection, 2u, u,v, 0u), 1e-6);
                }
                
                const float expectedAfterMerge = (u == 5) && (v == 7) ? 0.5 : 3.;
                const float actualAfterMerge = collection.get(1u)(u, v, 0u);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedAfterMerge, actualAfterMerge, 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedAfterMerge, getWeightViaConstInterface(collection, 1u, u,v, 0u), 1e-6);
                const float expectedSrc = (u == 5) && (v == 7) ? 0.25 : 2.;
                const float actualSrc = collection2.get(1u)(u, v, 0u);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedSrc, actualSrc, 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedSrc, getWeightViaConstInterface(collection2, 1u, u,v, 0u), 1e-6);
           }
       }
  
   }

   void testUVWeightCollectionSerialisation() {
       UVWeightCollection collection;
       collection.add(3u, 5u,7u, 1u);

       casacore::Cube<float> wtCube(10,15,1, 1.);
       collection.add(1u, wtCube);
        
       // second collection, make it non-empty to test that the old content is removed
       UVWeightCollection collection2;
       casacore::Cube<float> wtCube2(10,15,1, 2.);
       collection2.add(1u, wtCube2);

       // now serialise the first collection into the blob and deserialise into the second one
       {
          LOFAR::BlobString b1(false);
          LOFAR::BlobOBufString bob(b1);
          LOFAR::BlobOStream bos(bob);
          collection.writeToBlob(bos);

          LOFAR::BlobIBufString bib(b1);
          LOFAR::BlobIStream bis(bib);
          collection2.readFromBlob(bis);
       }
       // reset the original collection (as a test, this also would test that we don't have any reference semantics causing problems)
       collection.clear();
       CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(0u), collection.indices().size());

       // now test the content of the second collection
       CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(2u), collection2.indices().size());
       CPPUNIT_ASSERT(collection2.exists(3u));
       CPPUNIT_ASSERT(collection2.exists(1u));

       CPPUNIT_ASSERT(collection2.get(1u).shape() == casacore::IPosition(3, 10, 15, 1));
       CPPUNIT_ASSERT(collection2.get(3u).shape() == casacore::IPosition(3, 5, 7, 1));

       for (casacore::uInt u = 0; u < wtCube.nrow(); ++u) {
           for (casacore::uInt v = 0; v < wtCube.ncolumn(); ++v) {
                if ((u < 5u) && (v<7u)) {
                    const float value = collection2.get(3u)(u, v, 0u);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., value, 1e-6);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., getWeightViaConstInterface(collection2, 3u, u,v, 0u), 1e-6);
                }
                const float value = collection2.get(1u)(u, v, 0u);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1., value, 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1., getWeightViaConstInterface(collection2, 1u, u,v, 0u), 1e-6);
           }
       }
   }

   void testUVWeightCollectionIndices() {
       UVWeightCollection collection;
       collection.add(3u, 5u,7u, 1u);
       collection.add(1u, 3u, 10u, 1u);
       std::set<casacore::uInt> indices = collection.indices();
       // test sparse indices
       for (casa::uInt index = 0; index < 4u; ++index) {
            CPPUNIT_ASSERT(collection.exists(index) == (index % 2 == 1));
            CPPUNIT_ASSERT((indices.find(index) == indices.end()) == (index % 2 == 0));
       }
   }

   void testIndexTranslation() {
       GenericUVWeightIndexTranslator ttor(1u, 36u, 36u*100u); 
       boost::shared_ptr<GenericUVWeightIndexTranslator> pTtor(&ttor, utility::NullDeleter());
       UVWeightIndexTranslationHelper<VoidClass> hlp(pTtor);

       for (casacore::uInt beam = 0; beam < 36u; ++beam) {
            for (casacore::uInt field = 0; field < 2u; ++field) {
                 for (casacore::uInt facet = 0; facet < 36u; ++facet) {
                      const casacore::uInt expected = beam + 36u*(field + 100u*facet);
                      CPPUNIT_ASSERT_EQUAL(expected, ttor.indexOf(beam, field, facet));
                      CPPUNIT_ASSERT_EQUAL(expected, hlp.indexOf(beam, field, facet));
                 }
            }
       }

       boost::shared_ptr<GenericUVWeightIndexTranslator> pOtherTtor(new GenericUVWeightIndexTranslator(1u));
       hlp.setTranslator(pOtherTtor);

       for (casacore::uInt beam = 0; beam < 36u; ++beam) {
            for (casacore::uInt field = 0; field < 2u; ++field) {
                 for (casacore::uInt facet = 0; facet < 36u; ++facet) {
                      const casacore::uInt oldExpected = beam + 36u*(field + 100u*facet);
                      CPPUNIT_ASSERT_EQUAL(oldExpected, ttor.indexOf(beam, field, facet));
                      CPPUNIT_ASSERT_EQUAL(beam, pOtherTtor->indexOf(beam, field, facet));
                      CPPUNIT_ASSERT_EQUAL(beam, hlp.indexOf(beam, field, facet));
                 }
            }
       }
   }

   void testGenericUVWeightAccessor() {
       casacore::Cube<float> wtCube(10,15,1, 1.);
       UVWeightCollection collection;
       collection.add(0u, wtCube);

       // check reference semantics
       wtCube(5,7,0) = 0.5;
       const boost::shared_ptr<GenericUVWeightAccessor> genericAcc(new GenericUVWeightAccessor(collection, 0u, 0u, 0u));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightAccessor> acc = genericAcc;
       CPPUNIT_ASSERT(acc);

       // get UVWeight object via the accessor interface for some random beam, field and source/facet (to check translation)
       const UVWeight wt = acc->getWeight(35u, 2u, 1u);
       // check shape
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.nrow()), wt.uSize());
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.ncolumn()), wt.vSize());
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.nplane()), wt.nPlane());
       // test values
       for (casacore::uInt u = 0; u < wtCube.nrow(); ++u) {
           for (casacore::uInt v = 0; v < wtCube.ncolumn(); ++v) {
                const float expected = (u == 5) && (v == 7) ? 0.5 : 1.;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wt(u, v, 0u), 1e-6);
           }
       }
       
   }

   void testGenericUVWeightAccessorViaPtr() {
       casacore::Cube<float> wtCube(10,15,1, 1.);
       boost::shared_ptr<UVWeightCollection> collection(new UVWeightCollection());
       collection->add(0u, wtCube);

       // check reference semantics
       wtCube(5,7,0) = 0.5;
       const boost::shared_ptr<GenericUVWeightIndexTranslator> translator(new GenericUVWeightIndexTranslator(0u, 0u, 0u));
       const boost::shared_ptr<GenericUVWeightAccessor> genericAcc(new GenericUVWeightAccessor(collection, translator));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightAccessor> acc = genericAcc;
       CPPUNIT_ASSERT(acc);

       // get UVWeight object via the accessor interface for some random beam, field and source/facet (to check translation)
       const UVWeight wt = acc->getWeight(35u, 2u, 1u);
       // check shape
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.nrow()), wt.uSize());
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.ncolumn()), wt.vSize());
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.nplane()), wt.nPlane());
       // test values
       for (casacore::uInt u = 0; u < wtCube.nrow(); ++u) {
           for (casacore::uInt v = 0; v < wtCube.ncolumn(); ++v) {
                const float expected = (u == 5) && (v == 7) ? 0.5 : 1.;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wt(u, v, 0u), 1e-6);
           }
       }
   }

   void testGenericUVWeightBuilder() {
       const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder(new GenericUVWeightBuilder(0u, 0u, 0u));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightBuilder> builder = genericBuilder;
       CPPUNIT_ASSERT(builder);
       builder->initialise(10u,15u,1u);
       // get UVWeight object via the builder interface for some arbitrary beam, field and source/facet (to check translation)
       UVWeight wt = builder->addWeight(35u, 2u, 1u);
       CPPUNIT_ASSERT(!wt.empty());
       // check that dimensions are as setup by the initialise method
       CPPUNIT_ASSERT_EQUAL(10u, wt.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt.nPlane());
       wt(5,7,0) += 0.5;
       // getUVWeight object for different (but still rather arbitrary) beam, field and source/facet to check that it maps to the same grid
       UVWeight wt2 = builder->addWeight(5u, 1u, 2u);
       // check that dimensions are as setup by the initialise method
       CPPUNIT_ASSERT_EQUAL(10u, wt2.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt2.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt2.nPlane());
       wt2(5,7,0) += 0.5;

       // test values 
       for (casacore::uInt u = 0; u < wt.uSize(); ++u) {
           for (casacore::uInt v = 0; v < wt.vSize(); ++v) {
                const float expected = (u == 5) && (v == 7) ? 1. : 0.;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wt(u, v, 0u), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wt2(u, v, 0u), 1e-6);
           }
       }
   }

   void testGenericUVWeightBuilderMerge() {
       // first builder is beam aware
       const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder1(new GenericUVWeightBuilder(1u, 0u, 0u));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightBuilder> builder1 = genericBuilder1;
       CPPUNIT_ASSERT(builder1);
       builder1->initialise(10u,15u,1u);

       // another builder, but this one puts everything on the same grid
       const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder2(new GenericUVWeightBuilder(0u, 0u, 0u));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightBuilder> builder2 = genericBuilder2;
       CPPUNIT_ASSERT(builder2);
       builder2->initialise(10u,15u,1u);

       // get UVWeight object via the builder interface for some arbitrary field and source/facet, beam index will translate directly into weight grid index
       UVWeight wt1 = builder1->addWeight(0u, 2u, 1u);
       // check that dimensions are as setup by the initialise method
       CPPUNIT_ASSERT_EQUAL(10u, wt1.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt1.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt1.nPlane());
       wt1(5,7,0) += 0.5;
       // same for different beam, this should refer to an independent grid
       UVWeight wt2 = builder1->addWeight(5u, 2u, 1u);
       CPPUNIT_ASSERT_EQUAL(10u, wt2.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt2.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt2.nPlane());
       wt2(5,7,0) += 1.5;
       // now use the other builder which ignores all indices and uses index = 0 for all
       UVWeight wt3 = builder2->addWeight(35u, 2u, 1u);
       CPPUNIT_ASSERT_EQUAL(10u, wt3.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt3.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt3.nPlane());
       wt3(5,7,0) += 2.;
       // now merge 
       builder1->merge(*builder2);

       // check the results (use read/write interface as one day we may remove read-only one for builders)
 
       // access beam=0 grid (note, builder1 is beam-aware)
       UVWeight wt4 = builder1->addWeight(0u, 3u, 5u);
       // access beam=5 grid
       UVWeight wt5 = builder1->addWeight(5u, 5u, 3u);
       CPPUNIT_ASSERT_EQUAL(wt4.uSize(), wt5.uSize());
       CPPUNIT_ASSERT_EQUAL(wt4.vSize(), wt5.vSize());
       CPPUNIT_ASSERT_EQUAL(wt4.nPlane(), wt5.nPlane());
       CPPUNIT_ASSERT_EQUAL(10u, wt4.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt4.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt4.nPlane());

       for (casacore::uInt u = 0; u < wt4.uSize(); ++u) {
           for (casacore::uInt v = 0; v < wt4.vSize(); ++v) {
                const float expected4 = (u == 5) && (v == 7) ? 2.5 : 0.;
                const float expected5 = (u == 5) && (v == 7) ? 1.5 : 0.;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4, wt4(u, v, 0u), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5, wt5(u, v, 0u), 1e-6);
           }
       }
   }

   void testGenericUVWeightBuilderFinalise() {
       // builder is beam aware, use it to set up two weight grids
       const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder(new GenericUVWeightBuilder(1u, 0u, 0u));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightBuilder> builder = genericBuilder;
       CPPUNIT_ASSERT(builder);
       builder->initialise(10u,15u,1u);
       // don't need the result, the call will setup new grid with zero weights with flat index=35u.
       builder->addWeight(35u, 2u, 1u);
       // same for different index
       builder->addWeight(5u, 2u, 1u);
  
       TestFinaliseClass tfc;
       UVWeightCollection& collection = builder->finalise(tfc);

       // check results via the collection reference
       std::set<casacore::uInt> indices = collection.indices();
       CPPUNIT_ASSERT(indices.find(5u) != indices.end());
       CPPUNIT_ASSERT(indices.find(35u) != indices.end());
       CPPUNIT_ASSERT_EQUAL(size_t(2u), indices.size());
       const UVWeight wt1 = collection.get(5u);
       const UVWeight wt2 = collection.get(35u);

       CPPUNIT_ASSERT_EQUAL(wt1.uSize(), wt2.uSize());
       CPPUNIT_ASSERT_EQUAL(wt1.vSize(), wt2.vSize());
       CPPUNIT_ASSERT_EQUAL(wt1.nPlane(), wt2.nPlane());
       CPPUNIT_ASSERT_EQUAL(10u, wt1.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt1.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt1.nPlane());

       for (casacore::uInt u = 0; u < wt1.uSize(); ++u) {
           for (casacore::uInt v = 0; v < wt1.vSize(); ++v) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(2., wt1(u, v, 0u), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(2., wt2(u, v, 0u), 1e-6);
           }
       }
       
   }

   void testGenericUVWeightBuilderUninitialised() {
       // setup the builder, parameters don't matter really
       const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder(new GenericUVWeightBuilder(0u, 0u, 0u));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightBuilder> builder = genericBuilder;
       CPPUNIT_ASSERT(builder);
       // don't care about the result, but calling addWeight for uninitialised builder which doesn't have the appropriate grid (e.g. as a result of
       // merge) should cause an exception
       builder->addWeight(35u, 2u, 1u);
   }

   void testGenericUVWeightBuilderSerialisation() {
       GenericUVWeightBuilder genericBuilder(1u, 2u, 3u);
       genericBuilder.initialise(10u,15u,1u);
       // get UVWeight object via the builder interface for some arbitrary beam, field and source/facet (to check translation)
       UVWeight wt = genericBuilder.addWeight(35u, 2u, 1u);
       CPPUNIT_ASSERT(!wt.empty());
       // check that dimensions are as setup by the initialise method
       CPPUNIT_ASSERT_EQUAL(10u, wt.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt.nPlane());
       wt(5,7,0) += 0.5;
       // another builder to be replaced in deserialisation, with some different parameters
       GenericUVWeightBuilder builder2(10u, 20u, 30u);
       // now serialise the first one and deserialise into the second one   
       {
          LOFAR::BlobString b1(false);
          LOFAR::BlobOBufString bob(b1);
          LOFAR::BlobOStream bos(bob);
          genericBuilder.writeToBlob(bos);

          LOFAR::BlobIBufString bib(b1);
          LOFAR::BlobIStream bis(bib);
          builder2.readFromBlob(bis);
       }

       // check that the index translation is as in the original builder
       CPPUNIT_ASSERT_EQUAL(1u, builder2.indexOf(1u, 0u, 0u));
       CPPUNIT_ASSERT_EQUAL(2u, builder2.indexOf(0u, 1u, 0u));
       CPPUNIT_ASSERT_EQUAL(3u, builder2.indexOf(0u, 0u, 1u));

       // add weight with the new index to test that parameters passed to initialise method of the original builder still apply
       // (i.e. the brand new weight object was initialised with correct dimensions)
       const UVWeight wt2 = builder2.addWeight(5u, 0u, 0u);
       CPPUNIT_ASSERT(!wt2.empty());
       CPPUNIT_ASSERT_EQUAL(10u, wt2.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt2.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt2.nPlane());

       // get weight for indices setup with the original builder and test values
       const UVWeight wt3 = genericBuilder.addWeight(35u, 2u, 1u);
       CPPUNIT_ASSERT(!wt3.empty());
       CPPUNIT_ASSERT_EQUAL(10u, wt3.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt3.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt3.nPlane());

       for (casacore::uInt u = 0; u < wt3.uSize(); ++u) {
            for (casacore::uInt v = 0; v < wt3.vSize(); ++v) {
                 const float value = (u == 5u) && (v == 7u) ? 0.5f : 0.f;
                 CPPUNIT_ASSERT_DOUBLES_EQUAL(value, wt3(u,v,0u), 1e-6);
            }
       }
   }

   /// @brief test of merge operation via EstimatorAdapter
   /// @note The actual test is very basic (because the merge method has already been tested), do it largely to test
   /// compilation of the adapter with the given template argument.
   void testGenericUVWeightBuilderMergeViaAdapter() {
       const boost::shared_ptr<GenericUVWeightBuilder> builder1(new GenericUVWeightBuilder(0u, 0u, 0u));
       scimath::EstimatorAdapter<GenericUVWeightBuilder> adapter1(builder1);
       builder1->initialise(10u,15u,1u);
       // get UVWeight object via the builder interface for some arbitrary beam, field and source/facet (translated to index zero anyway)
       UVWeight wt = builder1->addWeight(35u, 2u, 1u);
       CPPUNIT_ASSERT(!wt.empty());
       CPPUNIT_ASSERT_EQUAL(10u, wt.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt.nPlane());
       wt(5,7,0) += 0.5;
       // now setup a second builder and adapter
       const boost::shared_ptr<GenericUVWeightBuilder> builder2(new GenericUVWeightBuilder(0u, 0u, 0u));
       builder2->initialise(10u,15u,1u);
       // get UVWeight object via the builder interface for some arbitrary beam, field and source/facet (translated to index zero anyway)
       UVWeight wt2 = builder2->addWeight(5u, 12u, 2u);
       CPPUNIT_ASSERT(!wt2.empty());
       CPPUNIT_ASSERT_EQUAL(10u, wt2.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt2.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt2.nPlane());
       wt2(5,7,0) += 1.5;
       scimath::EstimatorAdapter<GenericUVWeightBuilder> adapter2(builder2);
       // cast to normal equations interface (just for sake of the test, not really needed) and perform merge
       const scimath::INormalEquations& ne1 = adapter1;
       scimath::INormalEquations& ne2 = adapter2;
       ne2.merge(ne1);
       // get fresh pointer out of the adapter just for test, although the original one could have been used too
       const boost::shared_ptr<GenericUVWeightBuilder> builder3 = dynamic_cast<scimath::EstimatorAdapter<GenericUVWeightBuilder>&>(ne2).get();
       CPPUNIT_ASSERT(builder3);
       // get the weight plane for index 0 and test the content
       const UVWeight wt3 = builder3->addWeight(0u, 0u, 0u);
       CPPUNIT_ASSERT(!wt3.empty());
       CPPUNIT_ASSERT_EQUAL(10u, wt3.uSize());
       CPPUNIT_ASSERT_EQUAL(15u, wt3.vSize());
       CPPUNIT_ASSERT_EQUAL(1u, wt3.nPlane());
       for (casacore::uInt u = 0; u < wt3.uSize(); ++u) {
            for (casacore::uInt v = 0; v < wt3.vSize(); ++v) {
                 const float value = (u == 5u) && (v == 7u) ? 2.f : 0.f;
                 CPPUNIT_ASSERT_DOUBLES_EQUAL(value, wt3(u,v,0u), 1e-6);
            }
       }
   }
   
protected:

   float getWeightViaConstInterface(const UVWeightCollection &collection, casacore::uInt index, casacore::uInt u, casacore::uInt v, casacore::uInt plane) {
      const UVWeight wtGrid = collection.get(index);
      return wtGrid(u,v,plane);
   }

   struct VoidClass {
   };

   struct TestFinaliseClass : public IUVWeightCalculator {
       /// @brief override method to do processing for the given weight 
       /// @param[in] wt weight to work with (it is modified in situ).
       virtual void process(casacore::Matrix<float> &wt) const override final {
           wt.set(2.);
       };
   };

};
    
} // namespace synthesis

} // namespace askap

