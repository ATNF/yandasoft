/// @file
/// 
/// @brief Unit tests for VisMetaDataStats class
/// @details VisMetaDataStats accumulates certain statistics of the visibility data.
/// It is used to provide advise capability for the parameters used in imager and
/// calibrator.
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

#ifndef VIS_META_DATA_STATS_ME_TEST_H
#define VIS_META_DATA_STATS_ME_TEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <measurementequation/VisMetaDataStats.h>
#include <dataaccess/DataAccessorStub.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVDirection.h>

#include <Blob/BlobString.h>
#include <Blob/BlobOBufString.h>
#include <Blob/BlobIBufString.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>


namespace askap
{
  namespace synthesis
  {

    class VisMetaDataStatsTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(VisMetaDataStatsTest);
      CPPUNIT_TEST(testProcessModified);
      CPPUNIT_TEST(testProcess);
      CPPUNIT_TEST(testMerge);
      CPPUNIT_TEST(testBlobStream);
      CPPUNIT_TEST(testSnapShot);  
      CPPUNIT_TEST_EXCEPTION(testTangentCheck,AskapError);    
      CPPUNIT_TEST_EXCEPTION(testToleranceCheck,AskapError);    
      CPPUNIT_TEST(testReset);  
      //CPPUNIT_TEST(testDirectionMerge);
      CPPUNIT_TEST(testDirOffsets);
      CPPUNIT_TEST_SUITE_END();
    protected:
      static void modifyStubbedData(accessors::DataAccessorStub &acc) {
         for (casa::uInt row=0; row<acc.nRow(); ++row) {
            ++acc.itsFeed1[row];
            ++acc.itsFeed2[row];
            ++acc.itsAntenna1[row];
            ++acc.itsAntenna2[row];
            for (casa::uInt dim=0; dim<3; ++dim) {
                 acc.itsUVW[row](dim) *= 10.;
            }
            acc.itsPointingDir1[row].shift(-0.001,0.001,casa::True);
            acc.itsPointingDir2[row].shift(-0.001,0.001,casa::True);            
         }
         for (casa::uInt chan=0; chan<acc.nChannel(); ++chan) {
              acc.itsFrequency[chan] += 10e6;
         }
      }  
    public:
      void testProcessModified() {
         accessors::DataAccessorStub acc(true);
         modifyStubbedData(acc);
         VisMetaDataStats stats;
         CPPUNIT_ASSERT_EQUAL(0ul, stats.nVis());
         stats.process(acc);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.41e9,stats.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.27e9,stats.minFreq(),1.);
         
         const double freqFactor = 1.41/1.4; // we have different maximum frequency now
         CPPUNIT_ASSERT_DOUBLES_EQUAL(73120.88*freqFactor,stats.maxU(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(47906.8*freqFactor,stats.maxV(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(72008.7*freqFactor,stats.maxW(),1.);
         CPPUNIT_ASSERT_EQUAL(31u, stats.nAntennas());
         CPPUNIT_ASSERT_EQUAL(2u, stats.nBeams());
         CPPUNIT_ASSERT_EQUAL(3480ul, stats.nVis());
         // test directions for "offset" beam
         casa::MVDirection expectedDir(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         expectedDir.shift(-0.001,0.001,casa::True);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., expectedDir.separation(stats.centre()), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,stats.maxOffsets().first,1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,stats.maxOffsets().second,1e-6);         
      }

      void checkCombined(const VisMetaDataStats &stats) {
         // verify results after merge (we have the same end result in various tests)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.41e9,stats.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.260e9,stats.minFreq(),1.);
         
         const double freqFactor = 1.41/1.4; // we have different maximum frequency now
         CPPUNIT_ASSERT_DOUBLES_EQUAL(73120.88 * freqFactor,stats.maxU(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(47906.8 * freqFactor,stats.maxV(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(72008.7 * freqFactor,stats.maxW(),1.);
         CPPUNIT_ASSERT_EQUAL(31u, stats.nAntennas());
         CPPUNIT_ASSERT_EQUAL(2u, stats.nBeams());
         CPPUNIT_ASSERT_EQUAL(6960ul, stats.nVis());
         casa::MVDirection shiftedDir(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         shiftedDir.shift(-0.0005,0.0005,casa::True);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., shiftedDir.separation(stats.centre()), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0005,stats.maxOffsets().first,1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0005,stats.maxOffsets().second,1e-6);         
      }
      
      void testProcess() {
         accessors::DataAccessorStub acc(true);
         VisMetaDataStats stats;
         CPPUNIT_ASSERT_EQUAL(0ul, stats.nVis());
         stats.process(acc);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4e9,stats.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.260e9,stats.minFreq(),1.);

         // note, we didn't independently verify the following uvw 
         // values, but the magnitude make sense for the stubbed layout
         CPPUNIT_ASSERT_DOUBLES_EQUAL(7312.088,stats.maxU(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(4790.68,stats.maxV(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(7200.87,stats.maxW(),1.);
         bool noResidualW = false;
         try {
             stats.maxResidualW(); // this should fail!
         } 
         catch (const askap::AskapError &) {
             noResidualW = true;
         } 
         CPPUNIT_ASSERT(noResidualW);
         //
         CPPUNIT_ASSERT_EQUAL(30u, stats.nAntennas());
         CPPUNIT_ASSERT_EQUAL(1u, stats.nBeams());
         CPPUNIT_ASSERT_EQUAL(3480ul, stats.nVis());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,stats.maxOffsets().first,1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,stats.maxOffsets().second,1e-6);
         const casa::MVDirection expectedDir(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., expectedDir.separation(stats.centre()), 1e-6);
         // now change the data and accumulate again
         modifyStubbedData(acc);
         stats.process(acc);
         checkCombined(stats);
      } 
      
      void testMerge() {
         accessors::DataAccessorStub acc(true);
         const casa::MVDirection tangent(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         VisMetaDataStats stats1(tangent);
         stats1.process(acc);         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4e9,stats1.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.260e9,stats1.minFreq(),1.);
         modifyStubbedData(acc);
         VisMetaDataStats stats2(tangent);
         stats2.process(acc);         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.41e9,stats2.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.27e9,stats2.minFreq(),1.);
         // now merge
         stats1.merge(stats2);
         checkCombined(stats1);
         VisMetaDataStats stats3(tangent);
         stats1.merge(stats3);
         checkCombined(stats1);
         stats3.merge(stats1);
         checkCombined(stats3);         
      }
      
      void testBlobStream() {
         const casa::MVDirection tangent(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         accessors::DataAccessorStub acc(true);
         VisMetaDataStats stats1(tangent);
         stats1.process(acc);
         modifyStubbedData(acc);
         VisMetaDataStats stats2(tangent);
         stats2.process(acc);
         // serialise
         LOFAR::BlobString b1(false);
         LOFAR::BlobOBufString bob(b1);
         LOFAR::BlobOStream bos(bob);
         bos << stats2;
         // de-serialise
         LOFAR::BlobIBufString bib(b1);
         LOFAR::BlobIStream bis(bib);
         VisMetaDataStats stats3;
         bis >> stats3;
         // compare two versions         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.41e9,stats3.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.27e9,stats3.minFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(stats2.maxU(),stats3.maxU(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(stats2.maxV(),stats3.maxV(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(stats2.maxW(),stats3.maxW(),1.);
         CPPUNIT_ASSERT_EQUAL(stats2.nAntennas(),stats3.nAntennas());
         CPPUNIT_ASSERT_EQUAL(stats2.nBeams(),stats3.nBeams());
         CPPUNIT_ASSERT_EQUAL(stats2.nVis(),stats3.nVis());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,stats2.centre().separation(stats3.centre()),1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(stats2.maxOffsets().first,stats3.maxOffsets().first,1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(stats2.maxOffsets().second,stats3.maxOffsets().second,1e-6);
         // merge in a different way than before
         stats3.merge(stats1);
         checkCombined(stats3);
      }
      
      void testSnapShot() {
         accessors::DataAccessorStub acc(true);
         const casa::MVDirection tangent(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         VisMetaDataStats stats(tangent,1);
         stats.process(acc);         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4e9,stats.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.260e9,stats.minFreq(),1.);
         // the accessor stub doesn't do uvw-rotation, so the values are the same
         CPPUNIT_ASSERT_DOUBLES_EQUAL(7312.088,stats.maxU(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(4790.68,stats.maxV(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(7200.87,stats.maxW(),1.);
         // now can check residual w-term
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7953,stats.maxResidualW(),1e-4);         
         CPPUNIT_ASSERT_EQUAL(30u, stats.nAntennas());
         CPPUNIT_ASSERT_EQUAL(1u, stats.nBeams());
         CPPUNIT_ASSERT_EQUAL(3480ul, stats.nVis());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,stats.maxOffsets().first,1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,stats.maxOffsets().second,1e-6);
         const casa::MVDirection expectedDir(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., expectedDir.separation(stats.centre()), 1e-6);         
      }
      
      void testTangentCheck() {
         accessors::DataAccessorStub acc(true);
         const casa::MVDirection tangent(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         VisMetaDataStats stats(tangent);
         stats.process(acc);   
         VisMetaDataStats stats1;
         stats1.merge(stats);      
      }
      
      void testToleranceCheck() {
         accessors::DataAccessorStub acc(true);
         const casa::MVDirection tangent(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         VisMetaDataStats stats(tangent,1.);
         stats.process(acc);   
         VisMetaDataStats stats1(tangent,700.);
         stats1.merge(stats);
      }
      
      void testReset() {
         accessors::DataAccessorStub acc(true);
         const casa::MVDirection tangent(casa::Quantity(0, "deg"), casa::Quantity(0, "deg"));
         VisMetaDataStats stats1(tangent);
         stats1.process(acc);         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4e9,stats1.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.260e9,stats1.minFreq(),1.);
         VisMetaDataStats stats2(tangent);
         stats2.process(acc);         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4e9,stats2.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.260e9,stats2.minFreq(),1.);
         CPPUNIT_ASSERT_EQUAL(3480ul, stats2.nVis());
         stats2.reset();
         CPPUNIT_ASSERT_EQUAL(0ul, stats2.nVis());         
         modifyStubbedData(acc);
         stats2.process(acc);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.41e9,stats2.maxFreq(),1.);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.27e9,stats2.minFreq(),1.);
         CPPUNIT_ASSERT_EQUAL(3480ul, stats2.nVis());
         // now merge
         stats1.merge(stats2);
         checkCombined(stats1);
      }
   
      void testDirectionMerge() {
         // Unit test written while debugging ASKAPSDP-1689
         
         // set up the tree
         std::vector<boost::shared_ptr<VisMetaDataStats> > tree(7);
         for (size_t i = 0; i<tree.size(); ++i) {
              const casa::MVDirection tangent(asQuantity("12:30:00.00"), asQuantity("-045.00.00.00"));
              
              accessors::DataAccessorStub acc(true);
              const casa::MVDirection testDir1(asQuantity("12:35:39.36"), asQuantity("-044.59.28.59"));
              acc.itsPointingDir1.set(testDir1);
              acc.itsPointingDir2.set(testDir1);

              if (i == 5) {
                  const casa::MVDirection testDir2(asQuantity("12:35:39.34"), asQuantity("-044.59.28.59"));
                  acc.itsPointingDir1.set(testDir2);
                  acc.itsPointingDir2.set(testDir2);
              }
              
              VisMetaDataStats stats(tangent);
              stats.process(acc);         

              // serialise
              LOFAR::BlobString b1(false);
              LOFAR::BlobOBufString bob(b1);
              LOFAR::BlobOStream bos(bob);
              bos << stats;

              tree[i].reset(new VisMetaDataStats);
              CPPUNIT_ASSERT(tree[i]);

              // de-serialise
              LOFAR::BlobIBufString bib(b1);
              LOFAR::BlobIStream bis(bib);
              bis >> *tree[i];
         }

         // manual tree-reduction
         CPPUNIT_ASSERT_EQUAL(size_t(7u), tree.size());
         tree[1]->merge(*tree[3]);
         tree[1]->merge(*tree[4]);
         tree[2]->merge(*tree[5]);
         tree[2]->merge(*tree[6]);
         tree[0]->merge(*tree[1]);
         tree[0]->merge(*tree[2]);
      }

      void testDirOffsets() {
         const casa::MVDirection tangent(asQuantity("12:30:00.00"), asQuantity("-045.00.00.00"));
         VisMetaDataStats stats1(tangent);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., stats1.maxOffsets().first, 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., stats1.maxOffsets().second, 1e-6);

         
         const casa::MVDirection testDir1(asQuantity("12:35:39.36"), asQuantity("-044.59.28.59"));

         const std::pair<double, double> offsets1 = stats1.getOffsets(testDir1);
         const casa::MVDirection processedDir1 = stats1.getOffsetDir(offsets1);

         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., testDir1.separation(processedDir1), 1e-10);
      }
      
    };
  
  } // namespace synthesis

} // namespace askap

#endif // #ifndef VIS_META_DATA_STATS_ME_TEST_H

