/// @file
///
/// Unit test for the basis function
///
///
/// @copyright (c) 2007 CSIRO
/// Australia Telescope National Facility (ATNF)

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
/// @author Tim Cornwell <tim.cornwell@csiro.au>

#include <deconvolution/EntropyI.h>
#include <deconvolution/Emptiness.h>

#include <cppunit/extensions/HelperMacros.h>

#include <boost/shared_ptr.hpp>

using namespace casa;

namespace askap {
  
  namespace synthesis {
    
    class EntropyTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(EntropyTest);
      CPPUNIT_TEST(testEntropyI);
      CPPUNIT_TEST(testEmptiness);
      CPPUNIT_TEST_SUITE_END();
    public:
      
      void setUp() {
        itsModelShape=IPosition(3,3,3,1);
        itsModel.reset(new Array<Float>(itsModelShape));
        itsModel->set(0.0);
        itsPrior.reset(new Array<Float>(itsModelShape));
        itsPrior->set(0.0);
        itsMask.reset(new Array<Float>(itsModelShape));
        itsMask->set(1.0);
        itsResidual.reset(new Array<Float>(itsModelShape));
        itsResidual->set(10.0);
        itsStep.reset(new Array<Float>(itsModelShape));
        itsStep->set(1.0);
	itsEntropy=EntropyBase<Float>::ShPtr(new Emptiness<Float>());
        itsEntropy->setPrior(*itsPrior);
        itsEntropy->setMask(*itsMask);
      }
      void testEntropyI() {
	itsEntropy=EntropyBase<Float>::ShPtr(new EntropyI<Float>());
        itsPrior->set(3.0);
        itsEntropy->setPrior(*itsPrior);
        itsEntropy->setMask(*itsMask);
        itsModel->set(1.0);
        CPPUNIT_ASSERT(abs(itsEntropy->entropy(*itsModel)-2.19722)<1e-4);
        Matrix<Float> GDG=itsEntropy->formGDGStep(*itsModel, *itsResidual, *itsStep);
        CPPUNIT_ASSERT(abs(GDG(1,1)-3600.0)<1);
        CPPUNIT_ASSERT(abs((*itsStep)(IPosition(3,1,1,0))-1.09861)<1e-5);
        CPPUNIT_ASSERT(abs(itsEntropy->formGDS(*itsModel, *itsResidual, *itsStep)+10.8625)<1e-3);
      }

      void testEmptiness() {
        itsModel->set(1.0);
        CPPUNIT_ASSERT(abs(itsEntropy->entropy(*itsModel)-3.90403)<1e-3);
        CPPUNIT_ASSERT(abs(itsEntropy->formGDGStep(*itsModel, *itsResidual, *itsStep)(IPosition(3,1,1,0))-8571.95)<1e-1);
        CPPUNIT_ASSERT(abs((*itsStep)(IPosition(3,1,1,0))+1.81343)<1e-4);
        CPPUNIT_ASSERT(abs(itsEntropy->formGDS(*itsModel, *itsResidual, *itsStep)+12.4299<1e-3));
      }
      void tearDown() {
        itsEntropy.reset();
	itsModel.reset();
	itsPrior.reset();
	itsMask.reset();
	itsResidual.reset();
	itsStep.reset();
      }
      
    private:

      boost::shared_ptr<EntropyBase<Float> > itsEntropy;

      IPosition itsModelShape;

      boost::shared_ptr<Array<Float> > itsModel;

      boost::shared_ptr<Array<Float> > itsPrior;

      boost::shared_ptr<Array<Float> > itsMask;

      boost::shared_ptr<Array<Float> > itsResidual;

      boost::shared_ptr<Array<Float> > itsStep;

    };
    
  } // namespace synthesis
  
} // namespace askap

