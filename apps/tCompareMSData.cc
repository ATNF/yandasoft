/// @file
/// @brief test if two MSs have the same data
///
/// @copyright (c) 2023 CSIRO
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
#include <askap/askap_synthesis.h>

#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/Application.h>
#include <askap/askap/StatReporter.h>
#include <askap/askapparallel/AskapParallel.h>
#include <askap/dataaccess/ParsetInterface.h>

#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/IDataConverterImpl.h>
ASKAP_LOGGER(logger, ".tcomparemsdata");

// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/casa/Arrays/ArrayMath.h>

// std
#include <stdexcept>
#include <iostream>

using namespace askap;
using namespace askap::accessors;
using namespace casacore;

bool doCompare(LOFAR::ParameterSet parset, const std::string &name1, const std::string &name2, bool detail = false, float tol = 2.0e-6) {
    TableDataSource ds1(name1,TableDataSource::MEMORY_BUFFERS);
    IDataSelectorPtr sel1=ds1.createSelector();
    sel1 << parset;
    IDataConverterPtr conv1=ds1.createConverter();
    conv1->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"MHz");
    IDataSharedIter it1=ds1.createIterator(sel1,conv1);
    TableDataSource ds2(name2,TableDataSource::MEMORY_BUFFERS);
    IDataSelectorPtr sel2=ds2.createSelector();
    sel2 << parset;
    IDataConverterPtr conv2=ds2.createConverter();
    conv2->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"MHz");
    IDataSharedIter it2=ds2.createIterator(sel2,conv2);
    size_t cnt=0;
    bool match = true;
    bool nearMatch = true;
    ASKAPLOG_INFO_STR(logger, "Checking for exact and near match (tolerance = "<< tol <<")");

    for (it1.init(),it2.init();it1.hasMore()&&it2.hasMore();it1.next(),it2.next(),++cnt) {
        const bool equal = allEQ(it1->visibility(),it2->visibility());
        if (!equal) {
            const bool almostEqual = allNearAbs(it1->visibility(),it2->visibility(),tol);
            nearMatch &= almostEqual;
            if (!almostEqual) {
                const Cube<Complex> diff = it1->visibility() - it2->visibility();
                ASKAPLOG_WARN_STR(logger,"Difference in visibilities at iteration "<< cnt <<
                " : "<< real(max(abs(diff))) << " real : "<<max(real(diff))<< " imag : "<<max(imag(diff)));
                if (detail) {
                    const Cube<Complex> vis1 = it1->visibility();
                    const Cube<Complex> vis2 = it2->visibility();
                    const IPosition shape = vis1.shape();
                    ASKAPLOG_INFO_STR(logger, "data shape = "<<shape);
                    for (int i = 0; i < shape(0); i++) {
                        for (int j = 0; j < shape(1); j++) {
                            for (int k = 0; k  < shape(2); k++) {
                                if (vis1(i,j,k)!=vis2(i,j,k)) {
                                    ASKAPLOG_INFO_STR(logger," vis1("<<i<<","<<j<<","<<k<<") = "<<vis1(i,j,k) << " vis2 = "<<vis2(i,j,k)
                                        <<" diff = " << diff(i,j,k));
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }
        match &= equal;
    }
    ASKAPLOG_INFO_STR(logger, "Completed "<<cnt<<" iterations");
    if (it1.hasMore() || it2.hasMore()) {
        ASKAPLOG_WARN_STR(logger, "Inputs do not have the same number of integrations");
        return false;
    }
    if (match) {
        ASKAPLOG_INFO_STR(logger,"The data columns in the two MSs match");
    } else if (nearMatch){
        ASKAPLOG_INFO_STR(logger,"The data columns in the two MSs do not match exactly but do match to within "<< tol);
    } else {
        ASKAPLOG_WARN_STR(logger,"The data columns in the two MSs do not match");
        return false;
    }
    return true;
}

bool doCompareFlags(LOFAR::ParameterSet parset, const std::string &name1, const std::string &name2, bool detail = false) {
    TableDataSource ds1(name1,TableDataSource::MEMORY_BUFFERS);
    IDataSelectorPtr sel1=ds1.createSelector();
    sel1 << parset;
    IDataConverterPtr conv1=ds1.createConverter();
    conv1->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"MHz");
    IDataSharedIter it1=ds1.createIterator(sel1,conv1);
    TableDataSource ds2(name2,TableDataSource::MEMORY_BUFFERS);
    IDataSelectorPtr sel2=ds2.createSelector();
    sel2 << parset;
    IDataConverterPtr conv2=ds2.createConverter();
    conv2->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"MHz");
    IDataSharedIter it2=ds2.createIterator(sel2,conv2);
    size_t cnt=0;
    bool match = true;
    ASKAPLOG_INFO_STR(logger, "Checking for exact match of flags");

    for (it1.init(),it2.init();it1.hasMore()&&it2.hasMore()&&match;it1.next(),it2.next(),++cnt) {
        match &= allEQ(it1->flag(),it2->flag());
    }
    ASKAPLOG_INFO_STR(logger, "Completed "<<cnt<<" iterations");
    if (match && (it1.hasMore() || it2.hasMore())) {
        ASKAPLOG_WARN_STR(logger, "Inputs do not have the same number of integrations");
        return false;
    }
    if (match) {
        ASKAPLOG_INFO_STR(logger,"The flag columns in the two MSs match");
    } else {
        ASKAPLOG_WARN_STR(logger,"The flag columns in the two MSs do not match");
    }
    return match;
}

/// @brief application class
class CompareMSDataApp : public askap::Application
{
    public:
        int run(int argc, char* argv[]) final
        {
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;


                const std::vector<string> ms = config().getStringVector("dataset");
                if (ms.size()<2) {
                   ASKAPLOG_FATAL_STR(logger, "Need 2 entries in dataset parameter, the names of MSs to compare");
                   return 1;
                }
                const bool detail = config().getBool("detail",false);
                const float tolerance = config().getFloat("tolerance",2.0e-6);
                ASKAPLOG_INFO_STR(logger,"Comparing DATA column of "<<ms[0]<<" and "<<ms[1]);
                bool match = doCompare(config(), ms[0],ms[1],detail,tolerance);
                if (config().getBool("compareflag",false)) {
                    match &= doCompareFlags(config(), ms[0],ms[1],false);
                }
                stats.logSummary();
                return (match ? 0 : 1);

            } catch (const askap::AskapError& x) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
                std::cerr << "Askap error in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            } catch (const std::exception& x) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            }

            return 0;
        }

    private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

int main(int argc, char *argv[])
{
    CompareMSDataApp app;
    app.addParameter("profile", "p", "Write profiling output files", false);
    return app.main(argc, argv);
}
