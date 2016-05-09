/// @file tmssummary.cc
///
/// @brief inspect a measurement set using various MSSummary functions
///
/// @copyright (c) 2015 CSIRO
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
/// @author Daniel Mitchell <daniel.mitchell@csiro.au>

//# Includes
#include <string>
#include <casacore/casa/aips.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MSOper/MSSummary.h>
#include <casacore/ms/MSOper/MSLister.h>
#include <askap/AskapError.h>

//DAM testing stderr output
#include <askap/AskapLogging.h>
#include <askap/Log4cxxLogSink.h>
ASKAP_LOGGER(logger, "mslist");

void usage(void) {
  std::cerr << "Usage: " << "mslist [--verbose] MS_filename" << std::endl <<
               "  --brief               brief listing (default with verbose)" << std::endl <<
               "  --full                more extensive listing" << std::endl <<
               "  --what                what what observed? (fields, times, etc.)" << std::endl <<
               "  --how                 how was it observed? (antennas, frequencies, etc.)" << std::endl <<
               "  --tables              list tables in measurement set" << std::endl <<
               "  --data                display correlation data" << std::endl <<
               "  --verbose             verbose output" <<
               std::endl;
  std::cerr << "Data selection parameters used with --data. Standard CASA options." << std::endl <<
               "  --datacolumn=str      " << std::endl <<
               "  --field=str           " << std::endl <<
               "  --spw=str             default=all. e.g. --spw='0:5~10'" << std::endl <<
               "  --antenna=str         default=all. e.g. --antenna='3&&4;1'" << std::endl <<
               "  --timerange=str       default=all. e.g. --timerange='<2014/09/14/10:48:50.0'" << std::endl <<
               "  --correlation=str     default=all. e.g. --correlation=XX" << std::endl <<
               "  --scan=str            default=all. e.g. --scan=9" << std::endl <<
               //"  --feed=str            " << std::endl <<
               //"  --array=str           " << std::endl <<
               "  --uvrange=str         " << std::endl <<
               //"  --average=str         " << std::endl <<
               //"  --showflags=bool      1 or 0. default=0" << std::endl <<
               //"  --msselect=str        " << std::endl <<
               "  --pagerows=int        default=50" << std::endl <<
               "  --listfile=str        default=stdout" << std::endl <<
               std::endl;
}

casa::String passMsParameter(int argc, const char** argv)
{
  for(int i=1; i<argc; i++) {
      std::string arg = argv[i];
      if(arg.find("-")!=0) {
          return arg;
      }
  }
  return "";
}

casa::Bool passFlagParameter(int argc, const char** argv, const char* par)
{
  for(int i=1; i<argc; i++) {
      if(!std::strcmp(argv[i], par)) return casa::True;
  }
  return casa::False;
}

template <typename T>
T passFlaggedParameter(int argc, const char** argv, const char* par, T def)
{
  T val = def;
  for(int i=1; i<argc; i++) {
      std::string arg = argv[i];
      if(arg.find(par)==0) {
          casa::Int pos = arg.find("=")+1;
          val = arg.substr(pos,arg.size()-pos);
          break;
      }
  }
  return val;
}

template<>
casa::Long passFlaggedParameter(int argc, const char** argv, const char* par, casa::Long def)
{
  casa::Long val = def;
  for(int i=1; i<argc; i++) {
      std::string arg = argv[i];
      if(arg.find(par)==0) {
          casa::Int pos = arg.find("=")+1;
          val = atoi(arg.substr(pos,arg.size()-pos).c_str());
          break;
      }
  }
  return val;
}

template<>
casa::Bool passFlaggedParameter(int argc, const char** argv, const char* par, casa::Bool def)
{
  casa::Bool val = def;
  for(int i=1; i<argc; i++) {
      std::string arg = argv[i];
      if(arg.find(par)==0) {
          casa::Int pos = arg.find("=")+1;
          casa::Int val = atoi(arg.substr(pos,arg.size()-pos).c_str());
          switch(val) {
            case 0:
              return casa::False;
              break;
            case 1:
              return casa::True;
              break;
            default:
              ASKAPTHROW(askap::AskapError, "Unknown bool (use 0 or 1): " << arg);
              break;
          }
      }
  }
  return val;
}

int main(int argc, const char* argv[])
{

  try {

    if (argc < 2) {
        usage();
        ASKAPTHROW(askap::AskapError, "mslist requires a MS file");
    }

    // pass command line parameters
    casa::Bool brief = passFlagParameter(argc, argv, "--brief");
    casa::Bool full = passFlagParameter(argc, argv, "--full");
    casa::Bool what = passFlagParameter(argc, argv, "--what");
    casa::Bool how = passFlagParameter(argc, argv, "--how");
    casa::Bool tables = passFlagParameter(argc, argv, "--tables");
    casa::Bool data = passFlagParameter(argc, argv, "--data");
    casa::Bool verbose = passFlagParameter(argc, argv, "--verbose");
    // pass parameters for MSLister
    casa::String datacolumn  = passFlaggedParameter<casa::String>(argc, argv, "--datacolumn=","");
    casa::String field       = passFlaggedParameter<casa::String>(argc, argv, "--field=","");
    casa::String spw         = passFlaggedParameter<casa::String>(argc, argv, "--spw=","");
    casa::String antenna     = passFlaggedParameter<casa::String>(argc, argv, "--antenna=","");
    casa::String timerange   = passFlaggedParameter<casa::String>(argc, argv, "--timerange=","");
    casa::String correlation = passFlaggedParameter<casa::String>(argc, argv, "--correlation=","");
    casa::String scan        = passFlaggedParameter<casa::String>(argc, argv, "--scan=","");
    //casa::String feed        = passFlaggedParameter<casa::String>(argc, argv, "--feed=","");
    //casa::String array       = passFlaggedParameter<casa::String>(argc, argv, "--array=","");
    casa::String uvrange     = passFlaggedParameter<casa::String>(argc, argv, "--uvrange=","");
    //casa::String average     = passFlaggedParameter<casa::String>(argc, argv, "--average=","");
    //casa::Bool   showflags   = passFlaggedParameter<casa::Bool>(argc, argv, "--showflags=",casa::True);
    //casa::String msselect    = passFlaggedParameter<casa::String>(argc, argv, "--msselect=","");
    casa::Long   pagerows    = passFlaggedParameter<casa::Long>(argc, argv, "--pagerows=",50);
    casa::String listfile    = passFlaggedParameter<casa::String>(argc, argv, "--listfile=","");

    casa::String msfile = passMsParameter(argc, argv);

    // these parameters do not seem to make a difference, so leave them out
    casa::String feed        = "";
    casa::String array       = "";
    casa::String average     = "";
    casa::Bool showflags     = casa::False;
    casa::String msselect    = "";

    // This now needed as the interface to list has changed
    casa::String observation = "";

    // set default behaviour to brief
    if (!brief && !full && !what && !how && !tables && !data) {
        brief = casa::True;
        verbose = casa::True;
    }

    /*
    //DAM testing stderr output
    // the following is screwing up the formatting
    // adding \n before each mss call helps, but still get the annoying messages
    std::ifstream config("askap.log_cfg", std::ifstream::in);
    if (config) {
        ASKAPLOG_INIT("askap.log_cfg");
    } else {
        std::ostringstream ss;
        ss << argv[0] << ".log_cfg";
        ASKAPLOG_INIT(ss.str().c_str());
    }
    // Ensure that CASA log messages are captured
    casa::LogSinkInterface* globalSink = new askap::Log4cxxLogSink();
    casa::LogSink::globalSink(globalSink);
    */

    //DAM testing stderr output
    casa::LogSink globalSink(casa::LogMessage::NORMAL, &std::cout);
    casa::LogIO os(globalSink);

    //casa::LogIO os;

    casa::MeasurementSet ms(msfile, casa::Table::Old);
    casa::MSSummary mss(ms);

    if (full || verbose) {
        if (verbose) mss.listTitle(os);
    }

    if (brief) {
        mss.listMain(os, verbose);
    }
    else {

        if ( full ) {
            mss.listWhere(os, verbose);
        }
        if ( full || what ) {
            mss.listWhat(os, verbose);
        }
        if ( full || how ) {
            mss.listHow(os, verbose);
        } 
        if ( full ) {
            if (!ms.source().isNull()) {
                mss.listSource(os,verbose);
            }
            if (!ms.sysCal().isNull()) {
                mss.listSysCal(os,verbose);
            }
            if (!ms.weather().isNull()) {
                mss.listWeather(os,verbose);
            }
            mss.listHistory(os);
        }
        if ( tables ) {
            mss.listTables(os,verbose);
        }

    }

    if (data) {

        casa::MSLister msl(ms, os);

        std::cout << std::endl;

        //msl.setPage(width=120, height=20);
        //msl.setFormat(ndec=2);
        //msl.setPrecision( precTime=1, precUVDist=0, Int precAmpl=3, precPhase=1, precWeight=0 );

        //use !verbose to restrict some data output? miriad: "a maximum of 6 channels is printed"

        casa::String options = "";

        msl.list(options, datacolumn, field, spw,
                 antenna, timerange, correlation,
                 scan, feed, array,observation,uvrange,
                 average, showflags, msselect,
                 pagerows, listfile);

    }

  }
  catch (const askap::AskapError& e) {
    std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
    exit(1);
  }
  catch (const std::exception& e) {
     std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what()
            << std::endl;
    exit(1);
  }

  return 0;
}

