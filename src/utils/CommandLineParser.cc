///////////////////////////////////////////////////////////////////////////////
//
//  CommandLineParser - a set of classes to work with the command line
//                      parameters
//
//

#include "CommandLineParser.h"

using namespace cmdlineparser;
using namespace std;

// XParser
// a special exception class to distinguish easily a wrongly set command line
// from a more serious stuff
XParser::XParser() :
    runtime_error("Invalid string of command line arguments") {};

XParser::XParser(const string &s) :
    runtime_error(s) {};

// BaseParameter
// a common abstract base class to benefit from the container of pointers
BaseParameter::~BaseParameter() {}

BaseParameter::BaseParameter() throw() :
    isDefined(false) {}

bool BaseParameter::defined() const throw()
{
  return isDefined;
}

void BaseParameter::setDefinedFlag(bool is_defined) throw()
{
  isDefined=is_defined;
}

// FlagParameter  - a flag parameter, e.g. "-C"
FlagParameter::FlagParameter(const std::string &iflag)
    throw(std::exception) : flag(iflag) {}

size_t FlagParameter::howManyArgs2Define() const throw()
{
  return 1;
}

// fill itself from a vector of parameters
// (relevant parameters start from the 0th element)
void FlagParameter::fill(const char **argv) throw(std::exception)
{
  if (flag!=argv[0])
      throw runtime_error("Assert in FlagParamter::fill, that shouldn't happend");
  setDefinedFlag(true);
}

//  test whether the specified vector of parameters
// comforms to a particular instance of this class
bool FlagParameter::test(const char **argv) const throw(std::exception)
{
  return flag==argv[0];
}

// Parser - main class which does the job of filling parameters

// a set of handy constructors to add a few parameters straight away
// in all cases all parameters are considered mandatory
Parser::Parser() {}

Parser::Parser(BaseParameter &par) throw(std::exception)
{
  add(par);
}

Parser::Parser(BaseParameter &par1, BaseParameter &par2)
     throw(std::exception)
{
  add(par1);
  add(par2);
}

Parser::Parser(BaseParameter &par1, BaseParameter &par2, BaseParameter &par3)
       throw(std::exception)
{
  add(par1);
  add(par2);
  add(par3);
}

// add one more parameter
void Parser::add(BaseParameter &par,int absence_action)
     throw(std::exception)
{
  params.push_back(&par);
  absence_actions.push_back(absence_action);
}

// process a command line, after this all mandatory parameters
// have to be filled
void Parser::process(int argc, const char **argv) const throw(std::exception)
{
  curarg=1;
  argc_buf=argc;
  argv_buf=argv;
  list<BaseParameter*>::const_iterator it=params.begin();
  list<int>::const_iterator ci=absence_actions.begin();
  while (it!=params.end()) {
       if (testParameter(*it)) {
      	  (*it)->fill(argv_buf+curarg);
	  curarg+=(*it)->howManyArgs2Define();	  
       } else {
          (*it)->setDefinedFlag(false);
	  if (*ci==throw_exception)
	     throw XParser();
       }
       ++ci; ++it; // next parameter
  }
  if (curarg!=argc)
      throw XParserExtra(); // something extra is present
}

bool Parser::testParameter(BaseParameter *par) const throw(std::exception)
{
  if (curarg+int(par->howManyArgs2Define())>argc_buf)
      return false; // this parameter requires too many arguments
  return par->test(argv_buf+curarg); // give only subset of parameters  
}

//
///////////////////////////////////////////////////////////////////////////////
