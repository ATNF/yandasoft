///////////////////////////////////////////////////////////////////////////////
//
//  CommandLineParser - a set of classes to work with the command line
//                      parameters
//
//

#ifndef __COMMANDLINEPARSER_H
#define __COMMANDLINEPARSER_H

#include <stdexcept>
#include <string>
#include <list>
#include <sstream>

namespace cmdlineparser {

// XParser
// a special exception class to distinguish easily a wrongly set command line
// from a more serious stuff
struct XParser : public std::runtime_error {
  XParser();
  XParser(const std::string &s);
};

// unused argument in the command line
struct XParserExtra : public XParser {};

class Parser;

// BaseParameter
// a common abstract base class to benefit from the container of pointers
class BaseParameter {
    bool isDefined; // true if this parameter is defined
    friend class Parser;
public:
    BaseParameter() throw();
    virtual ~BaseParameter();
    bool defined() const throw();
    void setDefinedFlag(bool is_defined) throw();
protected:
    // should return a number of arguments required to define this parameter
    virtual size_t howManyArgs2Define() const throw() = 0;
    
    // should fill itself from a vector of parameters
    // (relevant parameters start from the 0th element)
    virtual void fill(const char **argv) = 0;

    // should test whether the specified vector of parameters
    // comforms to a particular instance of this class
    virtual bool test(const char **argv) const = 0;
};

// FlagParameter    - a flag parameter, e.g. "-C"
class FlagParameter : virtual public BaseParameter  {
   std::string flag;
public:
   FlagParameter(const std::string &iflag);
protected:
   // return a number of arguments are required to define this parameter
   virtual size_t howManyArgs2Define() const throw();

   // fill itself from a vector of parameters
   // (relevant parameters start from the 0th element)
   virtual void fill(const char **argv);

   // test whether the specified vector of parameters
   // comforms to a particular instance of this class
   virtual bool test(const char **argv) const;
};

// Generic input parameter
template<class T>
class GenericParameter : virtual public BaseParameter {
   T value;
public:

   // the argument is a default value
   GenericParameter(const T &default_value=T()) throw() :
               value(default_value) {}
   const T& getValue() const throw() { return value;}
   // more convenient access
   operator const T&() const throw() { return value;}
   void setValue(const T& in) throw() { value=in;}
protected:
   // return a number of arguments are required to define this parameter
   virtual size_t howManyArgs2Define() const throw() { return 1;}
   
   // fill itself from a vector of parameters
   // (relevant parameters start from the 0th element)
   virtual void fill(const char **argv) {
       std::istringstream is(argv[0]);
       is>>value;
       if (!is)
           throw XParser(std::string("Generic command line parameter: can't extract from ")+
	              argv[0]);
       setDefinedFlag(true);		      
   }

   // test whether the specified vector of parameters
   // comforms to a particular instance of this class
   virtual bool test(const char **) const
   {
     return true;
   }
};

// Input parameter, which follows a flag (e.g. -f file_name)
template<class T>
struct FlaggedParameter : public GenericParameter<T>,
                                public FlagParameter
{
  FlaggedParameter(const std::string &flag,
                          const T &default_value=T()):
	  GenericParameter<T>(default_value), FlagParameter(flag) {}
protected:
   // return a number of arguments are required to define this parameter
   virtual size_t howManyArgs2Define() const throw() { return 2;}

   // fill itself from a vector of parameters
   // (relevant parameters start from the 0th element)
   virtual void fill(const char **argv) {
       GenericParameter<T>::fill(argv+1);
   }
   // test whether the specified vector of parameters
   // comforms to a particular instance of this class
   virtual bool test(const char **argv) const {
      return FlagParameter::test(argv);  
   }
};

// Main class which does the job of filling parameters
class Parser {
std::list<BaseParameter *>  params; // this class doesn't have the ownership
                              // of parameters, therefore no destructor
std::list<int> absence_actions;  // what to do if a parameter is missing
public:			      
   // actions taken if a parameter is not defined
   static const int throw_exception = 0;
   static const int return_default = 1; // i.e. the parameter is optional

   // a set of handy constructors to add a few parameters straight away
   // in all cases all parameters are considered mandatory
   Parser(); // an empty object
   Parser(BaseParameter &par);
   Parser(BaseParameter &par1, BaseParameter &par2);
   Parser(BaseParameter &par1, BaseParameter &par2, BaseParameter &par3);

   // add one more parameter
   void add(BaseParameter &par,int absence_action = throw_exception);
	
   // process a command line, after this all mandatory parameters
   // have to be filled
   void process(int argc, const char **argv) const;

   // overloaded version of the process with different constness
   inline void process(int argc, char **argv) const
    { process(argc, const_cast<const char **>(argv));}

protected:
   bool testParameter(BaseParameter *par) const;
private:
   mutable int curarg; // current argument
   mutable int argc_buf;   //buffers for command line parameters
   mutable const char **argv_buf;
};

}
#endif // #ifndef __COMMANDLINEPARSER_H

//
///////////////////////////////////////////////////////////////////////////////
