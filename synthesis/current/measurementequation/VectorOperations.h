/// @file
///
/// @brief Helper methods to assist operations with Vectors
/// @details It is often necessary to copy vectors of different types
/// or do other operations like subtracting them. This file contains
/// templates, which make this operation more clearly visible in the code 
/// and assist structuring of the program.
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
/// @author Max Voronkov <maxim.voronkov@csiro.au
///

#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H

#include <scimath/Mathematics/AutoDiff.h>
#include <fitting/ComplexDiff.h>


namespace askap {

namespace synthesis {

/// this namespace contains a set of helper classes to manipulate containers
/// of an arbitrary type. The main functionality, which is intended for
/// external use is still in the askap::synthesis namespace 
namespace vector_operations {

/// @brief A helper template to extract type of the template parameter
/// @details To be able to provide template specialization for various
/// types stored in the container, one needs to know the actual type,
/// which is a parameter of (templated) container type. STL containers
/// provide a typedef (called value_type) for this purpose, but CASA
/// containers lack the same template (but have it for iterators). A 
/// general solution for the problem is to use partial template 
/// specialization (metaprogramming technique). The default version is
/// intended for STL-like containers, which just uses ::value_type.
/// @ingroup measurementequation
template<typename T>
struct ValueTypeExtractor {
  /// type of the value of the container is given by value_type by default
  typedef typename T::value_type type;
  // there is deliberately no container_type defined in this template
  // it prevent wrong code implying by mistake a reference semantics for
  // STL containers from compilation. 
};


/// Specialization of ValueTypeExtractor for casa::Vector
/// @ingroup measurementequation
template<typename Y>
struct ValueTypeExtractor<casa::Vector<Y> > {
 /// type of the value of the container is a template argument for casa::Vector
 typedef Y type;
 /// type of the container. Defining this allowing reference semantics
 typedef casa::Vector<Y> container_type;
};

/// Specialization of ValueTypeExtractor for casa::Array
/// @ingroup measurementequation
template<typename Y>
struct ValueTypeExtractor<casa::Array<Y> > {
  /// type of the value of the container is a template argument for casa::Array
  typedef Y type;
  /// type of the container. Defining this allowing reference semantics
  typedef casa::Array<Y> container_type;
};

/// @brief Simple increment of the iterator
/// @details Flattening of a container on-demand means that one element of
/// one vector may correspond to many elements of the other. Therefore, there is
/// no, in general, a one-to-one correspondence between the increment operations 
/// for the input and output iterators. This class does a simple increment of 
/// the given iterator, but can be overridden for other behavior
/// @ingroup measurementequation
struct BasicIncrementor {
  
  /// @brief go to next element
  /// @details Proceed to the next element, iterator is incremented if necessary.
  /// @param[in] it iterator to work with 
  template<typename Iter>
  inline static void increment(Iter &it) {++it;}
};

/// @brief Increment of the iterator appropriate for complex values
/// @details Complex numbers are represented as two doubles in the flattened
/// vector. This class increment iterator once for every second call as well
/// as tracks whether the current value should correspond to real or
/// imaginary part. 
/// @note There is no need to derive this class from 
/// BaseIncrementor, because we're not using run-time polymorphism here.
/// Templates ensure the most efficient operation where the correct function is
/// chosen at the compile time.
/// @ingroup measurementequation
struct ComplexNumberIncrementor {

  /// @brief construct the object
  /// @details The constructor is required to fix the order of two double
  /// numbers representing a Complex number. The first value is always
  /// the real part, which is followed by the imaginary part.
  ComplexNumberIncrementor() : itsRealNow(true) {}
   
   /// @brief increment operation
   /// @details This method increment iterator twice less frequently
   /// and changes between real and imaginary part on every call.
   /// @param[in] it iterator to work with 
   template<typename Iter>
   inline void increment(Iter &it) const {
     if (itsRealNow) {
         itsRealNow=false;
     } else {
         itsRealNow=true;
         ++it;
     }
   }
protected:
  /// @brief is it a real or imaginary part
  /// @return true if the current output is a real part of the complex value
  inline bool isRealNow() const throw() { return itsRealNow;}       
private:
   /// true if currently this accessor returns a real part
   mutable bool itsRealNow;

};


/// @brief Helper class to wrap around input iterator
/// @details The purpose of this class is to extract the required value 
/// regardless on the actual type returned by the iterator. The default 
/// implementation returns the value itself.
/// However, it can return, for example casa::AutoDiff<ValueType>, which has
/// to be converted to ValueType by calling the value() method. This is 
/// taken care of by the partial specialization.
/// @ingroup measurementequation
template<typename ValueType>
struct InputValueAccessor : public BasicIncrementor {
  /// @brief access to the value
  /// @details Default action is to return the value itself
  /// @param[in] in the value to work with
  inline ValueType operator()(const ValueType &in) const throw() {return in;}
};

/// @brief InputValueAccessor specialization for casa::AutoDiff
/// @ingroup measurementequation
template<typename T>
struct InputValueAccessor<casa::AutoDiff<T> > : public BasicIncrementor  {
  /// @brief access to the value
  /// @details This specialization just calls value() method of its input AutoDiff
  /// @param[in] in the AutoDiff to work with
  inline T operator()(const casa::AutoDiff<T> &in) const {return in.value();}
};

/// @brief InputValueAccessor specialization for casa::Complex
/// @ingroup measurementequation
template<>
struct InputValueAccessor<casa::Complex> : public ComplexNumberIncrementor {
   
   /// @brief access to the value
   /// @details This access method of this specialization returns real or
   /// imaginary part of the input complex number, depending on whether it is
   /// an even or odd call to the increment operator.
   /// @param[in] in input complex number
   inline casa::Double operator()(const casa::Complex &in) const 
      { return isRealNow()? real(in) : imag(in); }
};

/// @brief InputValueAccessor specialization for ComplexDiff
/// @ingroup measurementequation
template<>
struct InputValueAccessor<scimath::ComplexDiff> : public ComplexNumberIncrementor {
   
   /// @brief access to the value
   /// @details This access method of this specialization returns real or
   /// imaginary part of the value of the input ComplexDiff object, depending on 
   /// whether it is an even or odd call to the increment operator.
   /// @param[in] in input complex number
   inline casa::Double operator()(const scimath::ComplexDiff &in) const 
      { return isRealNow()? real(in.value()) : imag(in.value()); }
};


/// @brief Helper class to wrap around output iterator
/// @details The purpose of this class is to store the required value 
/// regardless on the actual type expected by the iterator. The default 
/// implementation assumes that no type change is required. The idea is
/// very similar to that for InputValueAccessor
/// @ingroup measurementequation
template<typename ValueType>
struct OutputValueAccessor : public BasicIncrementor {
  /// @brief write the value
  /// @details Default action is to copy the value without type conversion
  /// @param[in] in the value to write
  /// @param[out] out destination
  inline void write(const ValueType &in, ValueType& out) const throw() 
         {out=in;}
  
  /// @brief subtract the value
  /// @details To form residuals we need to perform '-=' operation.
  /// Default action is to subtract the value without any type conversion
  /// @param[in] in subtrahend
  /// @param[in,out] out minuend for input and result for output
  inline void subtract(const ValueType &in, ValueType& out) const throw() 
         {out-=in;} 
         
  /// @brief add the value
  /// @details To do the prediction we need to perform '+=' operation.
  /// Default action is to add the value without any type conversion
  /// @param[in] in a value to add
  /// @param[in,out] out another item and result for output
  inline void add(const ValueType &in, ValueType& out) const throw() 
         {out+=in;} 

  /// @brief add the value with scaling
  /// @details To do the prediction we need to perform '+=' operation.
  /// Default action is to add the value without any type conversion
  /// @param[in] in a value to add
  /// @param[in,out] out another item and result for output
  /// @param[in] factor scaling factor
  inline void addScaled(const ValueType &in, ValueType& out, const ValueType &factor) const throw() 
         {out += in * factor;} 
               
};

/// @brief OutputValueAccessor specialization for casa::Complex
/// @ingroup measurementequation
template<>
struct OutputValueAccessor<casa::Complex> : public ComplexNumberIncrementor {
   /// @brief default constructor
   /// @details strictly speaking, this method is not required for correct
   /// operations. However, the compiler gives a
   /// warning that itsRealPart data member can be uninitialized, if this
   /// constructor is not present. 
   OutputValueAccessor() : itsRealPart(0.) {}

   /// @brief write the value
   /// @details This method of this particular specialization stores real or
   /// imaginary part of the input complex number, depending on whether it is
   /// even or odd call to the increment operator.
   /// @param[in] in the value to write
   /// @param[out] out destination
   inline void write(casa::Double in, casa::Complex& out) const throw() 
      {
        if (isRealNow()) {
            itsRealPart=in;
        } else {
            out=casa::Complex(itsRealPart,in);
        }
      }
      
  /// @brief subtract the value
  /// @details To form residuals we need to perform '-=' operation.
  /// This method of this particular specialization stores real or
  /// imaginary part of the input complex number, depending on whether it is
  /// even or odd call to the increment operator.
  /// @param[in] in subtrahend
  /// @param[in,out] out minuend for input and result for output
  inline void subtract(casa::Double in, casa::Complex& out) const throw() 
      {
        if (isRealNow()) {
            itsRealPart=in;
        } else {
            out-=casa::Complex(itsRealPart,in);
        }
      }      

  /// @brief add the value
  /// @details To do the prediction we need to perform '+=' operation.
  /// This method of this particular specialization stores real or
  /// imaginary part of the input complex number, depending on whether it is
  /// even or odd call to the increment operator.
  /// @param[in] in a value to add
  /// @param[in,out] out another item and result for output
  inline void add(casa::Double in, casa::Complex& out) const throw() 
      {
        if (isRealNow()) {
            itsRealPart=in;
        } else {
            out+=casa::Complex(itsRealPart,in);
        }
      }      

  /// @brief add the value with scaling
  /// @details To do the prediction we need to perform '+=' operation.
  /// Default action is to add the value without any type conversion
  /// @param[in] in a value to add
  /// @param[in,out] out another item and result for output
  /// @param[in] factor scaling factor
  inline void addScaled(casa::Double in, casa::Complex& out, const casa::Complex& factor) const throw() 
  {
    if (isRealNow()) {
        itsRealPart=in;
    } else {
        out += casa::Complex(itsRealPart,in) * factor;
    }
  }      
       
private:
   /// a buffer for the real part of the value (write happens when
   /// the imaginary part is known
   mutable casa::Double itsRealPart;
};

/// @brief OutputValueAccessor specialization for ComplexDiff
/// @details this is an empty class to prevent compilation with the output to
/// a container of ComplexDiffs (because splitting into real and imaginary values
/// would loose derivative information) 
/// @ingroup measurementequation
template<>
struct OutputValueAccessor<scimath::ComplexDiff>  {
};


} // namespace vector_operations

/// @brief copy 1D-vector and flatten it on demand
/// @details A number of situation requires copying 1D vectors of different
/// types (e.g. composing vector<double> from casa::Vector<casa::Complex> or
/// casa::Vector<casa::AutoDiff<double> >). This templated function is intended 
/// to make this copy operation more clearly visibile in the code. The fitting
/// process deals with complex parameters as with two double-valued parameters.
/// This method does flattening if necessary, i.e. if one of the arguments has
/// complex type, it will be copied to/from two adjacent elements of a 
/// double-typed vector. All functionality is achieved via C++ templates.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @note Input and output vectors can, in principle, be of any type which
/// have begin() and end() method delivering iterators, which behave in the STL
/// sense (only operator!=, operator++ and operator* of these iterators are used).
/// However, ValueTypeExtractor should be specialized for this type, unless the
/// default implementation is adequate. The default implementation of 
/// ValueTypeExtractor is valid for all STL containers.
/// @ingroup measurementequation
template<typename InType,typename OutType>
inline void copyVector(const InType& inVec, OutType &outVec) 
{
  typename InType::const_iterator ci = inVec.begin();
  typename OutType::iterator it = outVec.begin();
  vector_operations::InputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<InType>::type> iva;
  vector_operations::OutputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<OutType>::type> ova;
  for (; ci!=inVec.end() && it != outVec.end(); 
       iva.increment(ci),ova.increment(it)) {
       ova.write(iva(*ci),*it);
  }
}

/// @brief copy 1D-vector and flatten it on demand
/// @details This version is intended for slices of casa containers, which
/// use reference semantics and can be temporary objects.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType, typename OutType>
inline void copyVector(const InType& inVec, const OutType &outVec)
{
  // use this to ensure that no STL containers are called with this
  // syntax. Passing const reference is necessary for casa containers 
  // only, which use reference semantics. If we allow the following code
  // for STL containers, it would modify the copy and the changes would not
  // propagate to the destination container.
  typename vector_operations::ValueTypeExtractor<OutType>::container_type 
                                                        OutVecRef(outVec);
  copyVector(inVec,OutVecRef); 
}

/// @brief copy 1D-vector of derivatives
/// @details This is a special case of copyVector, which extracts a
/// derivative instead of the value from AutoDiff. We could have parameterized
/// copyVector with one more parameter, but it would not make the syntax cleaner.
/// @param[in] par a number of derivative of interest
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType,typename OutType>
inline void copyDerivativeVector(size_t par, const InType& inVec, OutType& outVec) 
{
  typename InType::const_iterator ci = inVec.begin();
  typename OutType::iterator it = outVec.begin();
  vector_operations::InputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<InType>::type> iva;
  vector_operations::OutputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<OutType>::type> ova;
  for (; ci!=inVec.end() && it != outVec.end(); 
       iva.increment(ci),ova.increment(it)) {
       ova.write(ci->derivative(par),*it);
  }
}

/// @brief copy 1D-vector of derivatives
/// @details This version is intended for slices of casa containers, which
/// use reference semantics and can be temporary objects.
/// @param[in] par a number of derivative of interest
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType, typename OutType>
inline void copyDerivativeVector(size_t par, const InType& inVec, 
                                 const OutType &outVec)
{
  // use this to ensure that no STL containers are called with this
  // syntax. Passing const reference is necessary for casa containers 
  // only, which use reference semantics. If we allow the following code
  // for STL containers, it would modify the copy and the changes would not
  // propagate to the destination container.
  typename vector_operations::ValueTypeExtractor<OutType>::container_type 
                                                        OutVecRef(outVec);
  copyDerivativeVector(par, inVec, OutVecRef); 
}

/// @brief copy 1D-vector of derivatives by real parts
/// @details This is a special case of copyVector, which extracts a
/// derivative by real part of paramters instead of the value from ComplexDiff. 
/// @param[in] par a name of the parameter/derivative of interest
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType,typename OutType>
inline void copyReDerivativeVector(const std::string &par, const InType& inVec, 
                                   OutType& outVec) 
{
  typename InType::const_iterator ci = inVec.begin();
  typename OutType::iterator it = outVec.begin();
  vector_operations::InputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<InType>::type> iva;
  vector_operations::OutputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<OutType>::type> ova;
  for (; ci!=inVec.end() && it != outVec.end(); 
       iva.increment(ci),ova.increment(it)) {
       ova.write(iva((*ci).derivRe(par)),*it);
  }
}

/// @brief copy 1D-vector of derivatives by real parts
/// @details This version is intended for slices of casa containers, which
/// use reference semantics and can be temporary objects.
/// @param[in] par a name of the parameter/derivative of interest
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType, typename OutType>
inline void copyReDerivativeVector(const std::string &par, const InType& inVec, 
                                 const OutType &outVec)
{
  // use this to ensure that no STL containers are called with this
  // syntax. Passing const reference is necessary for casa containers 
  // only, which use reference semantics. If we allow the following code
  // for STL containers, it would modify the copy and the changes would not
  // propagate to the destination container.
  typename vector_operations::ValueTypeExtractor<OutType>::container_type 
                                                        OutVecRef(outVec);
  copyReDerivativeVector(par, inVec, OutVecRef); 
}


/// @brief copy 1D-vector of derivatives by imaginary parts
/// @details This is a special case of copyVector, which extracts a
/// derivative by imaginary part of paramters instead of the value 
/// from ComplexDiff. 
/// @param[in] par a name of the parameter/derivative of interest
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType,typename OutType>
inline void copyImDerivativeVector(const std::string &par, const InType& inVec, 
                                   OutType& outVec) 
{
  typename InType::const_iterator ci = inVec.begin();
  typename OutType::iterator it = outVec.begin();
  vector_operations::InputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<InType>::type> iva;
  vector_operations::OutputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<OutType>::type> ova;
  for (; ci!=inVec.end() && it != outVec.end(); 
       iva.increment(ci),ova.increment(it)) {
       ova.write(iva((*ci).derivIm(par)),*it);
  }
}

/// @brief copy 1D-vector of derivatives by imaginary parts
/// @details This version is intended for slices of casa containers, which
/// use reference semantics and can be temporary objects.
/// @param[in] par a name of the parameter/derivative of interest
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType, typename OutType>
inline void copyImDerivativeVector(const std::string &par, const InType& inVec, 
                                 const OutType &outVec)
{
  // use this to ensure that no STL containers are called with this
  // syntax. Passing const reference is necessary for casa containers 
  // only, which use reference semantics. If we allow the following code
  // for STL containers, it would modify the copy and the changes would not
  // propagate to the destination container.
  typename vector_operations::ValueTypeExtractor<OutType>::container_type 
                                                        OutVecRef(outVec);
  copyImDerivativeVector(par, inVec, OutVecRef); 
}
 

/// @brief subtract one 1D-vector from another, flatten on demand
/// @details See copyVector for more information. This method subtracts
/// the input vector from the output one instead of copying to it.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
template<typename InType,typename OutType>
inline void subtractVector(const InType& inVec, OutType& outVec) 
{
  typename InType::const_iterator ci = inVec.begin();
  typename OutType::iterator it = outVec.begin();
  vector_operations::InputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<InType>::type> iva;
  vector_operations::OutputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<OutType>::type> ova;
  for (; ci!=inVec.end() && it != outVec.end(); 
       iva.increment(ci),ova.increment(it)) {
       ova.subtract(iva(*ci),*it);
  }
}

/// @brief subtract one 1D-vector from another, flatten on demand
/// @details This version is intended for slices of casa containers, which
/// use reference semantics and can be temporary objects.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType, typename OutType>
inline void subtractVector(const InType& inVec, const OutType &outVec)
{
  // use this to ensure that no STL containers are called with this
  // syntax. Passing const reference is necessary for casa containers 
  // only, which use reference semantics. If we allow the following code
  // for STL containers, it would modify the copy and the changes would not
  // propagate to the destination container.
  typename vector_operations::ValueTypeExtractor<OutType>::container_type 
                                                        OutVecRef(outVec);
  subtractVector(inVec,OutVecRef); 
}


/// @brief add one 1D-vector to another, flatten on demand
/// @details See copyVector for more information. This method adds
/// the input vector to the output one instead of copying to it.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
template<typename InType,typename OutType>
inline void addVector(const InType& inVec, OutType& outVec) 
{
  typename InType::const_iterator ci = inVec.begin();
  typename OutType::iterator it = outVec.begin();
  vector_operations::InputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<InType>::type> iva;
  vector_operations::OutputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<OutType>::type> ova;
  for (; ci!=inVec.end() && it != outVec.end(); 
       iva.increment(ci),ova.increment(it)) {
       ova.add(iva(*ci),*it);
  }
}

/// @brief add one 1D-vector to another, flatten on demand
/// @details This version is intended for slices of casa containers, which
/// use reference semantics and can be temporary objects.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType, typename OutType>
inline void addVector(const InType& inVec, const OutType &outVec)
{
  // use this to ensure that no STL containers are called with this
  // syntax. Passing const reference is necessary for casa containers 
  // only, which use reference semantics. If we allow the following code
  // for STL containers, it would modify the copy and the changes would not
  // propagate to the destination container.
  typename vector_operations::ValueTypeExtractor<OutType>::container_type 
                                                        OutVecRef(outVec);
  addVector(inVec,OutVecRef); 
}

/// @brief add one 1D-vector to another after scaling, flatten on demand
/// @details See addVector for more information. This method scales values of the
/// input vector before addition.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @param[in] factor scaling factor (type of the output vector values)
template<typename InType,typename OutType>
inline void addScaledVector(const InType& inVec, OutType& outVec, 
            const typename vector_operations::ValueTypeExtractor<OutType>::type &factor) 
{
  typename InType::const_iterator ci = inVec.begin();
  typename OutType::iterator it = outVec.begin();
  vector_operations::InputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<InType>::type> iva;
  vector_operations::OutputValueAccessor<typename 
                vector_operations::ValueTypeExtractor<OutType>::type> ova;
  for (; ci!=inVec.end() && it != outVec.end(); 
       iva.increment(ci),ova.increment(it)) {
       ova.addScaled(iva(*ci), *it, factor);
  }
}

/// @brief add one 1D-vector to another after scaling, flatten on demand
/// @details This version is intended for slices of casa containers, which
/// use reference semantics and can be temporary objects.
/// @param[in] inVec input vector
/// @param[in] outVec output vector
/// @ingroup measurementequation
template<typename InType, typename OutType>
inline void addScaledVector(const InType& inVec, const OutType &outVec, 
            const typename vector_operations::ValueTypeExtractor<OutType>::type &factor)
{
  // use this to ensure that no STL containers are called with this
  // syntax. Passing const reference is necessary for casa containers 
  // only, which use reference semantics. If we allow the following code
  // for STL containers, it would modify the copy and the changes would not
  // propagate to the destination container.
  typename vector_operations::ValueTypeExtractor<OutType>::container_type 
                                                        OutVecRef(outVec);
  addScaledVector(inVec,OutVecRef,factor); 
}


} // namespace synthesis

} // namespace askap

#endif //#ifndef VECTOR_OPERATIONS_H
