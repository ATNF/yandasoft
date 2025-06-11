/// @file
/// @brief Class encapsulating uv-weights related interactions with the Params class
/// @details We ultimately store and ship around uv weights in the Params class (same as, for example, the  model). 
/// It is handy to gather all required operations in one place, so we're always consistent with parameter naming, etc.
/// Index translation is usually not possible to derive from parameters themselves and, besides, we can have different
/// mapping for research purposes (e.g. applying the weight determined for one beam to the other, etc) or use more
/// complex strategy in the future. This complexity should probably be handled through this class too.
/// @note The uv weights are always stored as single precision floats (this is sufficient for weights as they're in
/// the uv domain), but arrays in Params can be both double and single precision floats depending on the compile flags
/// (look for imtype typedef). As a result, having double precision option used would result in two copies of the weight
/// grid to exist (one in Params as double and one in the applicator as float) as opposed to reference semantics. If we
/// want to fix this, we'd probably have to implement run time coexistance of float and double parameters (as opposed to
/// compile time switch via imtype).
///
/// the parameters corresponding to uv-weights have the following style:
///     uvweight.name.index where name corresponds to the image parameter name the particular weight applies to and
///                         index is an unsigned integer index corresponding to the index used in the collection.
///     a special uvweight_indices.name parameter describes index translation coefficients (expected to be a 3-element vector). 
///     If absent, the trivial translation is assumed and only weights with index=0 should be present. Note, this is the area
///     where future changes are likely and low-level details are hidden from the end user on purpose. The index translator is
///     returned or setup via the interface (IUVWeightIndexTranslator) and this class is responsible for describing it right in
///     the Params class.
/// 
///   allow wildcard in name?
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

// own includes
#include <askap/gridding/UVWeightParamsHelper.h>
#include <askap/gridding/GenericUVWeightIndexTranslator.h> 
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapUtil.h>

ASKAP_LOGGER(logger, ".gridding.uvweightparamshelper");

namespace askap {

namespace synthesis {

/// @brief initialise for the particular params
/// @details Making it explicit to highlight intentions to use the helper for the given Params object (as the compiler
/// wouldn't allow an implicit cast), although otherwise, it wouldn't be required for such an adapter/helper.
/// @note This class is cheap to construct and destroy (essentially it is an adapter), can be created multiple times if needed.
/// An exception is thrown if the shared pointer is empty.
/// @param[in] params shared pointer to the params class to work with
UVWeightParamsHelper::UVWeightParamsHelper(const boost::shared_ptr<scimath::Params> &params) : itsParams(params) 
{
   ASKAPCHECK(params, "Attempt to initialise UVWeightParamsHelper with an empty Params shared pointer");
}

/// @brief initialise for the particular params specified via reference
/// @details This version of the constructor assumes that the ownership of the reference/pointer is managed by the caller
/// (i.e. a temporary shared pointer is created under assumption that the supplied reference would never go out of scope - 
/// this is handy for operations within one code block/method). As before, it has been made explicit to make the intentions 
/// to wrap Params class more clear in the code.
/// @param[in] params non-const reference to the params class to work with
/// @note Although it may be handy to have both const and non-const versions, it is only possible if one splits this class into
/// two (const and non-const). At present, we have to always supply a non-const reference as this class has some non-const operations.
/// One could use const_cast or mutable flag to counteract this. Although this is a bit of the technical debt, it shouldn't lead to 
/// big consequences.
UVWeightParamsHelper::UVWeightParamsHelper(scimath::Params &params) : itsParams(&params, utility::NullDeleter())
{
}



/// @brief check that uv-weight exists for a particular name
/// @param[in] name corresponding image name to query
/// @return true if uv-weight exists
bool UVWeightParamsHelper::exists(const std::string &name) const
{
   ASKAPDEBUGASSERT(itsParams);
   // we cannot query parameter name directly because the trailing flat index is unknown, 
   // need to check all parameters with the right prefix. The second parameter passed to completions() ensures
   // that fixed parameters are also included (uv-weight related parameters are expected to be fixed as we're not solving for them)
   const std::vector<std::string> completions = itsParams->completions("uvweight."+name+".", true);
   return completions.size() > 0u;
}

/// @brief obtain index translator
/// @details This method returns shared pointer to the index translation class which can be used together with 
/// the weight collection to setup weight accessor.
/// @note we implicitly assume that trivial index translation applies if no appropriate parameters are present. 
/// If this is the case, this method makes and returns the appropriate translator, so the output is always a
/// non-empty shared pointer.
/// @param[in] name corresponding image name to query
/// @return shared pointer to the index translation object
boost::shared_ptr<IUVWeightIndexTranslator> UVWeightParamsHelper::getIndexTranslator(const std::string &name) const
{
   if (indexTranslationIsTrivial(name)) {
       // this is the trivial case
       boost::shared_ptr<GenericUVWeightIndexTranslator> translator(new GenericUVWeightIndexTranslator);
       return translator;
   } 
   // some coefficients are stored
   ASKAPDEBUGASSERT(itsParams);
   const casacore::Vector<float> coeffs = itsParams->valueF("uvweight_indices."+name);
   ASKAPCHECK(coeffs.nelements() == 3, "Expect exactly 3 elements in the parameter vector representing index translation coefficients. Logic error somewhere!");
   ASKAPCHECK(coeffs[0] >= 0.f && coeffs[1] >= 0.f && coeffs[2] >= 0.f, "Index translation coefficients are expected to be non-negative. Logic error somewhere!");
   boost::shared_ptr<GenericUVWeightIndexTranslator> translator(new GenericUVWeightIndexTranslator(casacore::uInt(coeffs[0]), 
          casacore::uInt(coeffs[1]), casacore::uInt(coeffs[2])));
   return translator;
}

/// @brief obtain weight collection
/// @details This method returns collection of uv-weights corresponding to the given name
/// @param[in] name corresponding image name to query
/// @return shared pointer to the new weight collection. 
/// @note Weight collection is designed to be a non-copyable class to ensure stricter control on the operations
/// which are either fundamentally slow or waste memory by creating another copy. Therefore, here a brand new
/// collection object is returned. Its lifetime should be managed by the user.
boost::shared_ptr<UVWeightCollection> UVWeightParamsHelper::getUVWeights(const std::string &name) const
{
   boost::shared_ptr<UVWeightCollection> result(new UVWeightCollection());
   // The second parameter passed to completions() ensures that fixed parameters are also included 
   // (uv-weight related parameters are expected to be fixed as we're not solving for them)
   const std::string baseKey = "uvweight."+name+".";
   ASKAPDEBUGASSERT(itsParams);
   const std::vector<std::string> completions = itsParams->completions(baseKey, true);
   for (std::vector<std::string>::const_iterator ci = completions.begin(); ci != completions.end(); ++ci) {
        const casacore::uInt index = utility::fromString<casacore::uInt>(*ci);
        casacore::Array<float> buf = itsParams->valueF(baseKey + *ci);
        ASKAPCHECK(buf.shape().nelements() == 3, "Each parameter corresponding to uv-weight object should be a 3-dimensional cube, "
                   "you have shape = "<<buf.shape());
        result->add(index, buf);
   }
   return result;
}

/// @brief store weight collection in params under the given name
/// @details This method stores uv-weight collection in the Params object together with
/// associated index translator. At this stage it is assumed that no parameters with the given
/// name exist. If it is necessary to overwrite existing data, one must explicitly remove the old
/// content for the given name (we could've implemented an update, but logic gets a bit non-trivial due to 
/// possibility of orphan parameters which can be used for one collection but not needed for another)
/// @param[in] name corresponding image name to form part of the parameter names
/// @param[in] wts weight collection to store
/// @param[in] translator index translation object. Empty shared pointer implies a trivial translation (i.e. one weight
/// applies to all indices) and cross-check is done that the weight collection contains only one weight at index zero.
void UVWeightParamsHelper::addUVWeights(const std::string &name, const UVWeightCollection &wts, 
                  const boost::shared_ptr<IUVWeightIndexTranslator> &translator) const
{
   ASKAPCHECK(!exists(name), "uv-weight already exists for the image parameter "<<name);
   const std::set<casacore::uInt> indices = wts.indices();
   ASKAPCHECK(indices.size() > 0u, "uv-weight collection given to UVWeightParamsHelper::addUVWeights is empty");
   const std::string baseKey = "uvweight."+name+".";
   ASKAPDEBUGASSERT(itsParams);
   for (std::set<casacore::uInt>::const_iterator ci = indices.begin(); ci != indices.end(); ++ci) {
        // note, we use const-cast here to get native cube interface which is hidden for read-only objects. However,
        // we don't change the cube (the native cube interface was provided for the read-write case because it is needed
        // to be able to write weights more flexibly, this use case is somewhat overlooked but perhaps ok for now to
        // use const-cast but obtain a const-reference instead of exposing this interface for everyone).
        // Although note, that we rely on reference semantics not to cause problems here (e.g. one could edit the array
        // after it is added to Params)
        const casacore::Cube<float> &curWt = const_cast<UVWeightCollection&>(wts).get(*ci);
        const std::string key = baseKey + utility::toString(*ci);
        itsParams->add(key, curWt);
        itsParams->fix(key);
   }
   
   // now sort out the index translation
   if (translator) {
       if (!addIndexTranslator(name, translator) && (indices.size() != 1u || (*indices.begin() != 0u))) {
           // this is essentially a trivial mapping case, but with unmapped data - just give a warning
           // (although may be throwing an exception is more appropriate)
           ASKAPLOG_WARN_STR(logger, "uv-weights index translation is essentially trivial: unnecessary / unaccessible data may be shipped around");
       }
   } else {
      // trivial index translation case, no need to create additional parameters,
      // but we need to do cross-checks that the collection contains only
      // one element with index zero
      ASKAPCHECK(indices.size() == 1u && (*indices.begin() == 0u), 
                 "Expected exactly one uv-weight with zero flat index in the collection passed to addUVWeights for trivial index translation");
   }
}

/// @brief remove all uv-weight parameters associated with the given name
/// @details This method removes all parameters related to uv-weight for the given image name. It
/// is needed if one wants to update uv-weight stored in the associated Params object in a clean way
/// (the number of items in Params differ depending on the content of the uv-weight collection).
/// @param[in] name corresponding image name to delete uv-weights for
void UVWeightParamsHelper::remove(const std::string &name) const
{
   ASKAPDEBUGASSERT(itsParams);
   // The second parameter passed to completions() ensures that fixed parameters are also included 
   // (uv-weight related parameters are expected to be fixed as we're not solving for them)
   const std::string baseKey = "uvweight."+name+".";
   const std::vector<std::string> completions = itsParams->completions(baseKey, true);
   for (std::vector<std::string>::const_iterator ci = completions.begin(); ci != completions.end(); ++ci) {
        itsParams->remove(baseKey + *ci);
   }
   // now remove index translation details if present 
   const std::string indexTranslationKey = "uvweight_indices."+name;
   if (itsParams->has(indexTranslationKey)) {
       itsParams->remove(indexTranslationKey);
   }
}
   
/// @brief check if the trivial index translation is used for the given name
/// @details The absence of "uvweight_indices.name" parameter implies trivial index translation for
/// the parameter with the given name. This method checks the condition.
/// @param[in] name corresponding image name to query
/// @return true if the index translation is trivial for the given name
/// @note this is a low-level detail hidden from the user (hence the interface is not public). 
/// We may need to alter the logic dealing with index translation in the future - only implementing simple case for now.
bool UVWeightParamsHelper::indexTranslationIsTrivial(const std::string &name) const
{
   ASKAPDEBUGASSERT(itsParams);
   return !itsParams->has("uvweight_indices."+name);
}

/// @brief helper method to add index translation parameters
/// @details This is the place where the behaviour of index translator class is decoded and stored in Params.    
/// If a new index translation class is implemented, both this method and getIndexTranslator should be updated.
/// @note we implicitly assume that absence of index translation parameters imply trivial index translation. If more
/// translation classes are implemented this may need to be changed.
/// @param[in] name corresponding image name to form part of the parameter names
/// @param[in] translator index translation object (must not be empty)
/// @return true if the index translation is non trivial and, therefore, if the special parameter with coefficients was added
bool UVWeightParamsHelper::addIndexTranslator(const std::string &name, const boost::shared_ptr<IUVWeightIndexTranslator> &translator) const
{
   ASKAPDEBUGASSERT(translator);
   // we only support GenericUVWeightTranslator at the moment. If more options are implemented in the future this code
   // needs to be expanded
   boost::shared_ptr<GenericUVWeightIndexTranslator> genTranslator = boost::dynamic_pointer_cast<GenericUVWeightIndexTranslator>(translator);
   ASKAPCHECK(genTranslator, "UVWeightParamsHelper::addIndexTranslator only supports GenericUVWeightTranslator at the moment");
   // now extract coefficients to store them in Params. It is a bit hacky way, but it has to be specific to the particular 
   // translator class + I am not sure we'd benefit in any way from having such operation to be a class member as it is Params-specific.
   // We could've implemented proper getter methods, though. But I hope this way of getting this info is not too confusing.
   const casacore::uInt coeffBeam = genTranslator->indexOf(1u,0u,0u);
   const casacore::uInt coeffField = genTranslator->indexOf(0u,1u,0u);
   const casacore::uInt coeffSource = genTranslator->indexOf(0u,0u,1u);
   if (coeffBeam == 0u && coeffField == 0u && coeffSource == 0u) {
       // this is a trivial case, no special parameter is needed
       return false; 
   }
   casacore::Vector<float> buf = {float(coeffBeam), float(coeffField), float(coeffSource)};
   const std::string indexTranslationKey = "uvweight_indices."+name;
   ASKAPDEBUGASSERT(itsParams);
   itsParams->add(indexTranslationKey, buf);
   itsParams->fix(indexTranslationKey);
   return true;
}

/// @brief copy given parameter to another Params object
/// @details This method encapsulates both read and write operations required to copy one weight-related parameter
/// (corresponding to a certain image name) to another Params. 
/// @param[in] dest non-const reference to the destination Params object
/// @param[in] name image name corresponding to the parameter to copy
/// @note Nothing is copied if the given parameter doesn't exist. If the destination has 
/// old data, the appropriate parameters are removed first even if the required parameter is
/// missing in the source and nothing would be copied. At this stage, we're trying to copy by
/// reference if we can (because of double vs. float options it is not always possible).
void UVWeightParamsHelper::copyTo(scimath::Params &dest, const std::string &name) const
{
   UVWeightParamsHelper destHlp(dest);
   // always remove first, if the parameter with the given name doesn't exist in the current object,
   // it won't exist in the destination either - this is a safer behaviour
   destHlp.remove(name);
   if (exists(name)) {
       // we could've used high-level methods to read and write, but doing it from first principles allows us to 
       // cut down on unnecessary operations.
       ASKAPDEBUGASSERT(itsParams);
       // The second parameter passed to completions() ensures that fixed parameters are also included 
       // (uv-weight related parameters are expected to be fixed as we're not solving for them)
       const std::string baseKey = "uvweight."+name+".";
       const std::vector<std::string> completions = itsParams->completions(baseKey, true);
       for (std::vector<std::string>::const_iterator ci = completions.begin(); ci != completions.end(); ++ci) {
            const std::string parName = baseKey + *ci;
            dest.add(parName, itsParams->valueF(parName));
            // we expect all parameters related to uv-weights to always be fixed, so just fix the copied one without checking
            // whether the original one was fixed (although it should be)
            dest.fix(parName);
       }
       // now add index translation details if present 
       const std::string indexTranslationKey = "uvweight_indices."+name;
       if (itsParams->has(indexTranslationKey)) {
           dest.add(indexTranslationKey, itsParams->valueF(indexTranslationKey));
           dest.fix(indexTranslationKey);
       }
   }
}

} // namespace synthesis

} // namespace askap

