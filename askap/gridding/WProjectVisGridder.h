/// @file WProjectVisGridder.h
///
/// WProjectVisGridder: W projection gridding

/// @copyright (c) 2007,2016 CSIRO
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///
#ifndef ASKAP_SYNTHESIS_WPROJECTVISGRIDDER_H_
#define ASKAP_SYNTHESIS_WPROJECTVISGRIDDER_H_

// ASKAPsoft includes
#include <askap/gridding/WDependentGridderBase.h>

// Local package includes
#include <askap/dataaccess/IConstDataAccessor.h>

namespace askap
{
    namespace synthesis
    {
        /// helper function
        /// @brief a helper method for a ref copy of casa arrays held in stl vector
        /// @param[in] in input array
        /// @param[out] out output array (will be resized)
        /// @return size of the cache in bytes (assuming Complex array elements)
        template<typename T>
        size_t deepRefCopyOfSTDVector(const std::vector<T> &in, std::vector<T> &out)
        {
            out.resize(in.size());
            size_t total = 0;
            const typename std::vector<T>::const_iterator inEnd = in.end();
            typename std::vector<T>::iterator outIt = out.begin();
            for (typename std::vector<T>::const_iterator inIt = in.begin();
                inIt != inEnd; ++inIt,++outIt) {
                outIt->reference(*inIt);
                total += outIt->nelements()*sizeof(casa::Complex)+sizeof(T);
            }
            return total;
        }


        /// @brief Visibility gridder using W projection
        /// @details The visibilities are gridded using a convolution
        /// function that implements a Fresnel transform. This corrects
        /// for the w term in the full synthesis measurement equation.
        ///
        /// The convolution function is calculated straightforwardly
        /// by constructing an image of the complex w-dependent
        /// phasor and Fourier transforming. The calculation is
        /// done using a coarse but large grid in image space so
        /// that it is sub-sampled in uv space.
        ///
        /// The scaling is slow in data points, fast in w planes.
        ///
        /// @ingroup gridding
        class WProjectVisGridder : public WDependentGridderBase
        {
            public:
                /// @brief Construct a gridder for W projection
                /// @param wmax Maximum baseline (wavelengths)
                /// @param nwplanes Number of w planes
                /// @param cutoff Cutoff in determining support e.g. 10^-3 of the peak
                /// @param overSample Oversampling (currently limited to <=1)
                /// @param maxSupport Maximum support to be allowed
                /// @param limitSupport Upper limit of support
                /// @param name Name of table to save convolution function into
                /// @param alpha spheroidal function alpha value
                /// @param useDouble Use Double precision for the convolution functions
                WProjectVisGridder(const double wmax, const int nwplanes,
                        const double cutoff, const int overSample,
                        const int maxSupport, const int limitSupport,
                        const std::string& name=std::string(""),
                        const float alpha=1.,
                        const bool shareCF=false);

                virtual ~WProjectVisGridder();

                /// @brief copy constructor
                /// @details It is required to decouple internal arrays in the input
                /// object and the copy.
                /// @param[in] other input object
                WProjectVisGridder(const WProjectVisGridder &other);

                /// Clone a copy of this Gridder
                virtual IVisGridder::ShPtr clone();

                /// @brief static method to get the name of the gridder
                /// @details We specify parameters per gridder type in the parset file.
                /// This method returns the gridder name which should be used to extract
                /// a subset of parameters for createGridder method.
                static inline std::string gridderName() { return "WProject";}

                /// @brief static method to create gridder
                /// @details Each gridder should have a static factory method, which is
                /// able to create a particular type of the gridder and initialise it with
                /// the parameters taken from the given parset. It is assumed that the
                /// method receives a subset of parameters where the gridder name is already
                /// taken out.
                /// @param[in] parset input parset file
                /// @return a shared pointer to the gridder instance
                static IVisGridder::ShPtr createGridder(const LOFAR::ParameterSet& parset);

            protected:
                /// @brief additional operations to configure gridder
                /// @details This method is supposed to be called from createGridder and could be
                /// used in derived classes to avoid too much duplication of the code. For this
                /// particular class it configures variable/offset support and cutoff behavior.
                /// @param[in] parset input parset file
                virtual void configureGridder(const LOFAR::ParameterSet& parset);

                /// @brief obtain buffer used to create convolution functions
                /// @return a reference to the buffer held as a shared pointer
                casacore::Matrix<imtypeComplex> getCFBuffer() const;

                /// @brief initialise buffer for full-sized convolution function
                /// @param[in] uSize size in U
                /// @param[in] vSize size in V
                inline void initCFBuffer(casacore::uInt uSize, casacore::uInt vSize)
                {
                    itsCFBuffer.reset(new casacore::Matrix<imtypeComplex>(uSize,vSize));
                }

                /// @brief initialise sum of weights
                /// @details We keep track the number of times each convolution function is used per
                /// channel and polarisation (sum of weights). This method is made virtual to be able
                /// to do gridder specific initialisation without overriding initialiseGrid.
                /// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
                /// been initialised by the time this method is called.
                virtual void initialiseSumOfWeights();

                /// @brief Initialise the indices
                /// @param[in] acc const data accessor to work with
                virtual void initIndices(const accessors::IConstDataAccessor& acc);

                /// Offset into convolution function
                /// @param row Row number
                /// @param pol Polarisation
                /// @param chan Channel number
                virtual int cIndex(int row, int pol, int chan);

                /// Initialize convolution function
                /// @param[in] acc const data accessor to work with
                virtual void initConvolutionFunction(const accessors::IConstDataAccessor& acc);

                /// @brief largest possible support size
                /// @details Gridders search for support starting from the size returned by this method.
                /// @return maximum possible support of convolution function
                inline int maxSupport() const { return itsMaxSupport; }

                /// @brief truncate support, if necessary
                /// @details This method encapsulates all usage of itsLimitSupport. It truncates the support
                /// if necessary and reports the new value back.
                /// @param[in] support support size to truncate according to itsLimitSupport
                /// @return support size to use (after possible truncation)
                int limitSupportIfNecessary(int support) const;

                /// @brief a structure describing the region of the CF
                /// @details Define the region of significant power in CF using support + offsets.
                struct CFSupport {
                    /// @brief support size
                    int itsSize;
                    /// @brief offset in u of the centre w.r.t. the centre of the array (i.e. centred gaussian would have 0)
                    int itsOffsetU;
                    /// @brief offset in v of the centre w.r.t. the centre of the array (i.e. centred gaussian would have 0)
                    int itsOffsetV;
                    /// @brief constructor to simplify making the class
                    /// @param[in] size support size
                    /// @param[in] u offset in u
                    /// @param[in] v offset in v
                    explicit CFSupport(int size, int u = 0, int v = 0) : itsSize(size), itsOffsetU(u), itsOffsetV(v) {}
                };

                /// @brief search for support parameters
                /// @details This method encapsulates support search operation, taking into account the
                /// cutoff parameter and whether or not an offset is allowed.
                /// @param[in] cfPlane const reference to 2D plane with the convolution function
                /// @return an instance of CFSupport with support parameters
                CFSupport extractSupport(const casacore::Matrix<casacore::DComplex> &cfPlane) const;

                /// @brief search for support parameters
                /// @details This method encapsulates support search operation, taking into account the
                /// cutoff parameter and whether or not an offset is allowed.
                /// @param[in] cfPlane const reference to 2D plane with the convolution function
                /// @return an instance of CFSupport with support parameters
                CFSupport extractSupport(const casacore::Matrix<casacore::Complex> &cfPlane) const;

                /// @brief support is plane-dependent?
                /// @return true, if support should be searched individually for every CF cache plane
                inline bool isSupportPlaneDependent() const { return itsPlaneDependentCFSupport; }

                /// @brief configure support search
                /// @param[in] flag true to search for plane-dependent support, false (default) otherwise
                inline void planeDependentSupport(bool flag) { itsPlaneDependentCFSupport = flag; }

                /// @brief support can be offset?
                inline bool isOffsetSupportAllowed() const { return itsOffsetSupportAllowed; }

                /// @brief configure offset support option
                /// @param[in] flag true to search for plane-dependent support, false (default) otherwise
                inline void offsetSupport(bool flag) { itsOffsetSupportAllowed = flag; }

                /// Mapping from row, pol, and channel to planes of convolution function
                casacore::Cube<int> itsCMap;

                /// @brief obtain absolute cutoff flag
                /// @return true if itsCutoff is an absolute cutoff rather than relative to the peak
                inline bool isCutoffAbsolute() const { return itsCutoffAbs;}

                /// @brief set absolute cutoff flag
                /// @param[in] flag true, if cutoff should be treated as an absolute value
                inline void setAbsCutoffFlag(const bool flag) { itsCutoffAbs = flag; }

                /// @brief Are we using the shared CF cache?
                /// @return true if the shared cache is used
                inline bool shareCF() const { return itsShareCF; }

                /// @brief setup and intialise the variables to be used by the helper methods below.
                /// @param[in] nx, ny - size of full size convolution function grid
                /// @param[in] qnx, qny - size of filled CF grid
                /// @param[in] ccellx, celly - UV cell size after oversampling
                /// @param[in] ccfx, ccfy - Spheroidal weighting for grid
                void setup(int& nx, int& ny, int& qnx, int& qny,
                           double& ccellx, double& ccelly,
                           casacore::Vector<float>& ccfx,
                           casacore::Vector<float>& ccfy);
                /// @brief generates the convolution cache for the PSF gridder
                void doPCFGridder();
                /// @brief generates the convolution cache for the non PSF gridder
                /// @param[in] startlane - the nwplane where the CF generation starts from
                /// @param[in] endPlane - the nwplane where the CF generation finishes
                void generate(int startPlane, int endPlane);
                /// @brief helper method. save/copy the itsConvFunc to/from the convoultion cache
                void save();
                /// @brief helper methods. calculate the support for a given plane
                WProjectVisGridder::CFSupport calcSupport(const casacore::Matrix<casacore::DComplex> &cfPlane,int iw);
                WProjectVisGridder::CFSupport calcSupport(const casacore::Matrix<casacore::Complex> &cfPlane,int iw);
                /// @brief helper method. calculate the convolution function for this plane.
                /// @param[in] cfPlane - buffer for CF calculation
                /// @param[in] qnx, qny - size of filled CF grid
                /// @param[in] nx, ny - size of full size convolution function grid
                /// @param[in] ccellx, celly - UV cell size after oversampling
                /// @param[in] w - w-value for current plane
                void populateThisPlane(casacore::Matrix<casacore::Complex> &cfPlane,
                                       const int qnx, const int qny,
                                       const int nx, const int ny,
                                       const double ccellx, const double ccelly,
                                       const double w, const casacore::Vector<float>& ccfx,
                                       const casacore::Vector<float>& ccfy);
                /// @brief helper method. store the convolution function to the itsConvFunc vector.
                void populateItsConvFunc(const casacore::Matrix<casacore::Complex> &thisPlane, const int iw, 
                                         const int support, const CFSupport& cfSupport, const int cSize, 
                                         const int nx, const int ny);
                void normalise(std::vector<casacore::Matrix<casacore::Complex> >& convFunc);
                /// @brief - helper method. check if the given support is overflowing
                void calcConvFuncOverflow(const int support, const CFSupport& cfSupport, 
                                          const int nx, const int ny, int& overflow);
                                

                /// @brief Specify if we are using the shared CF cache
                /// @param[in] true if the shared cache is used
                inline void setShareCF(bool flag) { itsShareCF = flag; }

                /// @brief Are we using the shared CF cache?
                bool itsShareCF;

                /// @brief cached CF
                static std::vector<casa::Matrix<casa::Complex> > theirCFCache;

                /// @brief cached CF offsets
                static std::vector<std::pair<int,int> > theirConvFuncOffsets;


            private:
                /// @brief assignment operator
                /// @details Defined as private, so it can't be called (to enforce usage of the
                /// copy constructor
                /// @param[in] other input object
                /// @return reference to itself
                WProjectVisGridder& operator=(const WProjectVisGridder &other);

            protected:
                /// Maximum support
                int itsMaxSupport;

                /// Threshold for cutoff of convolution function
                double itsCutoff;

                /// Upper limit of support
                int itsLimitSupport;

                /// @brief true to search for plane-dependent support
                /// @details itsSupport is the support for the first plane in this case (usually the largest)
                bool itsPlaneDependentCFSupport;

                /// @brief true if the support can be offset
                /// @details If this parameter is true, offset convolution functions will be built.
                bool itsOffsetSupportAllowed;

                /// @brief buffer for full-sized convolution function
                /// @details We have to calculate convolution functions on a larger grid and then cut out
                /// a limited support out of it. Mosaicing gridders may need to compute a significant number
                /// of convolution functions. To speed things up, the allocation of the buffer is taken
                /// outside initConvolutionFunction method. A shared pointer to this buffer is held as a
                /// data member as initialisation and usage happen in different methods of this class.
                boost::shared_ptr<casacore::Matrix<imtypeComplex> > itsCFBuffer;

                /// @brief itsCutoff is an absolute cutoff, rather than relative to the peak of a particular CF plane
                bool itsCutoffAbs;
        };
    }
}
#endif
