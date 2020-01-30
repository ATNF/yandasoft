. ./ModuleFile
PREFIX=/group/askap/sord/yanda
export CPPFLAGS=" -fPIC "
export CC=cc
export CXX=c++

#../general/build_all.sh -w ${PREFIX}/factory -d -p ${PREFIX} -r -C "-DCMAKE_CXX_FLAGS=-Dcasa=casacore -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/group/askap/sord/yanda -DCFITSIO_ROOT_DIR=/pawsey/cle60up05/devel/PrgEnv-gnu/6.0.4/gcc/7.2.0/ivybridge/cfitsio/3.420 -DWCSLIB_ROOT_DIR=/pawsey/cle60up05/devel/PrgEnv-gnu/6.0.4/gcc/7.2.0/ivybridge/wcslib/5.18 -DFFTW_ROOT=/group/askap/sord/software/cle60up05/apps/PrgEnv-gnu/6.0.4/gcc/7.2.0/ivybridge/fftw/3.3.8 -DCASACORE_ROOT_DIR=/group/askap/sord/yanda"
../general/build_all.sh -x cray -w ${PREFIX}/factory -d -p ${PREFIX} -a -O " -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/group/askap/sord/yanda -DCFITSIO_ROOT_DIR=${CFITSIO_ROOT} -DBUILD_SHARED_LIBS=ON -DWCSLIB_ROOT_DIR=${WCSLIB_ROOT} -DFFTW_ROOT=${FFTW_DIR}/../ -DCASACORE_ROOT_DIR=${CASACORE_ROOT} "
../general/build_all.sh -x cray -w ${PREFIX}/factory -d -p ${PREFIX} -y -O " -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/group/askap/sord/yanda -DCFITSIO_ROOT_DIR=${CFITSIO_ROOT} -DBUILD_SHARED_LIBS=ON -DWCSLIB_ROOT_DIR=${WCSLIB_ROOT} -DFFTW_ROOT=${FFTW_DIR}/../ -DCASACORE_ROOT_DIR=${CASACORE_ROOT} "
#../general/build_all.sh -w ${PREFIX}/factory -d -p ${PREFIX} -y -O "-DCMAKE_CXX_COMPILER=c++ -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS=-Dcasa=casacore -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/group/askap/sord/yanda -DCFITSIO_ROOT_DIR=/pawsey/cle60up05/devel/PrgEnv-gnu/6.0.4/gcc/7.2.0/ivybridge/cfitsio/3.420 -DWCSLIB_ROOT_DIR=/pawsey/cle60up05/devel/PrgEnv-gnu/6.0.4/gcc/7.2.0/ivybridge/wcslib/5.18 -DFFTW_ROOT=/group/askap/sord/software/cle60up05/apps/PrgEnv-gnu/6.0.4/gcc/7.2.0/ivybridge/fftw/3.3.8 -DCASACORE_ROOT_DIR=/group/askap/sord/yanda"
#../general/build_all.sh -d -C "-DCMAKE_CXX_COMPILER=c++ -DCMAKE_C_COMPILER=cc -DCFITSIO_ROOT_DIR=$CFITSIO_ROOT -DDATA_DIR=${PREFIX}/share/casacore/data -DCMAKE_BUILD_TYPE=Release " -r -O "-DCMAKE_CXX_COMPILER=c++ -DCMAKE_C_COMPILER=cc -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCFITSIO_ROOT_DIR=$CFITSIO_ROOT -DWCSLIB_ROOT_DIR=${MAALI_WCSLIB_HOME} -DFFTW_ROOT=${MAALI_FFTW_HOME}  -DCASACORE_ROOT_DIR=${PREFIX}"