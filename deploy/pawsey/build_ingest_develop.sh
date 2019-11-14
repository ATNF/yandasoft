. ./ModuleFile-ingest
export log4cxx_ROOT_DIR=${log4cxx_ROOT_DIR}
export CPPUNIT_ROOT_DIR=${MAALI_CPPUNIT_HOME}
export XERCES_LIB=${MAALI_XERCESC_HOME}/lib/libxerces-c.so
export XERCES_INC=${MAALI_XERCESC_HOME}/include
export LMDB_HOME=${MAALI_LMDB_HOME}

../general/build_all.sh -d -p /group/askap/sord/ingest/yanda -C "-DCMAKE_INSTALL_PREFIX=/group/askap/sord/ingest/yanda -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DBUILD_PYTHON=OFF -DCFITSIO_ROOT_DIR=${MAALI_CFITSIO_HOME} -DWCSLIB_ROOT_DIR=${MAALI_WCSLIB_HOME} -DFFTW_ROOT=${MAALI_FFTW_HOME}"
../general/build_all.sh -d -p /group/askap/sord/ingest/yanda -r -a -y -O "-DCMAKE_INSTALL_PREFIX=/group/askap/sord/ingest/yanda -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCFITSIO_ROOT_DIR=${MAALI_CFITSIO_HOME} -DWCSLIB_ROOT_DIR=${MAALI_WCSLIB_HOME} -DFFTW_ROOT=${MAALI_FFTW_HOME} -DXercesC_INCLUDE_DIR=${XERCES_INC} -DXercesC_LIBRARY=${XERCES_LIB} -DCASACORE_ROOT_DIR=/group/askap/sord/ingest/yanda"
../general/build_all.sh -d -p /group/askap/sord/ingest/yanda -e -O "-DCMAKE_INSTALL_PREFIX=/group/askap/sord/ingest/yanda -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCFITSIO_ROOT_DIR=${MAALI_CFITSIO_HOME} -DWCSLIB_ROOT_DIR=${MAALI_WCSLIB_HOME} -DFFTW_ROOT=${MAALI_FFTW_HOME} -DXercesC_INCLUDE_DIR=${XERCES_INC} -DXercesC_LIBRARY=${XERCES_LIB} -DCASACORE_ROOT_DIR=/group/askap/sord/ingest/yanda"
../general/build_all.sh -d -p /group/askap/sord/ingest/yanda -i -O "-DCMAKE_INSTALL_PREFIX=/group/askap/sord/ingest/yanda -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCFITSIO_ROOT_DIR=${MAALI_CFITSIO_HOME} -DWCSLIB_ROOT_DIR=${MAALI_WCSLIB_HOME} -DFFTW_ROOT=${MAALI_FFTW_HOME} -DXercesC_INCLUDE_DIR=${XERCES_INC} -DXercesC_LIBRARY=${XERCES_LIB} -DCASACORE_ROOT_DIR=/group/askap/sord/ingest/yanda"
