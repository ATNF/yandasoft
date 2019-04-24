pipeline {
  agent {
    docker {
      image 'sord/yanda:latest'
    }

  }
  stages {
    stage('Build casacore') {
      steps {
        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          sh '''if [ -d casacore ]; then
echo "casacore directory already exists"
cd casacore
git checkout working_copy
else
git clone https://github.com/casacore/casacore.git
cd casacore
git checkout -b working_copy
fi
'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/casacore') {
          sh 'git reset --hard d3dad4d'
          sh '''if [ -d build ]; then
echo "build directory already exists"
else
mkdir build
fi
'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/casacore/build') {
          sh '''cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX}
make all -j2
make all -j2 install

'''
        }

      }
    }
    stage('Build casarest') {
      steps {
        echo 'casarest ...'
        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          sh '''if [ -d casarest ]; then
echo "casarest directory already exists"
cd casarest
git checkout working_copy
else
git clone https://github.com/casacore/casarest.git
cd casarest
git checkout -b working_copy
fi'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/casarest') {
          sh '''git reset --hard fa01137
if [ -d build ]; then
echo "casarest build directory already exists"
else
mkdir build
fi
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCASACORE_ROOT_DIR=${PREFIX}
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Build LOFAR dependencies') {
      steps {
        echo 'LOFAR....'
        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          sh '''if [ -d lofar-common ]; then
echo "lofar-common already exists"
else
git clone https://bitbucket.csiro.au/scm/askapsdp/lofar-common.git
fi'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/lofar-common') {
          sh '''if [ -d build ]; then
echo "Build directory already exists"
else
mkdir build
fi
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=${PREFIX}
make -j2
make -j2 install
'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          sh '''if [ -d lofar-blob ]; then
echo "lofar-blob already exists"
else
git clone https://bitbucket.csiro.au/scm/askapsdp/lofar-blob.git
fi'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/lofar-blob') {
          sh '''if [ -d build ]; then
echo "Build directory already exists"
else
mkdir build
fi
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=${PREFIX}
make -j2
make -j2 install
'''
        }
        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          sh '''if [ -d base-askap ]; then
echo "base-askap already exists"
else
git clone https://bitbucket.csiro.au/scm/askapsdp/base-askap.git
fi'''
        }
      }
    }
    stage('Build base-askap') {
    steps { 
        
        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          sh '''if [ -d base-askap ]; then
echo "base-askap already exists"
else
git clone https://bitbucket.csiro.au/scm/askapsdp/base-askap.git
fi'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/base-askap') {
          sh '''if [ -d build ]; then
echo "Build directory already exists"
else
mkdir build
fi
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=${PREFIX}
make -j2
make -j2 install
'''

       }
      }
    }
  }
  environment {
    WORKSPACE = '/var/lib/jenkins/workspace'
    PREFIX = '/var/lib/jenkins/workspace/yandasoft_development/install'
  }
}
