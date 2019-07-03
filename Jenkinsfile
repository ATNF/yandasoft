pipeline {
  agent {
    docker {
      image 'sord/yanda:latest'
    }

  }
  stages {
    stage ('Prepare environment') {
      steps {
        dir(path: '.') {
          sh '''apt install -y libboost-regex-dev'''
        }
      }
    }

    stage('Building base-askap') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-askap ]; then
echo "base-askap directory already exists"
rm -rf base-askap
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-askap.git
cd base-askap
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Building base-logfilters') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-logfilters ]; then
echo "base-logfilters directory already exists"
rm -rf base-logfilters
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-logfilters.git
cd base-logfilters
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Building base-imagemath') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-imagemath ]; then
echo "base-imagemath directory already exists"
rm -rf base-imagemath
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-imagemath.git
cd base-imagemath
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Building base-scimath') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-scimath ]; then
echo "base-scimath directory already exists"
rm -rf base-scimath
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-scimath.git
cd base-scimath
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Building base-askapparallel') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-askapparallel ]; then
echo "base-askapparallel directory already exists"
rm -rf base-askapparallel
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-askapparallel.git
cd base-askapparallel
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Building base-accessors') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-accessors ]; then
echo "base-accessors directory already exists"
rm -rf base-accessors
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-accessors.git
cd base-accessors
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Building base-components') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-components ]; then
echo "base-components directory already exists"
rm -rf base-components
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-components.git
cd base-components
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }

    stage('Building askap-pipelinetasks') {
      steps {
        dir(path: '.') {
          sh '''if [ -d askap-pipelinetasks ]; then
echo "askap-pipelinetasks directory already exists"
rm -rf askap-pipelinetasks
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/askap-pipelinetasks.git
cd askap-pipelinetasks
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
    stage('Building askap-analysis') {
      steps {
        dir(path: '.') {
          sh '''if [ -d askap-analysis ]; then
echo "askap-analysis directory already exists"
rm -rf askap-analysis
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/askap-analysis.git
cd askap-analysis
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }

    stage('Building yandasoft') {
      steps {
        dir(path: '.') {
          sh '''if [ -d build ]; then
echo "yandasoft build directory already exists"
cd build
if [ -f install_manifest.txt ]; then
make uninstall
cd ..
fi
rm -rf build
mkdir build
else
mkdir build
fi'''
        }

        dir(path: 'build') {
          sh '''cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }

      }
    }
  }
  environment {
    WORKSPACE = pwd()
    PREFIX = "${WORKSPACE}/install"
  }
}

