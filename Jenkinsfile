def get_email() {
   dir(path: "${env.REPO}" ) {
      return sh (script:"git log -1 --pretty=format:'%ae'",returnStdout:true).trim()
   }
}  
pipeline {
  agent {
    docker {
      image 'sord/devops:lofar'
    }

  }
  stages {
      stage('Building Dependency (ASKAP)') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-askap ]; then
echo "base-askap directory already exists"
rm -rf base-askap
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-askap.git
cd base-askap
git checkout develop
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }
      }
      }
      stage('Building Dependency (IMAGEMATH)') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-imagemath ]; then
echo "base-imagemath directory already exists"
rm -rf base-imagemath
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-imagemath.git
cd base-imagemath
git checkout develop
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }
      }
      }

      stage('Building Dependency (PARALLEL)') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-askapparallel ]; then
echo "base-askapparallel directory already exists"
rm -rf base-askapparallel
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-askapparallel.git
cd base-askapparallel
git checkout develop
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }
      }
      }
      stage('Building Dependency (SCIMATH)') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-scimath ]; then
echo "base-scimath directory already exists"
rm -rf base-scimath
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-scimath.git
cd base-scimath
git checkout develop
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }
      }
      }
      stage('Building Dependency (ACCESSORS)') {
      steps {
        dir(path: '.') {
          sh '''if [ -d base-accessors ]; then
echo "base-accessors directory already exists"
rm -rf base-accessors
fi
git clone https://bitbucket.csiro.au/scm/askapsdp/base-accessors.git
cd base-accessors
git checkout develop
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''
        }
      }
      }


      stage('Building Debug') {
      steps {
        dir(path: '.') {
          sh '''git fetch --tags
if [ -d build ]; then
echo "build directory already exists"
rm -rf build
fi
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-coverage" ../
make 
'''
        }
      }
    }

    stage('Building Release)') {
      steps {
        dir(path: '.') {
          sh '''git fetch --tags
if [ -d build-release ]; then
echo "build-release directory already exists"
rm -rf build-release
fi
mkdir build-release
cd build-release
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=Release ../
make 
'''
        }
      }
    }    
    stage('Test') {
     steps {
        dir(path: '.') {
          sh '''cd build
ctest -T test --no-compress-output
../askap-cmake/ctest2junit > ctest.xml
          cp ctest.xml $WORKSPACE
'''     }
     }
    }

  }

post {
        always {
             junit 'ctest.xml'
        }
        success {        
             mail to: "${env.EMAIL_TO}",
             from: "jenkins@csiro.au",
             subject: "Succeeded Pipeline: ${currentBuild.fullDisplayName}",
             body: "Build ${env.BUILD_URL} succeeded"
 
        }

        failure {        
             mail to: "${env.EMAIL_TO}",
             from: "jenkins@csiro.au",
             subject: "Failed Pipeline: ${currentBuild.fullDisplayName}",
             body: "Something is wrong with ${env.BUILD_URL}"
 
        }
        unstable {
             mail to: "${env.EMAIL_TO}",
             from: "jenkins@csiro.au",
             subject: "Unstable Pipeline: ${currentBuild.fullDisplayName}",
             body: "${env.BUILD_URL} unstable"
        } 
        changed {
             mail to: "${env.EMAIL_TO}",
             from: "jenkins@csiro.au",
             subject: "Changed Pipeline: ${currentBuild.fullDisplayName}",
             body: "${env.BUILD_URL} changed"
        }
 }
 
  environment {
    
    WORKSPACE = pwd()
    PREFIX = "${WORKSPACE}/install"
    REPO = "${WORKSPACE}/yandasoft/"
    EMAIL_TO = get_email()

  }
}

