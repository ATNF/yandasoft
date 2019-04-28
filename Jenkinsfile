pipeline {
  agent {
    docker {
      image 'sord/yanda:latest'
    }

  }
  stages {
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
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
make -j2
make -j2 install
'''

       }
      }
    } 
    stage('Build base-logfilters') {
    steps { 
        
        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          sh '''if [ -d base-logfilters ]; then
echo "base-logfilters already exists"
else
git clone https://bitbucket.csiro.au/scm/askapsdp/base-logfilters.git
fi'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/base-logfilters') {
          sh '''if [ -d build ]; then
echo "Build directory already exists"
else
mkdir build
fi
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ../
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
