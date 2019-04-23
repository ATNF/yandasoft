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
        echo 'Testing....'
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
    stage('Deploy') {
      steps {
        echo 'Deploying....'
      }
    }
  }
  environment {
    WORKSPACE = '/var/lib/jenkins/workspace'
    PREFIX = '/var/lib/jenkins/workspace/yandasoft_development/install'
  }
}