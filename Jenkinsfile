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
? echo "casacore directory already exists"
else
??git clone https://github.com/casacore/casacore.git

fi
'''
        }

        dir(path: '/var/lib/jenkins/workspace/yandasoft_development/casacore') {
          sh 'git checkout -b working_copy'
          sh 'git reset --hard d3dad4d'
          sh '''if [ -d build ]; then
? echo "build directory already exists"
else
? mkdir build
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
    stage('Test') {
      steps {
        echo 'Testing....'
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