pipeline {
  agent {
    docker {
      image 'sord/yanda:latest'
      args 'docker run -ti -p 8080:8080 -p 50000:50000 -v /var/run/docker.sock:/var/run/docker.sock --privileged --name jenkins jenkins'
    }

  }
  stages {
    stage('Get Dependencies') {
      steps {
        deleteDir()
        sh '''git clone https://github.com/casacore/casacore.git


'''
        dir(path: 'casacore') {
          sh '''git checkout -b working_copy
git reset --hard d3dad4d
mkdir build
'''
        }

        dir(path: 'casacore/build') {
          sh '''cmake .. -DCMAKE_INSTALL_PREFIX=${PREFIX}
make all -j2
make install -j2
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