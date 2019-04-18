pipeline {
  agent {
    node {
      label 'ubuntu-bionic'
    }

  }
  stages {
    stage('Build CASACORE') {
      steps {
        deleteDir()
        sh '''#provision
sudo apt update;
sudo apt install -y cmake;
sudo apt install -y flex bison;
sudo apt install -y gfortran;
sudo apt install -y git;
sudo apt install -y g++;
sudo apt install libboost-dev;

'''
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
    PREFIX = '${WORKSPACE}/install'
  }
}