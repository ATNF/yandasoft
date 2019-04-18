pipeline {
  agent {
    node {
      label 'ubuntu-bionic'
    }

  }
  stages {
    stage('Get dependencies') {
      steps {
        deleteDir()
        sh '''#provision
sudo apt update;
sudo apt install -y cmake;
sudo apt install -y flex bison;
sudo apt install -y gfortran;
sudo apt install -y git;
sudo apt install -y g++;
sudo apt install -y libboost-dev;
sudo apt install -y libboost-python-dev;
sudo apt install -y libboost-filesystem-dev;
sudo apt install -y libboost-program-options-dev;
sudo apt install -y libboost-signals-dev;    
sudo apt install -y libboost-system-dev; 
sudo apt install -y libboost-thread-dev; 
sudo apt install -y libcfitsio-dev;   
sudo apt install -y libffi-dev;    
sudo apt install -y libfftw3-dev; 
sudo apt install -y libgsl-dev;  
sudo apt install -y liblog4cxx-dev;
sudo apt install -y libopenblas-dev; 
sudo apt install -y libopenmpi-dev; 
sudo apt install -y libpython-dev; 
sudo apt install -y make;         
sudo apt install -y patch;       
sudo apt install -y python-pip; 
sudo apt install -y subversion;
sudo apt install -y wcslib-dev; 
sudo apt-add-repository -s ppa:kernsuite/kern-5;
sudo apt update;
sudo apt install -y casacore
sudo apt install -y casarest
sudo apt install -y casadata

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