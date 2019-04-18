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
sudo apt update

sudo apt install -y cmake ? ? ? ? ? ? ? ? ? ? ? ? `# many` \\
sudo apt install -y?flex bison ? ? ? ? ? ? ? ? ? ?`# casacore` \\
sudo apt install -y?gfortran ? ? ? ? ? ? ? ? ? ? ?`# many` \\
sudo apt install -y?git ? ? ? ? ? ? ? ? ? ? ? ? ? `# many` \\
sudo apt install -y?g++ ? ? ? ? ? ? ? ? ? ? ? ? ? `# many` \\
sudo apt install -y?libboost-dev ? ? ? ? ? ? ? ? ?`# casacore` \\
sudo apt install -y?libboost-python-dev ? ? ? ? ? `# casacore` \\
sudo apt install -y?libboost-filesystem-dev ? ? ? `# yandasoft` \\
sudo apt install -y?libboost-program-options-dev ?`# base-askap` \\
sudo apt install -y?libboost-signals-dev ? ? ? ? ?`# base-askap` \\
sudo apt install -y?libboost-system-dev ? ? ? ? ? `# casarest` \\
sudo apt install -y?libboost-thread-dev ? ? ? ? ? `# casarest` \\
sudo apt install -y?libcfitsio-dev ? ? ? ? ? ? ? ?`# casacore` \\
sudo apt install -y?libffi-dev ? ? ? ? ? ? ? ? ? ?`# cryptography (python) -> paramiko -> daliuge` \\
sudo apt install -y?libfftw3-dev ? ? ? ? ? ? ? ? ?`# casacore` \\
sudo apt install -y?libgsl-dev ? ? ? ? ? ? ? ? ? ?`# many` \\
sudo apt install -y?liblog4cxx-dev ? ? ? ? ? ? ? ?`# yandasoft` \\
sudo apt install -y?libopenblas-dev ? ? ? ? ? ? ? `# casacore` \\
sudo apt install -y?libopenmpi-dev ? ? ? ? ? ? ? ?`# adios, casacore, oskar` \\
sudo apt install -y?libpython-dev ? ? ? ? ? ? ? ? `# casacore, few python packages` \\
sudo apt install -y?make ? ? ? ? ? ? ? ? ? ? ? ? ?`# many` \\
sudo apt install -y?python-pip ? ? ? ? ? ? ? ? ? ?`# so we can pip install virtualenv` \\
sudo apt install -y?subversion ? ? ? ? ? ? ? ? ? ?`# lofar-blob, lofar-common` \\
sudo apt install -y?wcslib-dev ? ? ? ? ? ? ? ? ? ?`# casacore`
sudo apt clean'''
        sh '''git clone https://github.com/casacore/casacore.git


'''
        dir(path: 'casacore') {
          sh '''git checkout -b working_copy
git reset --hard d3dad4d
mkdir build
'''
        }

        dir(path: '${WORKSPACE}/casacore/build') {
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