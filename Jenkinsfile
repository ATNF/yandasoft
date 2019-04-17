pipeline {
  agent any
  stages {
    stage('Build CASACORE') {
      steps {
        sh '''if [ -d casacore ]; then?
? echo "cleaning up"
? rm -rf casacore
fi
'''
        sh '''git clone https://github.com/casacore/casacore.git
git checkout -b working_copy
git reset --hard d3dad4d

'''
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
}