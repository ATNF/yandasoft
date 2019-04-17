pipeline {
  agent any
  stages {
    stage('Build CASACORE') {
      steps {
        sh '''if [ -d casacore ]; then
? echo "cleaning up"
? rm -rf casacore
fi
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