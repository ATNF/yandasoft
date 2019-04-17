pipeline {
  agent any
  stages {
    stage('Build CASACORE') {
      steps {
        deleteDir()
        sh '''git clone https://github.com/casacore/casacore.git
git checkout -b working_copy
git reset --hard COMMIT-d3dad4d

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