pipeline {
  agent any
  stages {
    stage('Build') {
      steps {
        sh '''if [ -d yandasoft-install ]; then

rm -rf yandasoft-install

fi
'''
        sh 'echo "Pulling the build script from repo...."'
        sh 'git clone https://ord006@bitbucket.csiro.au/scm/askapsdp/yandasoft-install.git'
        sh 'cd yandasoft-install'
        sh 'echo "Done"'
        sh './build_all.sh -s ubuntu -p ${WORKSPACE}/install?'
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