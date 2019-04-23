pipeline {
  agent {
    docker {
      image 'sord/yanda:latest'
    }

  }
  stages {
    stage('Get Dependencies') {
      steps {
        dir(path: '/var/lib/jenkins/workspace/yandasoft_development') {
          echo 'hello world I am in directory'
          sh 'pwd'
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