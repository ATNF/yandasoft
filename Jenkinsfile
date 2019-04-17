pipeline {
	agent any
	stages {
    	stage('Build') {
        	steps {
			echo 'Pulling the build script from repo....'
			sh 'git clone https://ord006@bitbucket.csiro.au/scm/askapsdp/yandasoft-install.git'
			cd yandasoft-install
			echo 'Done'
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
