pipeline {
	agent any
	stages {
    	stage('Build') {
        	steps {
			sh 'echo "Pulling the build script from repo...."'
			sh 'git clone https://ord006@bitbucket.csiro.au/scm/askapsdp/yandasoft-install.git'
			sh 'cd yandasoft-install'
			sh 'echo "Done"'
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
