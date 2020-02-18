#!groovy

pipeline {

  agent {
    label 'docker-gpu-host'
  }

  options {
    timeout(time: 2, unit: 'HOURS')
  }

  stages {
    stage('Build') {
      parallel {
        stage('gcc-8') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.2-gcc-8'
            }
          }
          steps {
            sh '''
              mkdir -p build-gcc-8
              cd build-gcc-8
              cmake -DGMX_BUILD_OWN_FFTW=ON ..
              make 2>&1 |tee make.out
            '''
          }
          post {
            always {
              recordIssues enabledForFailure: true, aggregatingResults: false,
                tool: gcc(id: 'gcc-8', pattern: 'build-gcc-8/make.out')
            }
          }
        }
        stage('clang-6') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.2-clang-6'
            }
          }
          steps {
            sh '''
              mkdir -p build-clang-6
              cd build-clang-6
              cmake -DGMX_BUILD_OWN_FFTW=ON ..
              make 2>&1 |tee make.out
            '''
          }
          post {
            always {
              recordIssues enabledForFailure: true, aggregatingResults: false,
                tool: gcc(id: 'clang-6', pattern: 'build-clang-6/make.out')
            }
          }
        }
        stage('clang-8') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.2-clang-8'
            }
          }
          steps {
            sh '''
              mkdir -p build-clang-8
              cd build-clang-8
              cmake -DGMX_BUILD_OWN_FFTW=ON ..
              make 2>&1 |tee make.out
            '''
          }
          post {
            always {
              recordIssues enabledForFailure: true, aggregatingResults: false,
                tool: gcc(id: 'clang-8', pattern: 'build-clang-8/make.out')
            }
          }
        }
      }
    }
    stage('Test') {
      parallel {
        stage('gcc-8') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.2-gcc-8'
            }
          }
          steps {
            sh 'cd build-gcc-8 && make check'
          }
          post {
            always {
              step([
                $class: 'XUnitPublisher',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-gcc-8/Testing/Temporary/*.xml']]
              ])
            }
          }
        }
        stage('clang-6') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.2-clang-6'
            }
          }
          steps {
            sh 'cd build-clang-6 && make check'
          }
          post {
            always {
              step([
                $class: 'XUnitPublisher',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-clang-6/Testing/Temporary/*.xml']]
              ])
            }
          }
        }
        stage('clang-8') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.2-clang-8'
            }
          }
          steps {
            sh 'cd build-clang-8 && make check'
          }
          post {
            always {
              step([
                $class: 'XUnitPublisher',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-clang-8/Testing/Temporary/*.xml']]
              ])
            }
          }
        }
      }
    }
  }
  post {
    success {
      mail to: 'bernd.doser@h-its.org', subject: "SUCCESS: ${currentBuild.fullDisplayName}", body: "Success: ${env.BUILD_URL}"
    }
    failure {
      mail to: 'bernd.doser@h-its.org', subject: "FAILURE: ${currentBuild.fullDisplayName}", body: "Failure: ${env.BUILD_URL}"
    }
  }
}
