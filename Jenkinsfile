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
        stage('gcc-7') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.1-cmake-3.15-gcc-7'
            }
          }
          steps {
            sh '''
              mkdir -p build-gcc-7
              cd build-gcc-7
              cmake -DCMAKE_BUILD_TYPE=release \
                    -DGMX_BUILD_OWN_FFTW=ON \
                    ..
              make 2>&1 |tee make.out
            '''
          }
          post {
            always {
              recordIssues enabledForFailure: true, aggregatingResults: false,
                tool: gcc(id: 'gcc-7', pattern: 'build-gcc-7/make.out')
            }
          }
        }
        stage('clang-6') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.1-cmake-3.15-clang-6'
            }
          }
          steps {
            sh '''
              mkdir -p build-clang-6
              cd build-clang-6
              cmake -DCMAKE_BUILD_TYPE=release \
                    -DGMX_BUILD_OWN_FFTW=ON \
                    ..
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
      }
    }
    stage('Test') {
      parallel {
        stage('gcc-7') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.1-cmake-3.15-gcc-7'
            }
          }
          steps {
            sh 'cd build-gcc-7 && make check'
          }
          post {
            always {
              step([
                $class: 'XUnitPublisher',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-gcc-7/Testing/Temporary/*.xml']]
              ])
            }
          }
        }
        stage('clang-6') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.1-cmake-3.15-clang-6'
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
