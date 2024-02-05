#!groovy

pipeline {

  agent {
    label 'docker-gpu-host'
  }

  options {
    timeout(time: 2, unit: 'HOURS')
    disableConcurrentBuilds()
  }

  stages {
    stage('Build') {
      parallel {
        stage('gcc') {
          agent {
            docker {
              reuseNode true
              image 'mcr.microsoft.com/vscode/devcontainers/cpp:0-ubuntu-22.04'
            }
          }
          steps {
            sh '''
              cmake -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -B build-gcc
              cmake --build build-gcc 2>&1 |tee build-gcc/make.out
            '''
          }
          post {
            always {
              recordIssues enabledForFailure: true, aggregatingResults: false,
                tool: gcc(id: 'gcc', pattern: 'build-gcc/make.out')
            }
          }
        }
        stage('clang') {
          agent {
            docker {
              reuseNode true
              image 'mcr.microsoft.com/vscode/devcontainers/cpp:0-ubuntu-20.04'
            }
          }
          steps {
            sh '''
              cmake -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -B build-clang
              cmake --build build-clang 2>&1 |tee build-clang/make.out
            '''
          }
          post {
            always {
              recordIssues enabledForFailure: true, aggregatingResults: false,
                tool: gcc(id: 'clang', pattern: 'build-clang/make.out')
            }
          }
        }
      }
    }
    stage('Test') {
      parallel {
        stage('gcc') {
          agent {
            docker {
              reuseNode true
              image 'mcr.microsoft.com/vscode/devcontainers/cpp:0-ubuntu-22.04'
            }
          }
          steps {
            sh 'cmake --build build-gcc --target check'
          }
          post {
            always {
              step([
                $class: 'XUnitPublisher',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-gcc/Testing/Temporary/*.xml']]
              ])
            }
          }
        }
        stage('clang') {
          agent {
            docker {
              reuseNode true
              image 'mcr.microsoft.com/vscode/devcontainers/cpp:0-ubuntu-20.04'
            }
          }
          steps {
            sh 'cmake --build build-clang --target check'
          }
          post {
            always {
              step([
                $class: 'XUnitPublisher',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-clang/Testing/Temporary/*.xml']]
              ])
            }
          }
        }
      }
    }
  }
  post {
    failure {
      mail to: 'bernd.doser@h-its.org', subject: "FAILURE: ${currentBuild.fullDisplayName}", body: "Failure: ${env.BUILD_URL}"
    }
  }
}
