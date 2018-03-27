#!groovy

pipeline {

  agent {
    label 'docker'
  }

  options {
    timeout(time: 2, unit: 'HOURS')
  }

  stages {
    stage('Build') {
      parallel {
        stage('gcc-5') {
          agent {
            docker {
              reuseNode true
              image 'bernddoser/docker-devel-cpp:ubuntu-16.04-cmake-3.10-gcc-5-gtest-1.8.0'
            }
          }
          steps {
            mkdir -p build-gcc-5
            cd build-gcc-5
            cmake -DCMAKE_BUILD_TYPE=release \
                  -DGMX_BUILD_OWN_FFTW=ON \
                  ..
            make
          }
          post {
            always {
              step([
                $class: 'WarningsPublisher', canComputeNew: false, canResolveRelativePaths: false,
                defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '',
                parserConfigurations: [[parserName: 'GNU Make + GNU C Compiler (gcc)', pattern: 'build-gcc-5/make.out']],
                unHealthy: ''
              ])
            }
          }
        }
        stage('clang-4') {
          agent {
            docker {
              reuseNode true
              image 'bernddoser/docker-devel-cpp:ubuntu-16.04-cmake-3.10-clang-4-gtest-1.8.0'
            }
          }
          steps {
            mkdir -p build-clang-4
            cd build-clang-4
            cmake -DCMAKE_BUILD_TYPE=release \
                  -DGMX_BUILD_OWN_FFTW=ON \
                  ..
            make
          }
          post {
            always {
              step([
                $class: 'WarningsPublisher', canComputeNew: false, canResolveRelativePaths: false,
                defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '',
                parserConfigurations: [[parserName: 'Clang (LLVM based)', pattern: 'build-clang-4/make.out']],
                unHealthy: ''
              ])
            }
          }
        }
      }
    }
    stage('Test') {
      parallel {
        stage('gcc-5') {
          agent {
            docker {
              reuseNode true
              image 'bernddoser/docker-devel-cpp:ubuntu-16.04-cmake-3.10-gcc-5-gtest-1.8.0'
            }
          }
          steps {
            sh 'cd build-gcc-5 && make test'
          }
          post {
            always {
              step([
                $class: 'XUnitBuilder',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-gcc-5/Testing/*.xml']]
              ])
            }
          }
        }
        stage('clang-4') {
          agent {
            docker {
              reuseNode true
              image 'bernddoser/docker-devel-cpp:ubuntu-16.04-cmake-3.10-clang-4-gtest-1.8.0'
            }
          }
          steps {
            sh 'cd build-clang-4 && make test'
          }
          post {
            always {
              step([
                $class: 'XUnitBuilder',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-clang-4/Testing/*.xml']]
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
