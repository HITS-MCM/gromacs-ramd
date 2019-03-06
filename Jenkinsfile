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
              image 'braintwister/ubuntu-18.04-cuda-10.0-cmake-3.13-gcc-7'
            }
          }
          steps {
            sh '''
              mkdir -p build-gcc-7
              cd build-gcc-7
              cmake -DCMAKE_BUILD_TYPE=release \
                    -DGMX_BUILD_OWN_FFTW=ON \
                    ..
              make
            '''
          }
          post {
            always {
              step([
                $class: 'WarningsPublisher', canComputeNew: false, canResolveRelativePaths: false,
                defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '',
                parserConfigurations: [[parserName: 'GNU Make + GNU C Compiler (gcc)', pattern: 'build-gcc-7/make.out']],
                unHealthy: ''
              ])
            }
          }
        }
        stage('clang-6') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.0-cmake-3.13-clang-6'
            }
          }
          steps {
            sh '''
              mkdir -p build-clang-6
              cd build-clang-6
              cmake -DCMAKE_BUILD_TYPE=release \
                    -DGMX_BUILD_OWN_FFTW=ON \
                    ..
              make
            '''
          }
          post {
            always {
              step([
                $class: 'WarningsPublisher', canComputeNew: false, canResolveRelativePaths: false,
                defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '',
                parserConfigurations: [[parserName: 'Clang (LLVM based)', pattern: 'build-clang-6/make.out']],
                unHealthy: ''
              ])
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
              image 'braintwister/ubuntu-18.04-cuda-10.0-cmake-3.13-gcc-7'
            }
          }
          steps {
            sh 'cd build-gcc-7 && make check'
          }
          post {
            always {
              step([
                $class: 'XUnitBuilder',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-gcc-7/Testing/*.xml']]
              ])
            }
          }
        }
        stage('clang-6') {
          agent {
            docker {
              reuseNode true
              image 'braintwister/ubuntu-18.04-cuda-10.0-cmake-3.13-clang-6'
            }
          }
          steps {
            sh 'cd build-clang-6 && make check'
          }
          post {
            always {
              step([
                $class: 'XUnitBuilder',
                thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
                tools: [[$class: 'GoogleTestType', pattern: 'build-clang-6/Testing/*.xml']]
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
