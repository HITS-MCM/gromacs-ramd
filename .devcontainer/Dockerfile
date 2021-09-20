# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.187.0/containers/cpp/.devcontainer/base.Dockerfile

#FROM nvidia/cuda:10.2-devel
FROM braintwister/ubuntu-20.04-cuda-devel-11.3-gcc-8

RUN conan config set general.revisions_enabled=1 \
 && conan profile new default --detect \
 && conan profile update settings.compiler=clang default \
 && conan profile update settings.compiler.libcxx=libstdc++11 default \
 && conan profile update env.CC=/usr/bin/clang default \
 && conan profile update env.CXX=/usr/bin/clang++ default

RUN mkdir /test
