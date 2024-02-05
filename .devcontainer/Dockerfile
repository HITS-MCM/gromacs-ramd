ARG VARIANT=ubuntu-22.04
FROM mcr.microsoft.com/vscode/devcontainers/cpp:0-${VARIANT}

LABEL maintainer="Bernd Doser <bernd.doser@h-its.org>"

RUN apt-get update \
 && export DEBIAN_FRONTEND=noninteractive \
 && apt-get -y install --no-install-recommends \
       libomp-14-dev \
       python3-pip \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip \
 && hash -r pip3

RUN pip install \
      dash \
      pandas
