FROM neurodebian:bionic-non-free

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# Install the validator
RUN apt-get install -y apt-utils curl && \
     curl https://raw.githubusercontent.com/nvm-sh/nvm/v0.35.2/install.sh | bash

RUN export NVM_DIR=$HOME/.nvm && \
     [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh" && \
     [ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion" && \
     nvm install 11.15.0 && \
     nvm use 11.15.0 && \
     nvm alias default 11.15.0 && \
     npm install -g bids-validator

ENV PATH=/root/.nvm/versions/node/v11.15.0/bin:$PATH

# Install Ubuntu dependencies and utilities
RUN apt-get install -y \
      build-essential \
      cmake \
      git \
      graphviz \
      graphviz-dev \
      gsl-bin \
      libcanberra-gtk-module \
      libexpat1-dev \
      libgiftiio-dev \
      libglib2.0-dev \
      libglu1-mesa \
      libglu1-mesa-dev \
      libjpeg-progs \
      libgl1-mesa-dri \
      libglw1-mesa \
      libxml2 \
      libxml2-dev \
      libxext-dev \
      libxft2 \
      libxft-dev \
      libxi-dev \
      libxmu-headers \
      libxmu-dev \
      libxpm-dev \
      libxslt1-dev \
      m4 \
      make \
      mesa-common-dev \
      mesa-utils \
      netpbm \
      pkg-config \
      rsync \
      tcsh \
      unzip \
      vim \
      xvfb \
      xauth \
      zlib1g-dev

# Install 16.04 dependencies
RUN apt-get install -y \
      dh-autoreconf \
      libgsl-dev \
      libmotif-dev \
      libtool \
      libx11-dev \
      libxext-dev \
      x11proto-xext-dev \
      x11proto-print-dev \
      xutils-dev
      
RUN mkdir -p /antsinstall
ENV workingDir=/antsinstall

# Clone the repo
RUN cd ${workingDir} && git clone https://github.com/ANTsX/ANTs.git && cd ANTs && git checkout v2.3.4 && cd -

# Where to build, should be an empty directory
ENV buildDir=${workingDir}/build
ENV installDir=/usr/lib/ants
ENV CMAKE_INSTALL_PREFIX=${installDir}

RUN mkdir -p $buildDir $installDir && cd $buildDir && cmake ${workingDir}/ANTs -DCMAKE_INSTALL_PREFIX=${installDir} && make 2>&1 | tee build.log && cd ANTS-build && make install 2>&1 | tee install.log

ENV ANTSPATH=/usr/lib/ants/bin/
ENV PATH=${ANTSPATH}:$PATH