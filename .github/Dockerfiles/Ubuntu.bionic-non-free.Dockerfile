# using neurodebian runtime as parent image
FROM neurodebian:bionic-non-free

ARG DEBIAN_FRONTEND=noninteractive

# create usergroup and user
RUN groupadd -r c-pac && \
    useradd -r -g c-pac c-pac_user && \
    mkdir -p /home/c-pac_user/ && \
    chown -R c-pac_user:c-pac /home/c-pac_user

RUN apt-get update

# Install the BIDS validator
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
      bzip2 \
      ca-certificates \
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
      locales \
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
      connectome-workbench \
      dh-autoreconf \
      libgsl-dev \
      libmotif-dev \
      libtool \
      libx11-dev \
      libxext-dev \
      x11proto-xext-dev \
      x11proto-print-dev \
      xutils-dev

# Install libpng12
RUN curl -sLo /tmp/libpng12.deb http://mirrors.kernel.org/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1_amd64.deb && \
    dpkg -i /tmp/libpng12.deb && \
    rm /tmp/libpng12.deb

# Install libxp from third-party repository
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository --yes ppa:zeehio/libxp && \
    apt-get update && apt-get install libxp6 libxp-dev && \
    add-apt-repository --remove --yes ppa:zeehio/libxp && \
    apt-get update

# install miniconda
RUN curl -sO https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh && \
    bash Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-py37_4.8.2-Linux-x86_64.sh

# update path to include conda
ENV PATH=/usr/local/miniconda/bin:$PATH

# install conda dependencies
RUN conda update conda -y && \
    conda install nomkl && \
    conda install -y  \
        blas \
        cython \
        matplotlib==3.1.3 \
        networkx==2.4 \
        nose==1.3.7 \
        numpy==1.16.4 \
        pandas==0.23.4 \
        scipy==1.4.1 \
        traits==4.6.0 \
        wxpython \
        pip

# install torch
RUN pip install torch==1.2.0 torchvision==0.4.0 -f https://download.pytorch.org/whl/torch_stable.html

# install python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip install --upgrade setuptools
RUN pip install --upgrade pip
RUN pip install -r /opt/requirements.txt
RUN pip install xvfbwrapper

# install PyPEER
RUN pip install git+https://github.com/ChildMindInstitute/PyPEER.git

# install cpac templates
ADD dev/docker_data/cpac_templates.tar.gz /

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs
RUN git lfs install

# Get atlases
RUN mkdir -p /ndmg_atlases/label && \
    GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/neurodata/neuroparc.git /tmp/neuroparc && \
    cd /tmp/neuroparc && \
    git lfs install --skip-smudge && \
    git lfs pull -I "atlases/label/Human/*" && \
    cp -r /tmp/neuroparc/atlases/label/Human /ndmg_atlases/label && \
    cd -

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user