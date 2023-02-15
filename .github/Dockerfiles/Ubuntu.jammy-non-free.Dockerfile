# Copyright (C) 2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
FROM ghcr.io/fcp-indi/c-pac_templates:latest as c-pac_templates
FROM neurodebian:jammy-non-free AS dcan-hcp


ARG DEBIAN_FRONTEND=noninteractive

# add DCAN dependencies & HCP code
RUN apt-get update && \
    apt-get install -y git && \
    mkdir -p /opt/dcan-tools && \
    git clone -b 'v2.0.0' --single-branch --depth 1 https://github.com/DCAN-Labs/DCAN-HCP.git /opt/dcan-tools/pipeline

# use neurodebian runtime as parent image
FROM neurodebian:jammy-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Ubuntu Bionic base image"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# install install dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      software-properties-common python3-pip python3-dev \
    # upgrade Python
    # && add-apt-repository ppa:deadsnakes/ppa \
    # && apt-get update \
    # && apt-get install -y --no-install-recommends python3.11 python3-pip python3.11-dev python3.11-venv python3.11-distutils python3.11-gdbm python3.11-tk python3.11-lib2to3 \
    # && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 10 \
    # && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 11 \
    # add default user
    && alias python=python3 \
    && groupadd -r c-pac \
    && useradd -r -g c-pac c-pac_user \
    && mkdir -p /home/c-pac_user/ \
    && chown -R c-pac_user:c-pac /home/c-pac_user \
    && chmod 777 / \
    && chmod ugo+w /etc/passwd \
    # install general dependencies
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      apt-utils \
      apt-transport-https \
      bc \
      build-essential \
      bzip2 \
      ca-certificates \
      cmake \
      curl \
      dh-autoreconf \
      git \
      gnupg \
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
      libgsl-dev \
      libmotif-dev \
      libtool \
      libx11-dev \
      libxext-dev \
      libxft2 \
      libxft-dev \
      libxi-dev \
      libxml2 \
      libxml2-dev \
      libxmu-dev \
      libxmu-headers \
      libxpm-dev \
      libxslt1-dev \
      locales \
      m4 \
      make \
      mesa-common-dev \
      mesa-utils \
      netpbm \
      ninja-build \
      openssh-client \
      pkg-config \
      rsync \
      tcsh \
      unzip \
      vim \
      wget \
      x11proto-dev \
      xauth \
      xutils-dev \
      xvfb \
      zlib1g-dev \
  && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
  && echo $TZ > /etc/timezone \
  && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
  && dpkg-reconfigure --frontend=noninteractive locales \
  && update-locale LANG="en_US.UTF-8" \
  && chmod 777 /opt \
  && chmod a+s /opt

# install Python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip install --upgrade pip setuptools \
    && pip install -r /opt/requirements.txt \
    # install git-lfs
    && curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
    && apt-get install -y --no-install-recommends git-lfs \
    && git lfs install

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user