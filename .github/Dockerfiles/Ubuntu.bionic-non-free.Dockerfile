# Copyright (C) 2022-2023  C-PAC Developers

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
FROM python:3.11-slim as Python
FROM ghcr.io/fcp-indi/c-pac_templates:latest as c-pac_templates
FROM neurodebian:bionic-non-free AS dcan-hcp


ARG DEBIAN_FRONTEND=noninteractive

# Adding DCAN dependencies & HCP code
RUN apt-get update && \
    apt-get install -y git && \
    mkdir -p /opt/dcan-tools && \
    git clone -b 'v2.0.0' --single-branch --depth 1 https://github.com/DCAN-Labs/DCAN-HCP.git /opt/dcan-tools/pipeline

# using neurodebian runtime as parent image
FROM neurodebian:bionic-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Ubuntu Bionic base image"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# Installing specified Python version
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      software-properties-common \
    && add-apt-repository 'deb http://ports.ubuntu.com/ubuntu-ports jammy main' \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      libc6 \
      libc-bin \
    && add-apt-repository --remove 'deb http://ports.ubuntu.com/ubuntu-ports jammy main' \
    && apt-get update
COPY --from=Python /lib/ /lib/
COPY --from=Python /usr/local/ /usr/local/
RUN ldconfig

# Creating usergroup and user
# Allow users to update / create themselves
# Installing system requirments & the BIDS validator
RUN groupadd -r c-pac && \
    useradd -r -g c-pac c-pac_user && \
    mkdir -p /home/c-pac_user/ && \
    chown -R c-pac_user:c-pac /home/c-pac_user && \
    chmod 777 / && \
    chmod ugo+w /etc/passwd && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
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
      x11proto-xext-dev \
      x11proto-print-dev \
      xauth \
      xutils-dev \
      xvfb \
      zlib1g-dev && \
    export XDG_CONFIG_HOME=/usr/bin && \
    curl https://raw.githubusercontent.com/nvm-sh/nvm/v0.35.2/install.sh | bash && \
    export NVM_DIR=/usr/bin/nvm && \
    [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh" && \
    [ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion" && \
    nvm install 12.12.0 && \
    nvm use 12.12.0 && \
    nvm alias default 12.12.0 && \
    npm install --global npm@^7 && \
    npm install -g bids-validator && \
    curl -sLo /tmp/libpng12.deb http://mirrors.kernel.org/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1_amd64.deb && \
    dpkg -i /tmp/libpng12.deb && \
    rm /tmp/libpng12.deb && \
    add-apt-repository --yes ppa:zeehio/libxp && \
    apt-get update && \
    apt-get install --no-install-recommends \
      libxp6 \
      libxp-dev && \
    add-apt-repository --remove --yes ppa:zeehio/libxp && \
    apt-get update && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG="en_US.UTF-8" && \
    chmod 777 /opt && \
    chmod a+s /opt

ENV PATH=/usr/bin/nvm/versions/node/v12.12.0/bin:$PATH

# Installing torch & Python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip install \
      torch==1.2.0 torchvision==0.4.0 -f https://download.pytorch.org/whl/torch_stable.html && \
    pip install --upgrade setuptools && \
    pip install --upgrade pip && \
    pip install -r /opt/requirements.txt && \
    pip install xvfbwrapper && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
    apt-get install -y --no-install-recommends git-lfs && \
    git lfs install

# Installing C-PAC templates and atlases
COPY --from=c-pac_templates /cpac_templates /cpac_templates
COPY --from=dcan-hcp /opt/dcan-tools/pipeline/global /opt/dcan-tools/pipeline/global
COPY --from=ghcr.io/fcp-indi/c-pac/neuroparc:v1.0-human /ndmg_atlases /ndmg_atlases

# Installing surface files for downsampling
COPY --from=c-pac_templates /opt/dcan-tools/pipeline/global/templates/standard_mesh_atlases/ /opt/dcan-tools/pipeline/global/templates/standard_mesh_atlases/
COPY --from=c-pac_templates /opt/dcan-tools/pipeline/global/templates/Greyordinates/ /opt/dcan-tools/pipeline/global/templates/Greyordinates/


ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user
