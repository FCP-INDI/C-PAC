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
RUN apt-get update \
    && apt-get install -y git \
    && mkdir -p /opt/dcan-tools \
    && git clone -b 'v2.0.0' --single-branch --depth 1 https://github.com/DCAN-Labs/DCAN-HCP.git /opt/dcan-tools/pipeline

# use neurodebian runtime as parent image
FROM neurodebian:jammy-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Ubuntu Jammy base image"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# add default user
RUN groupadd -r c-pac \
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
      dc \
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
      python-is-python3 \
      python3 \
      python3-dev \
      python3-pip \
      rsync \
      software-properties-common \
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
RUN pip install -r /opt/requirements.txt

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
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/pip/*

# Set user
USER c-pac_user