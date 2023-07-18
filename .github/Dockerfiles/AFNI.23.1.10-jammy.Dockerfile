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
FROM ghcr.io/fcp-indi/c-pac/fsl:6.0.6.5-jammy as FSL
FROM ghcr.io/fcp-indi/c-pac/ubuntu:jammy-non-free as AFNI
USER root
# To use the same Python environment to share common libraries
COPY --from=FSL /usr/share/fsl/6.0 /usr/share/fsl/6.0
ENV FSLDIR=/usr/share/fsl/6.0 \
    PATH=/usr/share/fsl/6.0/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# install AFNI
COPY dev/docker_data/required_afni_pkgs.txt /opt/required_afni_pkgs.txt
COPY dev/docker_data/checksum/AFNI.23.1.10.sha384 /tmp/AFNI.23.1.10.sha384
RUN apt-get update \
    && apt-get install -y \
      apt-transport-https \
      apt-utils \
      bc \
      build-essential \
      bzip2 \
      ca-certificates \
      cmake \
      curl \
      dh-autoreconf \
      eog \
      evince \
      firefox \
      gcc \
      gedit \
      git \
      gnome-terminal \
      gnome-tweaks \
      gnupg \
      graphviz \
      graphviz-dev \
      gsl-bin \
      libcanberra-gtk-module \
      libcurl4-openssl-dev \
      libexpat1-dev \
      libgdal-dev \
      libgfortran-11-dev \
      libgiftiio-dev \
      libgl1-mesa-dri \
      libglib2.0-dev \
      libglu1-mesa \
      libglu1-mesa-dev \
      libglw1-mesa \
      libglw1-mesa-dev \
      libgomp1 \
      libgsl-dev \
      libjpeg-progs \
      libjpeg62 \
      libmotif-dev \
      libnode-dev \
      libopenblas-dev \
      libssl-dev \
      libtool \
      libudunits2-dev \
      libx11-dev \
      libxext-dev \
      libxft-dev \
      libxft2 \
      libxi-dev \
      libxm4 \
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
      nautilus \
      netpbm \
      ninja-build \
      openssh-client \
      pkg-config \
      r-base-dev \
      rsync \
      software-properties-common \
      tcsh \
      unzip \
      vim \
      wget \
      x11proto-xext-dev \
      xauth \
      xfonts-100dpi \
      xfonts-base \
      xterm \
      xutils-dev \
      xvfb \
      zlib1g-dev && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG="en_US.UTF-8" \
    && ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.27 /usr/lib/x86_64-linux-gnu/libgsl.so.19 \
    && AFNI_VERSION="23.1.10" \
    && curl -LOJ https://github.com/afni/afni/archive/AFNI_${AFNI_VERSION}.tar.gz \
    && sha384sum --check /tmp/AFNI.23.1.10.sha384 \
    && mkdir /opt/afni \
    && tar -xvf afni-AFNI_${AFNI_VERSION}.tar.gz -C /opt/afni --strip-components 1 \
    && rm -rf afni-AFNI_${AFNI_VERSION}.tar.gz \
    # Fix GLwDrawA per https://github.com/afni/afni/blob/AFNI_23.1.10/src/other_builds/OS_notes.linux_fedora_25_64.txt
    && cd /usr/include/GL \
    && mv GLwDrawA.h GLwDrawA.h.orig \
    && sed 's/GLAPI WidgetClass/extern GLAPI WidgetClass/' GLwDrawA.h.orig > GLwDrawA.h \
    && cd /opt/afni/src \
    && sed '/^INSTALLDIR =/c INSTALLDIR = /opt/afni' other_builds/Makefile.linux_ubuntu_22_64 > Makefile \
    && make vastness && make cleanest \
    && cd /opt/afni \
    # filter down to required packages
    ls > full_ls \
    && sed 's/linux_openmp_64\///g' /opt/required_afni_pkgs.txt | sort > required_ls \
    && comm -2 -3 full_ls required_ls | xargs rm -rf full_ls required_ls \
    # get rid of stuff we just needed for building
    && apt-get remove -y \
      apt-transport-https \
      apt-utils \
      build-essential \
      bzip2 \
      ca-certificates \
      cmake \
      curl \
      dh-autoreconf \
      evince \
      firefox \
      gedit \
      git \
      gnome-terminal \
      gnome-tweaks \
      libglw1-mesa-dev \
      m4 \
      make \
      ninja-build \
      openssh-client \
      pkg-config \
      unzip \
      wget \
      xterm \
    && ldconfig \
    && rm -rf /opt/afni/src

# set up AFNI
ENV PATH=/opt/afni:$PATH

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/*

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
AFNI 23.1.10 (Publius Helvius Pertinax) stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=AFNI /lib/x86_64-linux-gnu/ld* /lib/x86_64-linux-gnu/
COPY --from=AFNI /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
COPY --from=AFNI /lib64/ld* /lib64/
COPY --from=AFNI /opt/afni/ /opt/afni/
