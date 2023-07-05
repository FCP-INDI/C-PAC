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
FROM ghcr.io/fcp-indi/c-pac/ubuntu:jammy-non-free as AFNI
USER root

# install AFNI
COPY dev/docker_data/required_afni_pkgs.txt /opt/required_afni_pkgs.txt
COPY dev/docker_data/checksum/AFNI.23.1.10.sha384 /tmp/AFNI.23.1.10.sha384
RUN apt-get update \
    && apt-get install -y \
      build-essential \
      cmake \
      curl \
      eog \
      evince \
      firefox \
      gedit \
      gnome-terminal \
      gnome-tweaks \
      gsl-bin \
      libcurl4-openssl-dev \
      libgdal-dev \
      libgfortran-11-dev \
      libglu1-mesa-dev \
      libglw1-mesa-dev \
      libgomp1 \
      libjpeg62 \
      libnode-dev \
      libopenblas-dev \
      libssl-dev \
      libudunits2-dev \
      libxm4 \
      libxml2-dev \
      nautilus \
      netpbm \
      python-is-python3 \
      python3-matplotlib \
      python3-numpy \
      python3-pil \
      r-base-dev \
      tcsh \
      vim \
      xfonts-100dpi \
      xfonts-base \
      xterm \
      xvfb \
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
    && apt-get remove -y libglw1-mesa-dev \
    && ldconfig

# set up AFNI
ENV PATH=/opt/afni:$PATH

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
AFNI 23.1.10 (Publius Helvius Pertinax) stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=AFNI /lib/x86_64-linux-gnu/ld* /lib/x86_64-linux-gnu/
COPY --from=AFNI /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
COPY --from=AFNI /lib64/ld* /lib64/
COPY --from=AFNI /opt/afni/ /opt/afni/
