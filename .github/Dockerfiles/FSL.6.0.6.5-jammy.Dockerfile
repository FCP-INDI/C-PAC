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
FROM ghcr.io/fcp-indi/c-pac/fsl:data as data
FROM ghcr.io/fcp-indi/c-pac/ubuntu:jammy-non-free AS FSL

USER root
COPY ./dev/docker_data/fsl/6.0.6.5 /tmp/fsl/manifest

# set up FSL environment
ENV FSLDIR=/usr/share/fsl/6.0
ENV FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    FSLTCLSH=$FSLDIR/bin/fsltclsh \
    FSLWISH=$FSLDIR/bin/fslwish \
    LD_LIBRARY_PATH=/usr/share/fsl/6.0:$LD_LIBRARY_PATH \
    PATH=/usr/share/fsl/6.0/bin:$PATH \
    TZ=America/New_York \
    USER=c-pac_user

# Installing and setting up FSL
# Full list of packages, included excluded packages, in dev/docker_data/fsl/6.0.6.5/fsl-release.yml
RUN apt-get update > /dev/null && \
    apt-get install --no-install-recommends --yes \
        wget bzip2 ca-certificates \
        git \
        tini \
        > /dev/null && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/23.1.0-4/Mambaforge-23.1.0-4-Linux-x86_64.sh -O /tmp/mambaforge.sh && \
    wget --no-hsts --quiet https://git.fmrib.ox.ac.uk/fsl/conda/installer/-/raw/3.5.3/fsl/installer/fslinstaller.py?inline=false -O /tmp/fslinstaller.py \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone \
    && mkdir -p /usr/share/fsl \
    && yes | python3 /tmp/fslinstaller.py \
      -V 6.0.6.5 \
      -d ${FSLDIR} \
      --exclude_package '*asl' \
      --exclude_package '*basil*' \
      --exclude_package '*baycest*' \
      --exclude_package '*bianca*' \
      --exclude_package '*cifti*' \
      --exclude_package '*cud*' \
      --exclude_package '*deface*' \
      --exclude_package '*eddy*' \
      --exclude_package '*fabber*' \
      --exclude_package '*fdt*' \
      --exclude_package '*feed*' \
      --exclude_package '*first*' \
      --exclude_package '*free*' \
      --exclude_package '*gifti*' \
      --exclude_package '*lesions*' \
      --exclude_package '*mist*' \
      --exclude_package '*mo*' \
      --exclude_package '*mrs*' \
      --exclude_package '*pack*' \
      --exclude_package '*possum*' \
      --exclude_package '*ptx*' \
      --exclude_package '*shapemodel*' \
      --exclude_package '*siena*' \
      --exclude_package '*surface*' \
      --exclude_package '*susan*' \
      --exclude_package '*swe*' \
      --exclude_package '*tbss*' \
      --exclude_package '*tirl*' \
      --exclude_package '*verbena*' \
      --exclude_package '*wxpython*' \
      --exclude_package '*xt*' \
      --manifest '/tmp/fsl/manifest/manifest.json' \
      --miniconda '/tmp/mambaforge.sh' \
    && cd $FSLDIR \
    && conda remove gcc make \
    && ldconfig

ENTRYPOINT ["/bin/bash"]

# # Only keep what we need
# FROM scratch
# LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
# FSL 6.0.6.5 stage"
# LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
# COPY --from=FSL /lib/x86_64-linux-gnu /lib/x86_64-linux-gnu
# COPY --from=FSL /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu
# COPY --from=FSL /usr/bin /usr/bin
# COPY --from=FSL /usr/local/bin /usr/local/bin
# COPY --from=FSL /usr/share/fsl /usr/share/fsl
# install C-PAC resources into FSL
COPY --from=data /fsl_data/standard /usr/share/fsl/6.0/data/standard
COPY --from=data /fsl_data/atlases /usr/share/fsl/6.0/data/atlases
