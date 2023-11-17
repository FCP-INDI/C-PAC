# Copyright (C) 2021-2023  C-PAC Developers

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
# using neurodebian runtime as parent image
FROM neurodebian:bionic-non-free as MSM
USER root
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      curl \
      libexpat1-dev \
    && apt-get autoremove -y \
    && apt-get autoclean -y

#---------------------
# Install MSM Binaries
#---------------------
COPY dev/docker_data/checksum/msm.2.0.sha384 /tmp/checksum.sha384
RUN mkdir /opt/msm \
    && curl -ksSL --retry 5 https://www.doc.ic.ac.uk/~ecr05/MSM_HOCR_v2/MSM_HOCR_v2-download.tgz -o msm.tgz \
    && sha384sum --check /tmp/checksum.sha384 \
    && tar zx -C /opt -f msm.tgz \
    && mv /opt/homes/ecr05/MSM_HOCR_v2/* /opt/msm/ \
    && rm -rf /opt/homes /opt/msm/MacOSX /opt/msm/Centos \
    && chmod +x /opt/msm/Ubuntu/* \
    && MSMBINDIR=/opt/msm/Ubuntu \
       PATH=$PATH:/opt/msm/Ubuntu

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig \
    && apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /var/tmp/* \
    && /bin/bash -O extglob -c 'rm -rfv /tmp/!("home/c-pac_user")'

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
msm v2.0 stage \
    Multimodal Surface Matching with Higher order Clique Reduction Version 2.00 (Feb 2017)"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=MSM /opt/msm/Ubuntu/msm /opt/msm/Ubuntu/msm