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
FROM ghcr.io/fcp-indi/c-pac/fsl:data as data
FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free AS FSL

USER root

# set up FSL environment
ENV FSLDIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/lib/fsl/5.0:$PATH \
    TZ=America/New_York

# Installing and setting up FSL
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
     echo $TZ > /etc/timezone && \
    apt-get update && \
    apt-get install -y tclsh wish && \
    echo "Downloading FSL ..." \
    && mkdir -p /usr/share/fsl/5.0 \
    && curl -sSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-5.0.10-centos6_64.tar.gz \
    | tar zx -C /usr/share/fsl/5.0 --strip-components=1 \
    --exclude=fsl/bin/mist \
    --exclude=fsl/bin/possum \
    --exclude=fsl/data/possum \
    --exclude=fsl/data/mist \
    --exclude=fsl/data/first

ENTRYPOINT ["/bin/bash"]

# set user
USER c-pac_user

# Only keep what we need
FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FSL 5.0.10 stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=FSL /usr/bin/tclsh /usr/bin/tclsh
COPY --from=FSL /usr/bin/wish /usr/bin/wish
COPY --from=FSL /usr/share/fsl/ /usr/share/fsl/
COPY --from=FSL /usr/lib/ /usr/lib/
COPY --from=FSL /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
# install C-PAC resources into FSL
COPY --from=data /fsl_data/standard /usr/share/fsl/5.0/data/standard
