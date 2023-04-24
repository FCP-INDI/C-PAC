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
FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114 as FSL
USER root

ARG DEBIAN_FRONTEND=noninteractive

# install FSL
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      fsl-core=5.0.9-5~nd16.04+1 \
      fsl-atlases \
      fsl-mni152-templates=5.0.7-2

# setup FSL environment
ENV FSLDIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    AFNI_MODELPATH="/usr/lib/afni/models" \
    AFNI_IMSAVE_WARNINGS="NO" \
    AFNI_TTATLAS_DATASET="/usr/share/afni/atlases" \
    AFNI_PLUGINPATH="/usr/lib/afni/plugins" \
    PATH=/usr/lib/fsl/5.0:$PATH

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FSL 5.0.9-5 stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=FSL /etc/fsl /etc/fsl
COPY --from=FSL /usr/lib/fsl /usr/lib/fsl
COPY --from=FSL /usr/lib/libnewmat.so.10 /usr/lib/libnewmat.so.10
COPY --from=FSL /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
COPY --from=FSL /usr/lib/libniftiio.so.2 /usr/lib/libniftiio.so.2
COPY --from=FSL /usr/lib/libznz.so.2 /usr/lib/libznz.so.2
COPY --from=FSL /usr/share/doc/fsl-core /usr/share/doc/fsl-core
COPY --from=FSL /usr/share/man/man1/fsl-5.0-core.1.gz /usr/share/man/man1/
COPY --from=FSL /usr/share/man/man1/fsl.1.gz /usr/share/man/man1/
COPY --from=FSL /usr/share/data/fsl-mni152-templates /usr/share/data/fsl-mni152-templates
COPY --from=FSL /usr/share/doc/fsl-mni152-templates /usr/share/doc/fsl-mni152-templates
COPY --from=FSL /usr/share/fsl /usr/share/fsl
# install C-PAC resources into FSL
COPY --from=data /fsl_data/standard /usr/share/fsl/5.0/data/standard
