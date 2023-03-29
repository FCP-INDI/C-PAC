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
FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free AS FSL

USER root

# set up FSL environment
ENV FSLDIR=/usr/share/fsl/6.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/6.0 \
    LD_LIBRARY_PATH=/usr/share/fsl/6.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/share/fsl/6.0/bin:$PATH \
    TZ=America/New_York

# Installing and setting up FSL
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    &&  echo $TZ > /etc/timezone \
    && apt-get update \
    && apt-get install -y tclsh wish \
    && echo "Downloading FSL ..." \
    && mkdir -p /usr/share/fsl/6.0 \
    && curl -sSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.4-centos6_64.tar.gz \
    | tar zx -C /usr/share/fsl/6.0 --strip-components=1 \
    --exclude=fsl/bin/mist \
    --exclude=fsl/bin/possum \
    --exclude=fsl/data/possum \
    --exclude=fsl/data/mist \
    --exclude=fsl/data/first \
    && ln -s /usr/share/fsl/6.0 /usr/share/fsl/5.0 \
    && if [ -f /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1.2.0]; then \
        ln -svf /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1.2.0 /usr/lib/x86_64-linux-gnu/libGL.so.1; \
    fi \
    && ldconfig \
    && curl -sL http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz -o /tmp/cpac_resources.tar.gz \
    && tar xfz /tmp/cpac_resources.tar.gz -C /tmp \
    && cp -n /tmp/cpac_image_resources/MNI_3mm/* $FSLDIR/data/standard \
    && cp -n /tmp/cpac_image_resources/MNI_4mm/* $FSLDIR/data/standard \
    && cp -n /tmp/cpac_image_resources/symmetric/* $FSLDIR/data/standard \
    && cp -n /tmp/cpac_image_resources/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford \
    && cp -nr /tmp/cpac_image_resources/tissuepriors/2mm $FSLDIR/data/standard/tissuepriors \
    && cp -nr /tmp/cpac_image_resources/tissuepriors/3mm $FSLDIR/data/standard/tissuepriors \
    && chmod -R ugo+r $FSLDIR/data/standard

ENTRYPOINT ["/bin/bash"]

# set user
USER c-pac_user

# Only keep what we need
FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FSL 6.0.4 stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=FSL /usr/bin/tclsh /usr/bin/tclsh
COPY --from=FSL /usr/bin/wish /usr/bin/wish
COPY --from=FSL /usr/share/fsl /usr/share/fsl
COPY --from=FSL /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
# install C-PAC resources into FSL
COPY --from=data /fsl_data/standard /usr/share/fsl/6.0/data/standard
