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
FROM ghcr.io/fcp-indi/c-pac/freesurfer:6.0.0-min.neurodocker-jammy as FreeSurfer

FROM ghcr.io/fcp-indi/c-pac/stage-base:lite-v1.8.7dev1
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Standard software dependencies for C-PAC standard images"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# Installing FreeSurfer
RUN apt-get update \
    && apt-get install --no-install-recommends -y bc \
    && yes | mamba install tcsh \
    && yes | mamba clean --all \
    && cp -l `which tcsh` /bin/tcsh \
    && cp -l `which tcsh` /bin/csh
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    NO_FSFAST=1
ENV PATH="$FREESURFER_HOME/bin:$PATH" \
    PERL5LIB="$FREESURFER_HOME/mni/share/perl5" \
    FSFAST_HOME="$FREESURFER_HOME/fsfast" \
    SUBJECTS_DIR="$FREESURFER_HOME/subjects" \
    MNI_DIR="$FREESURFER_HOME/mni"
ENV MINC_BIN_DIR="$MNI_DIR/bin" \
    MINC_LIB_DIR="$MNI_DIR/lib" \
    PATH="$PATH:$MINC_BIN_DIR"
COPY --from=FreeSurfer /usr/lib/freesurfer/ /usr/lib/freesurfer/
COPY dev/docker_data/license.txt $FREESURFER_HOME/license.txt

# link libraries & clean up
RUN apt-get autoremove -y \
    && apt-get autoclean -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/* \
    && find / -type f -print0 | sort -t/ -k2 | xargs -0 rdfind -makehardlinks true \
    && rm -rf results.txt \
    && ldconfig \
    && chmod 777 / /home/c-pac_user \
    && chmod 777 $(ls / | grep -v sys | grep -v proc)

# set user
USER c-pac_user
