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
FROM ghcr.io/fcp-indi/c-pac/afni:23.3.09-jammy as AFNI
FROM ghcr.io/fcp-indi/c-pac/ants:2.4.3-jammy as ANTs
FROM ghcr.io/fcp-indi/c-pac/c3d:1.0.0-jammy as c3d
FROM ghcr.io/fcp-indi/c-pac/connectome-workbench:1.5.0.neurodebian-jammy as connectome-workbench
FROM ghcr.io/fcp-indi/c-pac/fsl:6.0.6.5-jammy as FSL
FROM ghcr.io/fcp-indi/c-pac/ica-aroma:0.4.4-beta-jammy as ICA-AROMA

FROM ghcr.io/fcp-indi/c-pac/ubuntu:jammy-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Standard software dependencies for C-PAC standard and lite images"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# Installing connectome-workbench
COPY --from=connectome-workbench /lib/x86_64-linux-gnu /lib/x86_64-linux-gnu
COPY --from=connectome-workbench /lib64/ld-linux-x86-64.so.2 /lib64/ld-linux-x86-64.so.2
COPY --from=connectome-workbench /usr/bin/wb_* /usr/bin
COPY --from=connectome-workbench /usr/share/bash-completion/completions/wb_command /usr/share/bash-completion/completions/wb_command
COPY --from=connectome-workbench /usr/share/doc/connectome-workbench /usr/share/doc/connectome-workbench
COPY --from=connectome-workbench /usr/share/man/man1/wb_* /usr/share/man/man1
COPY --from=connectome-workbench /usr/share/bash-completion/completions/wb_shortcuts /usr/share/bash-completion/completions/wb_shortcuts

# Installing FSL
ENV FSLDIR=/usr/share/fsl/6.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    TZ=America/New_York
ENV FSLTCLSH=$FSLDIR/bin/fsltclsh \
    FSLWISH=$FSLDIR/bin/fslwish \
    LD_LIBRARY_PATH=${FSLDIR}/6.0:$LD_LIBRARY_PATH \
    PATH=${FSLDIR}/bin:$PATH \
    TZ=America/New_York \
    USER=c-pac_user
COPY --from=FSL /lib/x86_64-linux-gnu /lib/x86_64-linux-gnu
COPY --from=FSL /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu
COPY --from=FSL /usr/bin /usr/bin
COPY --from=FSL /usr/local/bin /usr/local/bin
COPY --from=FSL /usr/share/fsl /usr/share/fsl

# Installing C-PAC dependencies
COPY requirements.txt /opt/requirements.txt
RUN mamba install git -y \
  && pip install -r /opt/requirements.txt \
  && rm -rf /opt/requirements.txt \
  && yes | mamba clean --all \
  && rm -rf /usr/share/fsl/6.0/pkgs/cache/*

# Installing and setting up c3d
COPY --from=c3d /opt/c3d/ opt/c3d/
ENV C3DPATH /opt/c3d
ENV PATH $C3DPATH/bin:$PATH

# Installing AFNI
COPY --from=AFNI /lib/x86_64-linux-gnu/ld* /lib/x86_64-linux-gnu/
COPY --from=AFNI /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
COPY --from=AFNI /lib64/ld* /lib64/
COPY --from=AFNI /opt/afni/ /opt/afni/
# set up AFNI
ENV PATH=/opt/afni:$PATH

# Installing ANTs
ENV PATH=/usr/lib/ants/bin:$PATH
COPY --from=ANTs /usr/lib/ants/ /usr/lib/ants/
COPY --from=ANTs /ants_template/ /ants_template/

# Installing ICA-AROMA
COPY --from=ICA-AROMA /opt/ICA-AROMA/ /opt/ICA-AROMA/
ENV PATH=/opt/ICA-AROMA:$PATH

# link libraries & clean up
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/* \
    && find / -type f -print0 | sort -t/ -k2 | xargs -0 rdfind -makehardlinks true \
    && rm -rf results.txt \
    && ldconfig \
    && chmod 777 / /home/c-pac_user \
    && chmod 777 $(ls / | grep -v sys | grep -v proc)

# set user
USER c-pac_user
