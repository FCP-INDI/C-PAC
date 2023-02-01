# we need mri_vol2vol which is not included in neurodocker freesurfer 6.0.0-min
FROM freesurfer/freesurfer:6.0 as FreeSurfer

FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FreeSurfer 6.0.0-min stage"
USER root

# install FreeSurfer
# set shell to BASH
RUN mkdir -p /usr/lib/freesurfer
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH"
SHELL ["/bin/bash", "-c"]
RUN curl -fsSL --retry 5 https://dl.dropbox.com/s/nnzcfttc41qvt31/recon-all-freesurfer6-3.min.tgz \
    | tar -xz -C /usr/lib/freesurfer --strip-components 1 && \
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
RUN printf 'source $FREESURFER_HOME/SetUpFreeSurfer.sh' > ~/.bashrc
# restore shell to default (sh)
SHELL ["/bin/sh", "-c"]
COPY dev/docker_data/license.txt $FREESURFER_HOME/license.txt
COPY dev/docker_data/FreeSurfer_bash /bin/

COPY --from=FreeSurfer opt/freesurfer/bin/mri_vol2vol /usr/lib/freesurfer/bin/mri_vol2vol
COPY --from=FreeSurfer opt/freesurfer/bin/mri_vol2vol.bin /usr/lib/freesurfer/bin/mri_vol2vol.bin

ENTRYPOINT ["/bin/FreeSurfer_bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user