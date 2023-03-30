FROM ghcr.io/fcp-indi/c-pac/ubuntu:python3.10-bionic-non-free as FreeSurfer

USER root
RUN curl -sL https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.3.2/freesurfer_ubuntu18-7.3.2_amd64.deb -o /tmp/freesurfer_ubuntu18-7.3.2_amd64.deb \
  && dpkg -i /tmp/freesurfer_ubuntu18-7.3.2_amd64.deb

ENV 
    NO_FSFAST=1
RUN mkdir -p /opt/freesurfer \
  && curl -fsSL --retry 5 https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.3.2/freesurfer-linux-ubuntu18_amd64-7.3.2.tar.gz \
  | tar -xz -C /opt/freesurfer --strip-components 1

# Only keep what we need
FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FreeSurfer 7.3.2 stage"