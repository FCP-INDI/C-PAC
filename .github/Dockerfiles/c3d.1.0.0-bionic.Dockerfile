FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free as c3d
USER root

# Installing and setting up c3d
RUN mkdir -p /opt/c3d && \
    curl -sSL "http://downloads.sourceforge.net/project/c3d/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz" \
    | tar -xzC /opt/c3d --strip-components 1
ENV C3DPATH /opt/c3d/
ENV PATH $C3DPATH/bin:$PATH

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
c3d 1.0.0 (Bionic) stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=c3d /opt/c3d/ /opt/c3d/
