FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
c3d 1.0.0 (Xenial) stage"
USER root

# Installing and setting up c3d
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    convert3d=0.0.20190204-1~nd16.04+1
ENV C3DPATH /usr/bin/
ENV PATH $C3DPATH/bin:$PATH

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user