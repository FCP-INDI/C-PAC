FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114

USER root

# Installing and setting up c3d
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    convert3d && \
    ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENTRYPOINT ["/bin/bash"]

# set user
USER c-pac_user