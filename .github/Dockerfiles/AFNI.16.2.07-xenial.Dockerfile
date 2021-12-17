FROM ghcr.io/shnizzedy/c-pac/ubuntu:xenial-20200114

USER root

# set up AFNI
ENV PATH=/opt/afni:$PATH
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    afni=16.2.07~dfsg.1-5~nd16.04+1 && \
    ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENTRYPOINT ["/bin/bash"]

# set user
USER c-pac_user