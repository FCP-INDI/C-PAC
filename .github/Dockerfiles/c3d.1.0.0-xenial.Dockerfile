FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114 AS c3d

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

# Only keep what we need
FROM scratch
COPY --from=c3d /usr/bin/c3d/ /usr/bin/c3d/
COPY --from=c3d /usr/lib/c3d* /usr/lib/