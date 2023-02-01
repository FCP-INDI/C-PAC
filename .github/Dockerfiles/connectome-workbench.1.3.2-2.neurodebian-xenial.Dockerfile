FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114 as base

USER root

# install wb_command v1.3.2
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    connectome-workbench=1.3.2-2~nd16.04+1 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user

# only keep what we need
FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
connectome-workbench 1.3.2-2 stage"
COPY --from=base /lib64/ld-linux-x86-64.so.2 /lib64/
COPY --from=base /usr/bin/wb_* /usr/bin/
COPY --from=base /lib/x86_64-linux-gnu/ld-* /lib/x86_64-linux-gnu/
COPY --from=base /usr/lib/x86_64-linux-gnu/lib*.so* /usr/lib/x86_64-linux-gnu/
COPY --from=base /usr/share/applications/connectome-workbench.desktop /usr/share/applications/connectome-workbench.desktop
COPY --from=base /usr/share/bash-completion/completions/wb* /usr/share/bash-completion/completions/
COPY --from=base /usr/share/doc/connectome-workbench /usr/share/doc/connectome-workbench
COPY --from=base /usr/share/man/man1/wb_* /usr/share/man/man1/
COPY --from=base /usr/share/pixmaps/connectome-workbench.png /usr/share/pixmaps/connectome-workbench.png