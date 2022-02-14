FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free as base

USER root

RUN apt-get update && \
    apt-get install -y connectome-workbench=1.5.0-1~nd18.04+1 && \
    ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER c-pac_user

FROM scratch

COPY --from=base /usr/bin/wb_* /usr/bin/
COPY --from=base /usr/share/applications/connectome-workbench.desktop /usr/share/applications/connectome-workbench.desktop
COPY --from=base /usr/share/bash-completion/completions/wb* /usr/share/bash_completion/completions/
COPY --from=base /usr/share/doc/connectome-workbench /usr/share/doc/connectome-workbench
# manfile is missing from this version of wb_command
# COPY --from=base /usr/share/man/man1/wb_* /usr/share/man/man1/   
COPY --from=base /usr/share/pixmaps/connectome-workbench.png /usr/share/pixmaps/connectome-workbench.png
