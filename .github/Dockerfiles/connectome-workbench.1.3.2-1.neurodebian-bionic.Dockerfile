FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free as base

USER root

# install wb_command v1.3.2
RUN APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 04EE7237B7D453EC && \
    APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 648ACFD622F3D138 && \
    APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key adv --keyserver keyserver.ubuntu.com --recv-keys DCC9EFBF77E11517 && \
    printf '\ndeb http://httpredir.debian.org/debian/ buster main non-free' >> /etc/apt/sources.list && \
    apt-get update && \
    apt-get install connectome-workbench=1.3.2-1 -y && \
    strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

# set user
USER c-pac_user

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD: connectome-workbench 1.3.2-1 stage"
COPY --from=base /usr/bin/wb_* /usr/bin/
COPY --from=base /usr/share/applications/connectome-workbench.desktop /usr/share/applications/connectome-workbench.desktop
COPY --from=base /usr/share/bash-completion/completions/wb* /usr/share/bash_completion/completions/
COPY --from=base /usr/share/doc/connectome-workbench /usr/share/doc/connectome-workbench
COPY --from=base /usr/share/man/man1/wb_* /usr/share/man/man1/
COPY --from=base /usr/share/pixmaps/connectome-workbench.png /usr/share/pixmaps/connectome-workbench.png