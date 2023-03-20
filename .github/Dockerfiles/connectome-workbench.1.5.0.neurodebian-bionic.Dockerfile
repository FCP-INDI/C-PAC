FROM ghcr.io/fcp-indi/c-pac/ubuntu:python3.10-bionic-non-free as base

USER root

RUN curl -sSL "https://www.humanconnectome.org/storage/app/media/workbench/workbench-linux64-v1.5.0.zip" -o /opt/workbench.zip \
  && unzip /opt/workbench.zip -d /opt \
  && rm -rf /opt/workbench.zip
ENV PATH $PATH:/opt/workbench/bin_linux64

USER c-pac_user

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
connectome-workbench 1.5.0 stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=base /opt/workbench /opt/workbench
