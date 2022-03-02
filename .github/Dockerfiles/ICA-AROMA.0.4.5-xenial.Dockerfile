LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD: ICA-AROMA 0.4.5 stage"
FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114 AS ICA-AROMA

USER root

# Installing and setting up ICA_AROMA
RUN mkdir -p /opt/ICA-AROMA && \
  curl -sSL "https://github.com/oesteban/ICA-AROMA/archive/v0.4.5.tar.gz" \
  | tar -xzC /opt/ICA-AROMA --strip-components 1 && \
  chmod +x /opt/ICA-AROMA/ICA_AROMA.py
ENV PATH="/opt/ICA-AROMA:$PATH" \
    AROMA_VERSION="0.4.5"

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user

# Only keep what we need
FROM scratch
COPY --from=ICA-AROMA /opt/ICA-AROMA/ /opt/ICA-AROMA/