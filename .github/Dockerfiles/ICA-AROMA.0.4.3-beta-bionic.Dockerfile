FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free AS ICA-AROMA
USER root

# install ICA-AROMA
RUN mkdir -p /opt/ICA-AROMA
RUN curl -sL https://github.com/rhr-pruim/ICA-AROMA/archive/v0.4.3-beta.tar.gz | tar -xzC /opt/ICA-AROMA --strip-components 1
RUN chmod +x /opt/ICA-AROMA/ICA_AROMA.py
ENV PATH=/opt/ICA-AROMA:$PATH

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
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
ICA-AROMA 0.4.3-beta stage"
COPY --from=ICA-AROMA /opt/ICA-AROMA/ /opt/ICA-AROMA/