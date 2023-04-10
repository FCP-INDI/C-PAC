FROM ghcr.io/fcp-indi/c-pac/ubuntu:python3.10-bionic-non-free as ANTs

USER root
COPY dev/docker_data/checksum/ANTs.2.4.3.md5 /tmp/checksum.md5
RUN curl -sL https://github.com/ANTsX/ANTs/releases/download/v2.4.3/ants-2.4.3-ubuntu-18.04-X64-gcc.zip -o /tmp/ANTs.zip \
  && curl -sL https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/3133832/Oasis.zip -o /tmp/Oasis.zip \
  && md5sum --check /tmp/checksum.md5 \
  && unzip /tmp/ANTs.zip -d /tmp \
  && mkdir /usr/lib/ants \
  && mv /tmp/ants-2.4.3/* /usr/lib/ants \
  && mkdir /ants_template \
  && unzip /tmp/Oasis.zip -d /tmp \
  && mv /tmp/MICCAI2012-Multi-Atlas-Challenge-Data /ants_template/oasis

# Only keep what we need
FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
ANTs 2.4.3 stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=ANTs /usr/lib/ants/ /usr/lib/ants/
COPY --from=ANTs /ants_template/ /ants_template/
