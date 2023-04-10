# using neurodebian runtime as parent image
FROM neurodebian:bionic-non-free as MSM

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      curl \
      libexpat1-dev \
    && apt-get autoremove -y \
    && apt-get autoclean -y

#---------------------
# Install MSM Binaries
#---------------------
COPY dev/docker_data/checksum/msm.2.0.md5 /tmp/checksum.md5
RUN mkdir /opt/msm \
    && curl -ksSL --retry 5 https://www.doc.ic.ac.uk/~ecr05/MSM_HOCR_v2/MSM_HOCR_v2-download.tgz -o msm.tgz \
    && md5sum --check /tmp/checksum.md5 \
    && tar zx -C /opt -f msm.tgz \
    && mv /opt/homes/ecr05/MSM_HOCR_v2/* /opt/msm/ \
    && rm -rf /opt/homes /opt/msm/MacOSX /opt/msm/Centos \
    && chmod +x /opt/msm/Ubuntu/* \
    && MSMBINDIR=/opt/msm/Ubuntu \
       PATH=$PATH:/opt/msm/Ubuntu

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig \
    && apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
msm v2.0 stage \
    Multimodal Surface Matching with Higher order Clique Reduction Version 2.00 (Feb 2017)"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=MSM /opt/msm/Ubuntu/msm /opt/msm/Ubuntu/msm