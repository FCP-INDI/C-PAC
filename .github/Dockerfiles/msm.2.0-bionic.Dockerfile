# using neurodebian runtime as parent image
FROM neurodebian:bionic-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD: msm v2.0 stage \
    Multimodal Surface Matching with Higher order Clique Reduction Version 2.00 (Feb 2017)"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      curl \
      libexpat1-dev && \
    apt-get autoremove -y && \
    apt-get autoclean -y

#---------------------
# Install MSM Binaries
#---------------------
RUN mkdir /opt/msm
RUN curl -ksSL --retry 5 https://www.doc.ic.ac.uk/~ecr05/MSM_HOCR_v2/MSM_HOCR_v2-download.tgz | tar zx -C /opt
RUN mv /opt/homes/ecr05/MSM_HOCR_v2/* /opt/msm/
RUN rm -rf /opt/homes /opt/msm/MacOSX /opt/msm/Centos
RUN chmod +x /opt/msm/Ubuntu/*
ENV MSMBINDIR=/opt/msm/Ubuntu \
    PATH=$PATH:/opt/msm/Ubuntu

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
