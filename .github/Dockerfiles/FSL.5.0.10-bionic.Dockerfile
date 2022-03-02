LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD: FSL 5.0.10 stage"
FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free AS FSL

USER root

# set up FSL environment
ENV FSLDIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/lib/fsl/5.0:$PATH \
    TZ=America/New_York


# Installing and setting up FSL
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
     echo $TZ > /etc/timezone && \
    apt-get update && \
    apt-get install -y tclsh wish && \
    echo "Downloading FSL ..." \
    && mkdir -p /usr/share/fsl/5.0 \
    && curl -sSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-5.0.10-centos6_64.tar.gz \
    | tar zx -C /usr/share/fsl/5.0 --strip-components=1 \
    --exclude=fsl/bin/mist \
    --exclude=fsl/bin/possum \
    --exclude=fsl/data/possum \
    --exclude=fsl/data/mist \
    --exclude=fsl/data/first

ENTRYPOINT ["/bin/bash"]

# set user
USER c-pac_user

# Only keep what we need
FROM scratch
COPY --from=FSL /usr/bin/tclsh /usr/bin/tclsh
COPY --from=FSL /usr/bin/wish /usr/bin/wish
COPY --from=FSL /usr/share/fsl/ /usr/share/fsl/
COPY --from=FSL /usr/lib/ /usr/lib/