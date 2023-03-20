FROM ghcr.io/shnizzedy/c-pac/fsl:neurodebian-bionic as FSL-Neurodebian
FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free AS FSL

USER root

# set up FSL environment
ENV FSLDIR=/usr/share/fsl/6.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/6.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/6.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/lib/fsl/6.0:$PATH \
    TZ=America/New_York


# Installing and setting up FSL
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
     echo $TZ > /etc/timezone && \
    apt-get update && \
    apt-get install -y tclsh wish && \
    echo "Downloading FSL ..." \
    && mkdir -p /usr/share/fsl/6.0 \
    && curl -sSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.4-centos6_64.tar.gz \
    | tar zx -C /usr/share/fsl/6.0 --strip-components=1 \
    --exclude=fsl/bin/mist \
    --exclude=fsl/bin/possum \
    --exclude=fsl/data/possum \
    --exclude=fsl/data/mist \
    --exclude=fsl/data/first \
    && ln -s /usr/share/fsl/6.0 /usr/share/fsl/5.0

ENTRYPOINT ["/bin/bash"]

# set user
USER c-pac_user

# Only keep what we need
FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FSL 6.0.4 stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=FSL /usr/bin/tclsh /usr/bin/tclsh
COPY --from=FSL /usr/bin/wish /usr/bin/wish
COPY --from=FSL /usr/share/fsl /usr/share/fsl
COPY --from=FSL /usr/lib /usr/lib
COPY --from=FSL /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
COPY --from=FSL-Neurodebian /usr/share/fsl/5.0/data/standard/tissuepriors/*mm /usr/share/fsl/6.0/data/standard/tissuepriors/