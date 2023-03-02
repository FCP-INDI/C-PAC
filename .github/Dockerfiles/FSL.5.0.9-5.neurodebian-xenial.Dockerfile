FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114 as FSL
USER root

ARG DEBIAN_FRONTEND=noninteractive

# install FSL
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      fsl-core=5.0.9-5~nd16.04+1 \
      fsl-atlases \
      fsl-mni152-templates=5.0.7-2

# setup FSL environment
ENV FSLDIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    AFNI_MODELPATH="/usr/lib/afni/models" \
    AFNI_IMSAVE_WARNINGS="NO" \
    AFNI_TTATLAS_DATASET="/usr/share/afni/atlases" \
    AFNI_PLUGINPATH="/usr/lib/afni/plugins" \
    PATH=/usr/lib/fsl/5.0:$PATH

# install CPAC resources into FSL
RUN curl -sL http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz -o /tmp/cpac_resources.tar.gz && \
    tar xfz /tmp/cpac_resources.tar.gz -C /tmp && \
    cp -n /tmp/cpac_image_resources/MNI_3mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/MNI_4mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/symmetric/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/2mm $FSLDIR/data/standard/tissuepriors && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/3mm $FSLDIR/data/standard/tissuepriors && \
    chmod -R ugo+r $FSLDIR/data/standard

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FSL 5.0.9-5 stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=FSL /etc/fsl /etc/fsl
COPY --from=FSL /usr/lib/fsl /usr/lib/fsl
COPY --from=FSL /usr/lib/libnewmat.so.10 /usr/lib/libnewmat.so.10
COPY --from=FSL /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
COPY --from=FSL /usr/lib/libniftiio.so.2 /usr/lib/libniftiio.so.2
COPY --from=FSL /usr/lib/libznz.so.2 /usr/lib/libznz.so.2
COPY --from=FSL /usr/share/doc/fsl-core /usr/share/doc/fsl-core
COPY --from=FSL /usr/share/man/man1/fsl-5.0-core.1.gz /usr/share/man/man1/
COPY --from=FSL /usr/share/man/man1/fsl.1.gz /usr/share/man/man1/
COPY --from=FSL /usr/share/data/fsl-mni152-templates /usr/share/data/fsl-mni152-templates
COPY --from=FSL /usr/share/doc/fsl-mni152-templates /usr/share/doc/fsl-mni152-templates
COPY --from=FSL /usr/share/fsl /usr/share/fsl
