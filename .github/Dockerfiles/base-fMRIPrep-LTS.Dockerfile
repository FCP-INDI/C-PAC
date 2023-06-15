# Choose versions
FROM ghcr.io/fcp-indi/c-pac/afni:16.2.07.neurodocker-xenial as AFNI
FROM ghcr.io/fcp-indi/c-pac/ants:2.3.4.neurodocker-xenial as ANTs
FROM ghcr.io/fcp-indi/c-pac/c3d:1.0.0-xenial as c3d
FROM ghcr.io/fcp-indi/c-pac/connectome-workbench:1.3.2-2.neurodebian-xenial as connectome-workbench
FROM ghcr.io/fcp-indi/c-pac/freesurfer:6.0.1-min-xenial as FreeSurfer
FROM ghcr.io/fcp-indi/c-pac/fsl:5.0.9-5.neurodebian-xenial as FSL
FROM ghcr.io/fcp-indi/c-pac/ica-aroma:0.4.5-xenial as ICA-AROMA
FROM ghcr.io/fcp-indi/c-pac/msm:2.0-bionic as MSM

# using Ubuntu 16.04 LTS as parent image
FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Software dependencies version-matched to `fMRIPrep LTS <https://reproducibility.stanford.edu/fmriprep-lts#long-term-support-lts->`_"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# allow users to update / create themselves
RUN chmod ugo+w /etc/passwd

# set up FSL environment
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

# install ANTs
ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ANTSPATH="/usr/lib/ants/" \
    PATH="/usr/lib/ants:$PATH"
COPY --from=ANTs /usr/lib/ants/ /usr/lib/ants/
COPY --from=ANTs /ants_template /ants_template

# install ICA-AROMA
COPY --from=ICA-AROMA /opt/ICA-AROMA/ /opt/ICA-AROMA/
ENV PATH="/opt/ICA-AROMA:$PATH" \
    AROMA_VERSION="0.4.5"

# install FreeSurfer
COPY --from=FreeSurfer /opt/freesurfer/ /usr/lib/freesurfer/
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH" \
    NO_FSFAST=1 \
    FSL_DIR="$FSLDIR" \
    OS="Linux" \
    FS_OVERRIDE=0 \
    FIX_VERTEX_AREA="" \
    FSF_OUTPUT_FORMAT="nii.gz"
ENV SUBJECTS_DIR="$FREESURFER_HOME/subjects" \
    FUNCTIONALS_DIR="$FREESURFER_HOME/sessions" \
    MNI_DIR="$FREESURFER_HOME/mni" \
    LOCAL_DIR="$FREESURFER_HOME/local" \
    MINC_BIN_DIR="$FREESURFER_HOME/mni/bin" \
    MINC_LIB_DIR="$FREESURFER_HOME/mni/lib" \
    MNI_DATAPATH="$FREESURFER_HOME/mni/data"
ENV PERL5LIB="$MINC_LIB_DIR/perl5/5.8.5" \
    MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8.5" \
    PATH="$FREESURFER_HOME/bin:$FSFAST_HOME/bin:$FREESURFER_HOME/tktools:$MINC_BIN_DIR:$PATH"

# Installing and setting up AFNI
COPY --from=AFNI /lib/x86_64-linux-gnu/ld* /lib/x86_64-linux-gnu/
COPY --from=AFNI /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/ 
COPY --from=AFNI /lib64/ld* /lib64/
COPY --from=AFNI /usr/lib/afni/bin/ /usr/lib/afni/bin/
COPY --from=AFNI /usr/lib/x86_64-linux-gnu/lib*so* /usr/lib/x86_64-linux-gnu/
# set up AFNI
ENV PATH=/usr/lib/afni/bin:$PATH \
    AFNI_MODELPATH="/usr/lib/afni/models" \
    AFNI_IMSAVE_WARNINGS="NO" \
    AFNI_TTATLAS_DATASET="/usr/share/afni/atlases" \
    AFNI_PLUGINPATH="/usr/lib/afni/plugins"

# Intalling and setting up c3d
COPY --from=c3d /usr/lib/libITK* /usr/lib/
COPY --from=c3d /usr/lib/libitk* /usr/lib/
COPY --from=c3d /usr/lib/x86_64-linux-gnu/libgdcm* /usr/lib/x86_64-linux-gnu/
COPY --from=c3d /usr/lib/x86_64-linux-gnu/libfftw3* /usr/lib/x86_64-linux-gnu/
COPY --from=c3d /lib/x86_64-linux-gnu/lib*.so* /lib/x86_64-linux-gnu/
COPY --from=c3d /usr/lib/x86_64-linux-gnu/lib*.so* /usr/lib/x86_64-linux-gnu/
COPY --from=c3d /lib64/ld-linux-x86-64.so.2 /lib64/ld-linux-x86-64.so.2
COPY --from=c3d /usr/bin/c*d /usr/bin/
COPY --from=c3d /usr/bin/c3d_* /usr/bin/
COPY --from=c3d /usr/share/doc/convert3d /usr/share/doc/convert3d
COPY --from=c3d /usr/lib/c3d_gui-1.1.0/Convert3DGUI /usr/lib/c3d_gui-1.1.0/Convert3DGUI

# Installing and setting up FSL 
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

# install Multimodal Surface Matching
COPY --from=MSM /opt/msm/Ubuntu/msm /opt/msm/Ubuntu/msm
ENV MSMBINDIR=/opt/msm/Ubuntu \
    PATH=$PATH:/opt/msm/Ubuntu

# Installing and setting up wb_command
COPY --from=connectome-workbench /lib64/* /lib64/
COPY --from=connectome-workbench /lib/x86_64-linux-gnu/* /lib/x86_64-linux-gnu/
COPY --from=connectome-workbench /usr/* /usr/

# link libraries & clean up
RUN sed -i 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    ldconfig && \
    chmod 777 / && \
    chmod 777 $(ls / | grep -v sys | grep -v proc)

# set user
USER c-pac_user
