# Choose versions
FROM ghcr.io/fcp-indi/c-pac/afni:update.afni.binaries-bionic as AFNI
FROM ghcr.io/fcp-indi/c-pac/ants:2.3.5-bionic as ANTs
FROM ghcr.io/fcp-indi/c-pac/ants:2.3.4.neurodocker-xenial as oldANTS
FROM ghcr.io/fcp-indi/c-pac/c3d:1.0.0-bionic as c3d
FROM ghcr.io/fcp-indi/c-pac/connectome-workbench:1.5.0.neurodebian-bionic as connectome-workbench
FROM ghcr.io/fcp-indi/c-pac/freesurfer:6.0.0-min.neurodocker-bionic as FreeSurfer
FROM ghcr.io/fcp-indi/c-pac/fsl:neurodebian-bionic as FSL
FROM ghcr.io/fcp-indi/c-pac/ica-aroma:0.4.3-beta-bionic as ICA-AROMA
FROM ghcr.io/fcp-indi/c-pac/msm:2.0-bionic as MSM
FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Standard software dependencies for C-PAC standard and lite images"
USER root

# allow users to update / create themselves
RUN chmod ugo+w /etc/passwd

# install and set up c3d
COPY --from=c3d /opt/c3d/ /opt/c3d/
ENV C3DPATH /opt/c3d/
ENV PATH $C3DPATH/bin:$PATH

# install AFNI
COPY --from=AFNI /opt/afni/ /opt/afni/
# set up AFNI
ENV PATH=/opt/afni:$PATH

# install Connectome Workbench
COPY --from=connectome-workbench /usr/ /usr/

# install FSL
COPY --from=FSL /usr/bin/tclsh /usr/bin/wish /usr/bin/
COPY --from=FSL /usr/share/data/ /usr/share/data/
COPY --from=FSL /usr/share/fsl/ /usr/share/fsl/
COPY --from=FSL /usr/lib/ /usr/lib/
# set up FSL environment
ENV FSLDIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/lib/fsl/5.0:$PATH

# install ANTs
ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ANTSPATH="/opt/ants/bin" \
    PATH="/opt/ants/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/ants/lib:$LD_LIBRARY_PATH"
COPY --from=oldANTS /etc/locale.gen /etc/locale.gen
COPY --from=ANTS /ants_template/ /ants_template/
COPY --from=ANTs /opt/ants/ /opt/ants/
RUN dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG="en_US.UTF-8" \
    && chmod 777 /opt && chmod a+s /opt

# install ICA-AROMA
COPY --from=ICA-AROMA /opt/ICA-AROMA/ /opt/ICA-AROMA/
RUN curl -sL https://github.com/rhr-pruim/ICA-AROMA/archive/v0.4.3-beta.tar.gz | tar -xzC /opt/ICA-AROMA --strip-components 1
RUN chmod +x /opt/ICA-AROMA/ICA_AROMA.py
ENV PATH=/opt/ICA-AROMA:$PATH

# install FreeSurfer
COPY --from=FreeSurfer /usr/lib/freesurfer/ /usr/lib/freesurfer/
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH" \
    NO_FSFAST=1
ENV PERL5LIB="$FREESURFER_HOME/mni/share/perl5" \
    FSFAST_HOME="$FREESURFER_HOME/fsfast" \
    SUBJECTS_DIR="$FREESURFER_HOME/subjects" \
    MNI_DIR="$FREESURFER_HOME/mni"
ENV MINC_BIN_DIR="$MNI_DIR/bin" \
    MINC_LIB_DIR="$MNI_DIR/lib" \
    PATH="$PATH:$MINC_BIN_DIR"

# install Multimodal Surface Matching
COPY --from=MSM /opt/msm/Ubuntu/msm /opt/msm/Ubuntu/msm
ENV MSMBINDIR=/opt/msm/Ubuntu \
    PATH=$PATH:/opt/msm/Ubuntu

# link libraries & clean up
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    ldconfig && \
    chmod 777 / && \
    chmod 777 $(ls / | grep -v sys | grep -v proc)

# set user
USER c-pac_user
