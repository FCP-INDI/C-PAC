# Choose versions
FROM ghcr.io/fcp-indi/c-pac/afni:16.2.07.neurodocker-xenial as AFNI
FROM ghcr.io/fcp-indi/c-pac/fmriprep:20.2.7-xenial as fmriprep 
# ^ includes ANTs, c3d, Freesurfer, FSL, wb_command ^

# using Ubuntu 16.04 LTS as parent image
FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114

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
    PATH=/usr/lib/fsl/5.0:$PATH

# install ANTs
ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ANTSPATH="/usr/lib/ants" \
    PATH="/usr/lib/ants:$PATH"
COPY --from=fmriprep /usr/lib/ants/ /usr/lib/ants/
COPY --from=fmriprep /ants_template /ants_template

# install ICA-AROMA
COPY --from=fmriprep /opt/ICA-AROMA/ /opt/ICA-AROMA/
ENV PATH=/opt/ICA-AROMA:$PATH

# install FreeSurfer
COPY --from=fmriprep /opt/freesurfer/ /usr/lib/freesurfer/
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH" \
    NO_FSFAST=1

# Installing and setting up AFNI
COPY --from=AFNI /usr/lib/afni/bin/ /usr/lib/afni/bin/
# set up AFNI
ENV PATH=/usr/lib/afni/bin:$PATH \
    AFNI_MODELPATH="/usr/lib/afni/models" \
    AFNI_IMSAVE_WARNINGS="NO" \
    AFNI_TTATLAS_DATASET="/usr/share/afni/atlases" \
    AFNI_PLUGINPATH="/usr/lib/afni/plugins"

# Intalling and setting up c3d
COPY --from=fmriprep /usr/bin/c*d /usr/bin/
COPY --from=fmriprep /usr/share/doc/convert3d /usr/share/doc/convert3d
COPY --from=fmriprep /usr/lib/c3d_gui-1.1.0/Convert3DGUI /usr/lib/c3d_gui-1.1.0/Convert3DGUI
# Installing and setting up FSL 
COPY --from=fmriprep /etc/fsl /etc/fsl
COPY --from=fmriprep /usr/share/doc/fsl-core /usr/share/doc/fsl-core
COPY --from=fmriprep /usr/share/man/man1/fsl* /usr/share/man/man1/
COPY --from=fmriprep /usr/share/data/fsl-mni152-templates /usr/share/data/fsl-mni152-templates
COPY --from=fmriprep /usr/share/doc/fsl-mni152-templates /usr/share/doc/fsl-mni152-templates
COPY --from=fmriprep /usr/share/fsl /usr/share/fsl
ENV FSLDIR="/usr/share/fsl/5.0" \
    FSLOUTPUTTYPE="NIFTI_GZ" \
    FSLMULTIFILEQUIT="TRUE" \
    POSSUMDIR="/usr/share/fsl/5.0" \
    LD_LIBRARY_PATH="/usr/lib/fsl/5.0:$LD_LIBRARY_PATH" \
    FSLTCLSH="/usr/bin/tclsh" \
    FSLWISH="/usr/bin/wish"
# Installing and setting up wb_command
COPY --from=fmriprep /usr/bin/wb_* /usr/bin/
COPY --from=fmriprep /usr/share/applications/connectome-workbench.desktop /usr/share/applications/connectome-workbench.desktop
COPY --from=fmriprep /usr/share/bash-completion/completions/wb* /usr/share/bash_completion/completions/
COPY --from=fmriprep /usr/share/doc/connectome-workbench /usr/share/doc/connectome-workbench
COPY --from=fmriprep /usr/share/man/man1/wb_* /usr/share/man/man1/
COPY --from=fmriprep /usr/share/pixmaps/connectome-workbench.png /usr/share/pixmaps/connectome-workbench.png

# Allow users to install Python packages
RUN chmod -R ugo+w /usr/local/miniconda

# install C-PAC & set up runscript
COPY . /code
COPY dev/docker_data /code/docker_data
RUN pip install -e /code && \
    rm -Rf /code/docker_data/Dockerfiles && \
    mv /code/docker_data/* /code && \
    rm -Rf /code/docker_data && \
    chmod +x /code/run.py && \
    chmod +x /code/run-with-freesurfer.sh
ENTRYPOINT ["/code/run-with-freesurfer.sh"]

# link libraries & clean up
RUN sed -i 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    ldconfig && \
    chmod 777 $(ls / | grep -v sys | grep -v proc)

# set user
USER c-pac_user