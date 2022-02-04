# Choose versions
FROM ghcr.io/fcp-indi/c-pac/fmriprep:20.2.7-xenial as fmriprep 
# ^ includes AFNI, AFNI, c3d, Freesurfer, FSL, wb_command ^

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

# install ICA-AROMA
COPY --from=fmriprep /opt/ICA-AROMA/ /opt/ICA-AROMA/
ENV PATH=/opt/ICA-AROMA:$PATH

# install FreeSurfer
COPY --from=fmriprep /opt/freesurfer/ /usr/lib/freesurfer/
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH" \
    NO_FSFAST=1

# Installing and setting up AFNI, c3d, FSL & wb_command
COPY --from=fmriprep /usr/local/etc/neurodebian.gpg /usr/local/etc/
# set up AFNI
ENV PATH=/usr/lib/afni/bin:$PATH

# Allow users to install Python packages
RUN chmod -R ugo+w /usr/local/miniconda

# install C-PAC & set up runscript
COPY . /code
COPY dev/docker_data /code/docker_data
RUN pip install git+https://git@github.com/FCP-INDI/INDI-Tools.git@main && \
    pip install -e /code && \
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