# Choose versions
FROM nipreps/fmriprep:20.2.1 as fmriprep

FROM ghcr.io/fcp-indi/c-pac/afni:16.2.07-xenial as AFNI
FROM ghcr.io/fcp-indi/c-pac/ants:2.3.4.neurodocker-xenial as ANTs
FROM ghcr.io/fcp-indi/c-pac/c3d:1.0.0-xenial as c3d
FROM ghcr.io/fcp-indi/c-pac/fsl:5.0.9-5.neurodebian-xenial as FSL
FROM ghcr.io/fcp-indi/c-pac/freesurfer:6.0.1-xenial as FreeSurfer
FROM ghcr.io/fcp-indi/c-pac/ica-aroma:0.4.3-beta-bionic as ICA-AROMA

# using Ubuntu 16.04 LTS as parent image
FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114

USER root

# allow users to update / create themselves
RUN chmod ugo+w /etc/passwd

# install and set up c3d
COPY --from=c3d /usr/bin/c3d/ /usr/bin/c3d/
COPY --from=c3d /usr/lib/c3d* /usr/lib/
ENV C3DPATH /usr/bin

# install AFNI
COPY --from=AFNI /usr/lib/afni/ /usr/lib/afni/
# set up AFNI
ENV PATH=/usr/lib/afni/bin:$PATH

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
    ANTSPATH="/usr/lib/ants" \
    PATH="/usr/lib/ants:$PATH"
COPY --from=ANTs /ants_template/ /ants_template/
COPY --from=ANTs /etc/locale.gen /etc/
COPY --from=ANTs /usr/lib/ants/ /usr/lib/ants/
RUN dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG="en_US.UTF-8" \
    && chmod 777 /opt && chmod a+s /opt

# install ICA-AROMA
COPY --from=ICA-AROMA /opt/ICA-AROMA/ /opt/ICA-AROMA/
RUN curl -sL https://github.com/rhr-pruim/ICA-AROMA/archive/v0.4.3-beta.tar.gz | tar -xzC /opt/ICA-AROMA --strip-components 1 && \
    chmod +x /opt/ICA-AROMA/ICA_AROMA.py
ENV PATH=/opt/ICA-AROMA:$PATH

# install FreeSurfer
COPY --from=FreeSurfer /usr/lib/freesurfer/ /usr/lib/freesurfer/
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH" \
    NO_FSFAST=1

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
RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    ldconfig && \
    chown c-pac_user /usr/local/miniconda/lib/python3.7/site-packages \
    chmod ugo+w /usr/local/miniconda/lib/python3.7/site-packages \
    chmod 777 $(ls / | grep -v sys | grep -v proc)

# set user
USER c-pac_user