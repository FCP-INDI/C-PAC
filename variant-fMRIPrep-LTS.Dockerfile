# we need mri_vol2vol which is not included in neurodocker freesurfer 6.0.0-min
FROM freesurfer/freesurfer:6.0 AS freesurfer

FROM nipreps/fmriprep:20.2.1 as fmriprep

# using Ubuntu 16.04 LTS as parent image
FROM ubuntu:xenial-20200114

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# Install the validator
RUN apt-get install -y apt-utils curl && \
     curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.38.0/install.sh | bash

RUN export NVM_DIR=$HOME/.nvm && \
     [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh" && \
     [ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion" && \
     nvm install 16.8.0 && \
     nvm use 16.8.0 && \
     nvm alias default 16.8.0 && \
     mkdir /root/.npm-packages && \
     npm config set prefix /root/.npm-packages && \
     NPM_PACKAGES=/root/.npm-packages && \
     npm install -g bids-validator

ENV PATH=/root/.npm-packages/bin:$PATH

# Install Ubuntu dependencies and utilities
RUN apt-get install -y \
      build-essential \
      bzip2 \
      ca-certificates \
      cmake \
      git \
      graphviz \
      graphviz-dev \
      gsl-bin \
      libcanberra-gtk-module \
      libexpat1-dev \
      libgiftiio-dev \
      libglib2.0-dev \
      libglu1-mesa \
      libglu1-mesa-dev \
      libjpeg-progs \
      libgl1-mesa-dri \
      libglw1-mesa \
      libxml2 \
      libxml2-dev \
      libxext-dev \
      libxft2 \
      libxft-dev \
      libxi-dev \
      libxmu-headers \
      libxmu-dev \
      libxpm-dev \
      libxslt1-dev \
      locales \
      m4 \
      make \
      mesa-common-dev \
      mesa-utils \
      netpbm \
      pkg-config \
      rsync \
      tcsh \
      unzip \
      vim \
      xvfb \
      xauth \
      zlib1g-dev

# Install 16.04 dependencies
RUN apt-get install -y \
      connectome-workbench \
      dh-autoreconf \
      libgsl-dev \
      libmotif-dev \
      libtool \
      libx11-dev \
      libxext-dev \
      python3 \
      x11proto-xext-dev \
      x11proto-print-dev \
      xutils-dev


# Install libxp from third-party repository
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository --yes ppa:zeehio/libxp && \
    apt-get update && apt-get install libxp6 libxp-dev && \
    add-apt-repository --remove --yes ppa:zeehio/libxp && \
    apt-get update

# Installing and setting up c3d, AFNI, FSL, wb_command
COPY --from=fmriprep /usr/local/etc/neurodebian.gpg /usr/local/etc/

RUN curl -sSL "http://neuro.debian.net/lists/$( lsb_release -c | cut -f2 ).us-ca.full" >> /etc/apt/sources.list.d/neurodebian.sources.list && \
    APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key add /usr/local/etc/neurodebian.gpg && \
    (apt-key adv --refresh-keys --keyserver hkp://ha.pool.sks-keyservers.net 0xA5D32F012649A5A9 || true)

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    fsl-core=5.0.9-5~nd16.04+1 \
                    fsl-mni152-templates=5.0.7-2 \
                    afni=16.2.07~dfsg.1-5~nd16.04+1 \
                    convert3d \
                    connectome-workbench=1.3.2-2~nd16.04+1 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY --from=fcpindi/c-pac:nightly /opt/afni/3dTto1D /usr/lib/afni/bin/
COPY --from=fcpindi/c-pac:nightly /lib/x86_64-linux-gnu/libm.so.6 /lib/x86_64-linux-gnu/

# set up AFNI
ENV PATH=/usr/lib/afni/bin:$PATH

# setup FSL environment
ENV FSLDIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/lib/fsl/5.0:$PATH

# install CPAC resources into FSL
RUN curl -sL http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz -o /tmp/cpac_resources.tar.gz && \
    tar xfz /tmp/cpac_resources.tar.gz -C /tmp && \
    mkdir -p $FSLDIR/data/atlases && \
    cp -n /tmp/cpac_image_resources/MNI_3mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/MNI_4mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/symmetric/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/2mm $FSLDIR/data/standard/tissuepriors && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/3mm $FSLDIR/data/standard/tissuepriors

# install ANTs from Neurodocker
ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ANTSPATH="/usr/lib/ants" \
    PATH="/usr/lib/ants:$PATH"

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG="en_US.UTF-8" \
    && chmod 777 /opt && chmod a+s /opt

RUN echo "Downloading ANTs ..." \
    && mkdir -p /usr/lib/ants \
    && curl -fsSL --retry 5 https://dl.dropbox.com/s/gwf51ykkk5bifyj/ants-Linux-centos6_x86_64-v2.3.4.tar.gz \
    | tar -xz -C /usr/lib/ants --strip-components 1

# download OASIS templates for niworkflows-ants skullstripping
RUN mkdir /ants_template && \
    curl -sL https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/3133832/Oasis.zip -o /tmp/Oasis.zip && \
    unzip /tmp/Oasis.zip -d /tmp &&\
    mv /tmp/MICCAI2012-Multi-Atlas-Challenge-Data /ants_template/oasis && \
    rm -rf /tmp/Oasis.zip /tmp/MICCAI2012-Multi-Atlas-Challenge-Data

# install ICA-AROMA
RUN mkdir -p /opt/ICA-AROMA
RUN curl -sL https://github.com/oesteban/ICA-AROMA/archive/v0.4.5.tar.gz| tar -xzC /opt/ICA-AROMA --strip-components 1
RUN chmod +x /opt/ICA-AROMA/ICA_AROMA.py
ENV PATH=/opt/ICA-AROMA:$PATH

# Installing and setting up miniconda
RUN curl -sSLO https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh && \
    bash Miniconda3-4.5.11-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-4.5.11-Linux-x86_64.sh

# Set CPATH for packages relying on compiled libs (e.g. indexed_gzip)
ENV PATH="/usr/local/miniconda/bin:$PATH" \
    CPATH="/usr/local/miniconda/include/:$CPATH" \
    LANG="C.UTF-8" \
    LC_ALL="C.UTF-8" \
    PYTHONNOUSERSITE=1

# install conda dependencies
RUN conda update conda -y && \
    conda install nomkl && \
    conda install -y  \
        blas \
        cython \
        matplotlib==2.2.2 \
        networkx==2.4 \
        nose==1.3.7 \
        numpy==1.15.4 \
        pandas==0.23.4 \
        scipy==1.1.0 \
        traits==4.6.0 \
        pip

# install torch
RUN pip install torch==1.2.0 torchvision==0.4.0 -f https://download.pytorch.org/whl/torch_stable.html

# install python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip install --upgrade setuptools
RUN pip install --upgrade pip
RUN pip install -r /opt/requirements.txt
RUN pip install xvfbwrapper

# install PyPEER
RUN pip install git+https://github.com/ChildMindInstitute/PyPEER.git

# install cpac templates
ADD dev/docker_data/cpac_templates.tar.gz /

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs
RUN git lfs install

# Get atlases
COPY --from=ghcr.io/fcp-indi/c-pac/neuroparc:v1.0-human /ndmg_atlases /ndmg_atlases

COPY dev/docker_data/default_pipeline.yml /cpac_resources/default_pipeline.yml
COPY dev/circleci_data/pipe-test_ci.yml /cpac_resources/pipe-test_ci.yml

# install FreeSurfer
# set shell to BASH
SHELL ["/bin/bash", "-c"]
ENV FREESURFER_HOME=/usr/lib/freesurfer
RUN mkdir -p /usr/lib/freesurfer && \
    curl -sSL https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.1/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.1.tar.gz | tar zxv --no-same-owner -C /usr/lib \
    --exclude='freesurfer/diffusion' \
    --exclude='freesurfer/docs' \
    --exclude='freesurfer/fsfast' \
    --exclude='freesurfer/lib/cuda' \
    --exclude='freesurfer/lib/qt' \
    --exclude='freesurfer/matlab' \
    --exclude='freesurfer/mni/share/man' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/bert' \
    --exclude='freesurfer/subjects/lh.EC_average' \
    --exclude='freesurfer/subjects/rh.EC_average' \
    --exclude='freesurfer/subjects/sample-*.mgz' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/trctrain' && \
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
RUN printf 'source $FREESURFER_HOME/SetUpFreeSurfer.sh' > ~/.bashrc
# restore shell to default (sh)
SHELL ["/bin/sh", "-c"]
COPY dev/docker_data/license.txt $FREESURFER_HOME/license.txt

COPY dev/docker_data /code/docker_data
RUN mv /code/docker_data/* /code && rm -Rf /code/docker_data && chmod +x /code/run-with-freesurfer.sh

COPY --from=freesurfer opt/freesurfer/bin/mri_vol2vol /usr/lib/freesurfer/bin/mri_vol2vol
COPY --from=freesurfer opt/freesurfer/bin/mri_vol2vol.bin /usr/lib/freesurfer/bin/mri_vol2vol.bin

COPY . /code
RUN pip install -e /code

COPY dev/docker_data /code/docker_data
RUN mv /code/docker_data/* /code && rm -Rf /code/docker_data && chmod +x /code/run.py

ENTRYPOINT ["/code/run-with-freesurfer.sh"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
