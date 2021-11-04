# we need mri_vol2vol which is not included in neurodocker freesurfer 6.0.0-min
FROM freesurfer/freesurfer:6.0

# Matching software package versions in https://github.com/DCAN-Labs/abcd-hcp-pipeline/blob/e480a8f99534f1b05f37bf44c64827384b69b383/Dockerfile

# using neurodebian runtime as parent image
FROM neurodebian:bionic-non-free

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
      dh-autoreconf \
      libgsl-dev \
      libmotif-dev \
      libtool \
      libx11-dev \
      libxext-dev \
      x11proto-xext-dev \
      x11proto-print-dev \
      xutils-dev

# Install libpng12
RUN curl -sLo /tmp/libpng12.deb http://mirrors.kernel.org/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1_amd64.deb && \
    dpkg -i /tmp/libpng12.deb && \
    rm /tmp/libpng12.deb

# Install libxp from third-party repository
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository --yes ppa:zeehio/libxp && \
    apt-get update && apt-get install libxp6 libxp-dev && \
    add-apt-repository --remove --yes ppa:zeehio/libxp && \
    apt-get update

# install wb_command v1.3.2
RUN APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 04EE7237B7D453EC && \
    APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 648ACFD622F3D138 && \
    APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key adv --keyserver keyserver.ubuntu.com --recv-keys DCC9EFBF77E11517 && \
    printf '\ndeb http://httpredir.debian.org/debian/ buster main non-free' >> /etc/apt/sources.list && \
    apt-get update && \
    apt-get install connectome-workbench=1.3.2-1 -y && \
    strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

# Installing and setting up c3d
RUN mkdir -p /opt/c3d && \
    curl -sSL "http://downloads.sourceforge.net/project/c3d/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz" \
    | tar -xzC /opt/c3d --strip-components 1
ENV C3DPATH /opt/c3d/
ENV PATH $C3DPATH/bin:$PATH

# install AFNI
COPY dev/docker_data/required_afni_pkgs.txt /opt/required_afni_pkgs.txt
RUN if [ -f /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1.2.0]; then \
        ln -svf /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1.2.0 /usr/lib/x86_64-linux-gnu/libGL.so.1; \
    fi && \
    libs_path=/usr/lib/x86_64-linux-gnu && \
    if [ -f $libs_path/libgsl.so.23 ]; then \
        ln -svf $libs_path/libgsl.so.23 $libs_path/libgsl.so.19 && \
        ln -svf $libs_path/libgsl.so.23 $libs_path/libgsl.so.0; \
    elif [ -f $libs_path/libgsl.so.23.0.0 ]; then \
        ln -svf $libs_path/libgsl.so.23.0.0 $libs_path/libgsl.so.19 && \
        ln -svf $libs_path/libgsl.so.23.0.0 $libs_path/libgsl.so.0; \
    elif [ -f $libs_path/libgsl.so ]; then \
        ln -svf $libs_path/libgsl.so $libs_path/libgsl.so.0; \
    fi
ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
RUN curl -O https://afni.nimh.nih.gov/pub/dist/bin/linux_ubuntu_16_64/@update.afni.binaries && \
    tcsh @update.afni.binaries -package linux_openmp_64 -bindir /opt/afni -prog_list $(cat /opt/required_afni_pkgs.txt) && \
    ldconfig

# set up AFNI
ENV PATH=/opt/afni:$PATH

# install FSL from Neurodocker
RUN apt-get install -y \
    bc \
    dc \
    file \
    libfontconfig1 \
    libfreetype6 \
    libgl1-mesa-dri \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libgomp1 \
    libice6 \
    libxcursor1 \
    libxft2 \
    libxinerama1 \
    libxrandr2 \
    libxrender1 \
    libxt6 \
    wget \
    sudo && \
    echo "Downloading FSL ..." \
    && mkdir -p /usr/share/fsl/5.0 \
    && curl -sSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-5.0.10-centos6_64.tar.gz \
    | tar zx -C /usr/share/fsl/5.0 --strip-components=1 \
    --exclude=fsl/bin/mist \
    --exclude=fsl/bin/possum \
    --exclude=fsl/data/possum \
    --exclude=fsl/data/mist \
    --exclude=fsl/data/first

# setup FSL environment
ENV FSLDIR=/usr/share/fsl/5.0 \
    FSL_DIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/share/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/share/fsl/5.0/bin:$PATH

# install CPAC resources into FSL
RUN curl -sL http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz -o /tmp/cpac_resources.tar.gz && \
    tar xfz /tmp/cpac_resources.tar.gz -C /tmp && \
    cp -n /tmp/cpac_image_resources/MNI_3mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/MNI_4mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/symmetric/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/2mm $FSLDIR/data/standard/tissuepriors && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/3mm $FSLDIR/data/standard/tissuepriors

# install Multimodal Surface Matching
COPY --from=ghcr.io/fcp-indi/c-pac/msm:v2.0-bionic /opt/msm/Ubuntu/msm /opt/msm/Ubuntu/msm
ENV MSMBINDIR=/opt/msm/Ubuntu \
    PATH=$PATH:/opt/msm/Ubuntu

# install ANTs from Neurodocker
ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ANTSPATH="/usr/lib/ants" \
    PATH="/usr/lib/ants:$PATH"

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG="en_US.UTF-8" \
    && chmod 777 /opt && chmod a+s /opt

RUN apt-get install -y \
    cmake \
    g++ \
    gcc \
    git \
    make \
    zlib1g-dev && \
    echo "Downloading ANTs ..." \
    && mkdir -p /usr/lib/ants \
    && curl -fsSL --retry 5 https://dl.dropbox.com/s/2f4sui1z6lcgyek/ANTs-Linux-centos5_x86_64-v2.2.0-0740f91.tar.gz \
    | tar -xz -C /usr/lib/ants --strip-components 1

# download OASIS templates for niworkflows-ants skullstripping
RUN mkdir /ants_template && \
    curl -sL https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/3133832/Oasis.zip -o /tmp/Oasis.zip && \
    unzip /tmp/Oasis.zip -d /tmp &&\
    mv /tmp/MICCAI2012-Multi-Atlas-Challenge-Data /ants_template/oasis && \
    rm -rf /tmp/Oasis.zip /tmp/MICCAI2012-Multi-Atlas-Challenge-Data

# install ICA-AROMA
RUN mkdir -p /opt/ICA-AROMA
RUN curl -sL https://github.com/rhr-pruim/ICA-AROMA/archive/v0.4.3-beta.tar.gz | tar -xzC /opt/ICA-AROMA --strip-components 1
RUN chmod +x /opt/ICA-AROMA/ICA_AROMA.py
ENV PATH=/opt/ICA-AROMA:$PATH

# install miniconda
RUN curl -sO https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh && \
    bash Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-py37_4.8.2-Linux-x86_64.sh

# update path to include conda
ENV PATH=/usr/local/miniconda/bin:$PATH

# install conda dependencies
RUN conda update conda -y && \
    conda install -y  \
        blas \
        cython \
        matplotlib==3.1.3 \
        networkx==2.4 \
        nose==1.3.7 \
        numpy==1.16.4 \
        pandas==0.23.4 \
        scipy==1.4.1 \
        traits==4.6.0 \
        wxpython \
        pip

# install torch
RUN pip install torch==1.2.0 torchvision==0.4.0 -f https://download.pytorch.org/whl/torch_stable.html

# install python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip install --upgrade setuptools
RUN pip install --upgrade pip
RUN pip install --upgrade docutils
RUN pip install -r /opt/requirements.txt
RUN pip install xvfbwrapper

# install PyPEER
RUN pip install git+https://github.com/ChildMindInstitute/PyPEER.git

# install cpac templates
COPY --from=ghcr.io/fcp-indi/c-pac_templates:latest /cpac_templates /cpac_templates

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs
RUN git lfs install

# Get atlases
COPY --from=ghcr.io/fcp-indi/c-pac/neuroparc:v1.0-human /ndmg_atlases /ndmg_atlases

COPY dev/docker_data/default_pipeline.yml /cpac_resources/default_pipeline.yml
COPY dev/circleci_data/pipe-test_ci.yml /cpac_resources/pipe-test_ci.yml

# install FreeSurfer
# set shell to BASH
RUN mkdir -p /usr/lib/freesurfer
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH" \
    NO_FSFAST=1
SHELL ["/bin/bash", "-c"]
RUN apt-get install -y \
    bc \
    libgomp1 \
    libxmu6 \
    libxt6 \
    tcsh \
    perl && \
    curl -fsSL --retry 5 https://dl.dropbox.com/s/nnzcfttc41qvt31/recon-all-freesurfer6-3.min.tgz \
    | tar -xz -C /usr/lib/freesurfer --strip-components 1 && \
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
RUN printf 'source $FREESURFER_HOME/SetUpFreeSurfer.sh' > ~/.bashrc
# restore shell to default (sh)
SHELL ["/bin/sh", "-c"]
COPY dev/docker_data/license.txt $FREESURFER_HOME/license.txt

COPY dev/docker_data /code/docker_data
RUN mv /code/docker_data/* /code && rm -Rf /code/docker_data && chmod +x /code/run-with-freesurfer.sh

COPY --from=0 opt/freesurfer/bin/mri_vol2vol /usr/lib/freesurfer/bin/mri_vol2vol
COPY --from=0 opt/freesurfer/bin/mri_vol2vol.bin /usr/lib/freesurfer/bin/mri_vol2vol.bin

# add DCAN dependencies
RUN mkdir -p /opt/dcan-tools
# DCAN HCP code
RUN git clone -b 'v2.0.0' --single-branch --depth 1 https://github.com/DCAN-Labs/DCAN-HCP.git /opt/dcan-tools/pipeline

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
