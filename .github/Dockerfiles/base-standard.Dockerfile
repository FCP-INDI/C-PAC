FROM ghcr.io/fcp-indi/c-pac/freesurfer:6.0.0-min.neurodocker-bionic as FreeSurfer
FROM ghcr.io/fcp-indi/c-pac/ica-aroma:0.4.3-beta-bionic as ICA-AROMA
FROM ghcr.io/fcp-indi/c-pac/msm:2.0-bionic as MSM

FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Standard software dependencies for C-PAC standard and lite images"
USER root

# Installing connectome-workbench & FSL
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      connectome-workbench=1.5.0-1~nd18.04+1 \
      fsl-core \
      fsl-atlases \
      fsl-mni152-templates && \
    ldconfig && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Installing and setting up c3d
RUN mkdir -p /opt/c3d && \
    curl -sSL "http://downloads.sourceforge.net/project/c3d/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz" \
    | tar -xzC /opt/c3d --strip-components 1
ENV C3DPATH /opt/c3d
ENV PATH $C3DPATH/bin:$PATH

# Installing AFNI
# Installing C-PAC resources into FSL
# Installing ANTs
COPY dev/docker_data/required_afni_pkgs.txt /opt/required_afni_pkgs.txt
ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ANTSPATH=/usr/lib/ants \
    FSLDIR=/usr/share/fsl/5.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/5.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/lib/ants:/opt/afni:/usr/lib/fsl/5.0:$PATH
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
    fi && \
    LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH && \
    export LD_LIBRARY_PATH && \
    curl -O https://afni.nimh.nih.gov/pub/dist/bin/linux_openmp_64/@update.afni.binaries && \
    tcsh @update.afni.binaries -package linux_openmp_64 -bindir /opt/afni -prog_list $(cat /opt/required_afni_pkgs.txt) && \
    ldconfig && \
    curl -sL http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz -o /tmp/cpac_resources.tar.gz && \
    tar xfz /tmp/cpac_resources.tar.gz -C /tmp && \
    cp -n /tmp/cpac_image_resources/MNI_3mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/MNI_4mm/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/symmetric/* $FSLDIR/data/standard && \
    cp -n /tmp/cpac_image_resources/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/2mm $FSLDIR/data/standard/tissuepriors && \
    cp -nr /tmp/cpac_image_resources/tissuepriors/3mm $FSLDIR/data/standard/tissuepriors && \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null \
    | apt-key add - \
    && apt-add-repository -y 'deb https://apt.kitware.com/ubuntu/ bionic main'   \
    && apt-get update \
    && apt-get -y install --no-install-recommends \
      cmake=3.18.3-0kitware1 \
      cmake-data=3.18.3-0kitware1 \
    && mkdir -p /tmp/ants/build \
    && cd /tmp/ants/build \
    && mkdir -p /opt/ants \
    && git config --global url."https://".insteadOf git:// \
    && git clone -b v2.3.5 --depth 1 https://github.com/ANTsX/ANTs.git /tmp/ants/ANTs \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt/ants \
        -DBUILD_SHARED_LIBS=OFF \
        -DUSE_VTK=OFF \
        -DSuperBuild_ANTS_USE_GIT_PROTOCOL=OFF \
        -DBUILD_TESTING=OFF \
        -DRUN_LONG_TESTS=OFF \
        -DRUN_SHORT_TESTS=OFF \
      ../ANTs 2>&1 | tee cmake.log \
    && make -j 4 2>&1 | tee build.log \
    && cd ANTS-build \
    && make install 2>&1 | tee install.log \
    && mkdir /ants_template && \
    curl -sL https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/3133832/Oasis.zip -o /tmp/Oasis.zip && \
    unzip /tmp/Oasis.zip -d /tmp &&\
    mv /tmp/MICCAI2012-Multi-Atlas-Challenge-Data /ants_template/oasis && \
    rm -rf /tmp/Oasis.zip /tmp/MICCAI2012-Multi-Atlas-Challenge-Data

# Installing ICA-AROMA
COPY --from=ICA-AROMA /opt/ICA-AROMA/ /opt/ICA-AROMA/
ENV PATH=/opt/ICA-AROMA:$PATH

# Installing FreeSurfer
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
COPY --from=FreeSurfer /usr/lib/freesurfer/ /usr/lib/freesurfer/
COPY dev/docker_data/license.txt $FREESURFER_HOME/license.txt

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
