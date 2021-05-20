FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free

USER root

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
    fi && \
    LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH && \
    export LD_LIBRARY_PATH && \
    curl -O https://afni.nimh.nih.gov/pub/dist/bin/linux_openmp_64/@update.afni.binaries && \
    tcsh @update.afni.binaries -package linux_openmp_64 -bindir /opt/afni -prog_list $(cat /opt/required_afni_pkgs.txt) && \
    ldconfig

# set up AFNI
ENV PATH=/opt/afni:$PATH

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user