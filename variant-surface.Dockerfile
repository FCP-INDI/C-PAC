FROM fcpindi/c-pac:latest
MAINTAINER The C-PAC Team <cnl@childmind.org>

# install FreeSurfer
# set shell to BASH
ENV FREESURFER_HOME=/usr/lib/freesurfer
SHELL ["/bin/bash", "-c"]
RUN curl https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/dev/freesurfer-linux-centos6_x86_64-dev.tar.gz -o /usr/lib/freesurfer.tar.gz && \
    tar -xzvf /usr/lib/freesurfer.tar.gz -C /usr/lib && \
    source $FREESURFER_HOME/SetUpFreeSurfer.sh && \
    rm /usr/lib/freesurfer.tar.gz
RUN printf 'source $FREESURFER_HOME/SetUpFreeSurfer.sh' > ~/.bashrc
# restore shell to default (sh)
SHELL ["/bin/sh", "-c"]
COPY dev/docker_data/license.txt $FREESURFER_HOME/license.txt

COPY dev/docker_data /code/docker_data
RUN mv /code/docker_data/* /code && rm -Rf /code/docker_data && chmod +x /code/run-with-freesurfer.sh

ENTRYPOINT ["/code/run-with-freesurfer.sh"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*