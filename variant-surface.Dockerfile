FROM fcpindi/c-pac:latest
MAINTAINER The C-PAC Team <cnl@childmind.org>

# install FreeSurfer
# set shell to BASH
ENV FREESURFER_HOME="/usr/lib/freesurfer" \
    PATH="/usr/lib/freesurfer/bin:$PATH"
SHELL ["/bin/bash", "-c"]
RUN curl -fsSL --retry 5 https://dl.dropbox.com/s/c3earkfhhvdyuo4/freesurfer-7.1.1-min.tgz | \
    tar -xz -C /usr/lib/freesurfer --strip-components 1 && \
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
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