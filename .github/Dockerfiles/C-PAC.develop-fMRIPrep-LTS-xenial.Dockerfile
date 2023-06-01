FROM ghcr.io/fcp-indi/c-pac/stage-base:fmriprep-lts-v1.8.5.dev
LABEL org.opencontainers.image.description "Full C-PAC image with software dependencies version-matched to [fMRIPrep LTS](https://reproducibility.stanford.edu/fmriprep-lts#long-term-support-lts)"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# install C-PAC & set up runscript
COPY dev/circleci_data/pipe-test_ci.yml /cpac_resources/pipe-test_ci.yml
COPY . /code
RUN pip install -e /code
# set up runscript
COPY dev/docker_data /code/docker_data
RUN rm -Rf /code/docker_data/Dockerfiles && \
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
# USER c-pac_user
