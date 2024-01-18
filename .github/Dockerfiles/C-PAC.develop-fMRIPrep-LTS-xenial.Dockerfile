FROM ghcr.io/fcp-indi/c-pac/stage-base:fmriprep-lts-v1.8.7dev1
LABEL org.opencontainers.image.description "Full C-PAC image with software dependencies version-matched to [fMRIPrep LTS](https://reproducibility.stanford.edu/fmriprep-lts#long-term-support-lts)"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# install C-PAC & set up runscript
COPY dev/circleci_data/pipe-test_ci.yml /cpac_resources/pipe-test_ci.yml
COPY . /code
RUN pip cache purge && pip install -e /code
# set up runscript
COPY dev/docker_data /code/docker_data
RUN rm -Rf /code/docker_data/checksum && \
    mv /code/docker_data/* /code && \
    rm -Rf /code/docker_data && \
    chmod +x /code/run.py && \
    chmod +x /code/run-with-freesurfer.sh
ENTRYPOINT ["/code/run-with-freesurfer.sh"]

# link libraries & clean up
RUN sed -i 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && locale-gen \
    && apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/pip/* \
    && ldconfig \
    && chmod 777 / \
    && chmod -R 777 /home/c-pac_user \
    && chmod 777 $(ls / | grep -v sys | grep -v proc)
ENV PYTHONUSERBASE=/home/c-pac_user/.local
ENV PATH=$PATH:/home/c-pac_user/.local/bin \
    PYTHONPATH=$PYTHONPATH:$PYTHONUSERBASE/lib/python3.10/site-packages

# set user
WORKDIR /home/c-pac_user
# USER c-pac_user
