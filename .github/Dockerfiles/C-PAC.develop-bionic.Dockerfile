FROM ghcr.io/fcp-indi/c-pac/stage-base:standard-v1.8.6.dev1
LABEL org.opencontainers.image.description "Full C-PAC image"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# install C-PAC
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
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && ldconfig \
    && chmod 777 / \
    && chmod 777 $(ls / | grep -v sys | grep -v proc)
ENV PYTHONUSERBASE=/home/c-pac_user/.local
ENV PATH=$PATH:/home/c-pac_user/.local/bin \
    PYTHONPATH=$PYTHONPATH:$PYTHONUSERBASE/lib/python3.10/site-packages

# set user
# USER c-pac_user
