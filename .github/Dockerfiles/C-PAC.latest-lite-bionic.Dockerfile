# Choose versions
FROM ghcr.io/fcp-indi/c-pac:latest-bionic

USER root

# remove FreeSurfer
RUN rm -rf /usr/lib/freesurfer/ /code/run-with-freesurfer.sh
ENTRYPOINT ["/code/run.py"]

# link libraries
RUN ldconfig

# clean up
RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set user
USER c-pac_user
