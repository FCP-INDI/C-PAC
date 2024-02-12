FROM ghcr.io/fcp-indi/c-pac:latest
LABEL org.opencontainers.image.description "Full C-PAC image without FreeSurfer"
USER root
ENTRYPOINT ["/code/run.py"]

# remove FreeSurfer, link libraries & clean up
RUN rm -rf /usr/lib/freesurfer/ /code/run-with-freesurfer.sh /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    ln -svf /usr/lib/x86_64-linux-gnu/libgsl.so.23 /usr/lib/x86_64-linux-gnu/libgsl.so.0 && ldconfig && \
    chmod 777 $(ls / | grep -v sys | grep -v proc)

# set user
# USER c-pac_user
