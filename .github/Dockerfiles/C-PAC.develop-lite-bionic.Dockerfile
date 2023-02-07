FROM ghcr.io/fcp-indi/c-pac/stage-base:standard-lite-v1.8.5.dev
LABEL org.opencontainers.image.description "Full C-PAC image without FreeSurfer"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root
ENTRYPOINT ["/code/run.py"]
# set user
USER c-pac_user
