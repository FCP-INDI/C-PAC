FROM ghcr.io/fcp-indi/c-pac/freesurfer:6.0.1-xenial as FreeSurfer
USER root
COPY .github/scripts/freesurfer-prune /freesurfer-prune
COPY dev/docker_data/required_freesurfer_pkgs.txt /required_freesurfer_pkgs.txt
RUN /freesurfer-prune

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FreeSurfer 6.0.1-min stage"
COPY --from=FreeSurfer /opt/freesurfer /opt/freesurfer
COPY dev/docker_data/license.txt /opt/freesurfer/license.txt