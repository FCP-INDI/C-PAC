FROM ghcr.io/fcp-indi/c-pac/freesurfer:6.0.1-xenial as FreeSurfer

USER root

COPY .github/scripts/freesurfer-prune /freesurfer-prune
COPY dev/docker_data/required_freesurfer_pkgs.txt /required_freesurfer_pkgs.txt

RUN /freesurfer-prune

FROM scratch

COPY --from=FreeSurfer /opt/freesurfer /opt/freesurfer