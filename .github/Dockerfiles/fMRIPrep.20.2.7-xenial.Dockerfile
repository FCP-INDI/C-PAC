# Just what we need from fMRIPrep image
FROM nipreps/fmriprep:20.2.7 as fmriprep
FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114

USER root

COPY --from=fmriprep /usr/lib/ants/ /usr/lib/ants/
COPY --from=fmriprep /opt/ICA-AROMA/ /opt/ICA-AROMA/
COPY --from=fmriprep /usr/lib/freesurfer/ /usr/lib/freesurfer/
COPY --from=fmriprep /usr/local/etc/neurodebian.gpg /usr/local/etc/

USER c-pac_user