# Just what we need from fMRIPrep image
FROM nipreps/fmriprep:20.2.7 as fmriprep

USER root

# install AFNI
COPY dev/docker_data/required_afni_pkgs.txt /opt/required_afni_pkgs.txt
RUN rm -rf /usr/lib/afni \
    && echo "Downloading AFNI ..." \
    && mkdir -p /opt/afni-latest \
    && curl -fsSL --retry 5 https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz \
    | tar -xz -C /opt/afni-latest --strip-components 1 \
    --exclude "linux_openmp_64/*.gz" \
    --exclude "linux_openmp_64/funstuff" \
    --exclude "linux_openmp_64/shiny" \
    --exclude "linux_openmp_64/afnipy" \
    --exclude "linux_openmp_64/lib/RetroTS" \
    --exclude "linux_openmp_64/meica.libs" && \
    KEEPERS=$(while read LINE; do echo " -name ${LINE:16} -or "; done < /opt/required_afni_pkgs.txt) \
    && find /opt/afni-latest -type f -not \( ${KEEPERS:-4} \) -delete

FROM ghcr.io/fcp-indi/c-pac/ants:2.3.4.neurodocker-bionic as ANTs
FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114

USER root

COPY --from=fmriprep /opt/afni-latest/ /usr/lib/afni/bin/
COPY --from=fmriprep /usr/lib/ants/ /usr/lib/ants/
COPY --from=fmriprep /opt/ICA-AROMA/ /opt/ICA-AROMA/
COPY --from=fmriprep /opt/freesurfer/ /opt/freesurfer/
COPY --from=fmriprep /usr/local/etc/neurodebian.gpg /usr/local/etc/
COPY --from=ANTs /ants_template/ /ants_template/

USER c-pac_user