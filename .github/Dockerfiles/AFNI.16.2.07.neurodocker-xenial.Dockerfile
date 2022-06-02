FROM ubuntu:xenial-20200114 AS AFNI

USER root

# install AFNI
COPY dev/docker_data/required_afni_pkgs.txt /opt/required_afni_pkgs.txt
SHELL ["/bin/bash", "-c"]
RUN apt-get update && \
    apt-get install -y curl && \
    rm -rf /usr/lib/afni \
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
    && find /opt/afni-latest -type f -not \( ${KEEPERS::-4} \) -delete

FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
AFNI 16.2.07 stage"
USER root

COPY --from=AFNI /opt/afni-latest/ /usr/lib/afni/bin/

USER c-pac_user