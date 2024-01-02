FROM ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free AS FSL

USER root

# install CPAC resources into FSL
COPY dev/docker_data/checksum/FSL.data.sha384 /tmp/checksum.sha384
RUN mkdir -p /fsl_data/atlases/HarvardOxford fsl_data/standard/tissuepriors \
    && curl -sL http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz -o /tmp/cpac_resources.tar.gz \
    && sha384sum --check /tmp/checksum.sha384 \
    && tar xfz /tmp/cpac_resources.tar.gz -C /tmp \
    && cp -n /tmp/cpac_image_resources/MNI_3mm/* /fsl_data/standard \
    && cp -n /tmp/cpac_image_resources/MNI_4mm/* /fsl_data/standard \
    && cp -n /tmp/cpac_image_resources/symmetric/* /fsl_data/standard \
    && cp -n /tmp/cpac_image_resources/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz /fsl_data/atlases/HarvardOxford \
    && cp -nr /tmp/cpac_image_resources/tissuepriors/2mm /fsl_data/standard/tissuepriors \
    && cp -nr /tmp/cpac_image_resources/tissuepriors/3mm /fsl_data/standard/tissuepriors \
    && chmod -R ugo+r /fsl_data/standard \
    && chmod -R ugo+r /fsl_data/atlases

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
FSL data"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC

COPY --from=FSL /fsl_data/standard fsl_data/standard
COPY --from=FSL /fsl_data/atlases fsl_data/atlases
