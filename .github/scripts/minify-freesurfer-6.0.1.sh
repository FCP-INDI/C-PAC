#!/bin/bash
# This script takes a long time (over 24 hours) to run.
# Have $CPAC_DIR set to the C-PAC directory.
# Run from the parent of the data directory.
# with this tree:
# .
# └── data
#     └── fs601
#         ├── fs-subjects
#         └── input
#             ├── ds000239-sub-01_run-01
#             │   └── sub-01_run-01_T1w.nii.gz  # s3://openneuro.org/ds000239/sub-01/anat/sub-01_run-01_T1w.nii.gz
#             ├── ds000239-sub-01_run-02
#             │   └── sub-01_run-02_T1w.nii.gz  # s3://openneuro.org/ds000239/sub-01/anat/sub-01_run-02_T1w.nii.gz
#             └── studyforrest-sub-01
#                 ├── sub-01_T1w.json
#                 ├── sub-01_T1w.nii.gz
#                 ├── sub-01_T2w.json
#                 └── sub-01_T2w.nii.gz

# Set up local env vars
DATA_DIR=${PWD}/data/fs601
INPUT_DIR=${DATA_DIR}/input
SUBJECTS_DIR=${DATA_DIR}/fs-subjects
mkdir -p ${DATA_DIR}
mkdir -p ${SUBJECTS_DIR}
mkdir -p ${INPUT_DIR}

# Delete any FreeSurfer generated files from previous runs
rm -rf ${SUBJECTS_DIR}/*

# Minify
## Invoke the `freesurfer:6.0.1` container:
docker run --rm -dit --security-opt seccomp:unconfined \
  -v ${INPUT_DIR}:/data/input \
  -v ${SUBJECTS_DIR}:/data/fs-subjects \
  -e SUBJECTS_DIR='/data/fs-subjects' \
  -v ${CPAC_DIR}/dev/docker_data/license.txt:/data/license.txt \
  -e FS_LICENSE='/data/license.txt' \
  -e FS_TIME_ALLOW=0 \
  --entrypoint /bin/bash \
  --name fs601 \
  ghcr.io/fcp-indi/c-pac/freesurfer:6.0.1-xenial

## Trace the following commands:
cmd1="recon-all -s studyforrest-sub-01 -all -i /data/input/studyforrest-sub-01/sub-01_T1w.nii.gz  -T2 /data/input/studyforrest-sub-01/sub-01_T2w.nii.gz -T2pial"
cmd2="recon-all -s ds000239-sub-01_run-01 -all -i /data/input/ds000239-sub-01_run-01/sub-01_run-01_T1w.nii.gz"
cmd3="recon-all -s ds000239-sub-01_run-02 -all -i /data/input/ds000239-sub-01_run-02/sub-01_run-02_T1w.nii.gz"
cmd4="python2 $FREESURFER_HOME/bin/asegstats2table --subjects studyforrest-sub-01 --segno 11 17 18 --meas mean --tablefile aseg.mean-intensity.table"
cmd5="python2 $FREESURFER_HOME/bin/aparcstats2table --hemi lh --subjects studyforrest-sub-01 --tablefile lh.aparc.area.table"
cmd6="recon-all -base long-base-ds000239-sub-01 -tp ds000239-sub-01_run-01 -tp ds000239-sub-01_run-02 -all"
cmd7="recon-all -long ds000239-sub-01_run-01 long-base-ds000239-sub-01 -all"
cmd8="mri_vol2vol --all-info"

## With `neurodocker-minify`:
yes | neurodocker-minify --container fs601 \
  --dirs-to-prune /opt \
  --commands "$cmd1" "$cmd2" "$cmd3" "$cmd4" "$cmd5" "$cmd6" "$cmd7" "cmd8" && \
docker export fs601 | docker import - ghcr.io/fcp-indi/c-pac/freesurfer:6.0.1-min-xenial

### Careful with the `yes` command.  This bypasses user verification of the files `neurodocker-minify` has identified for deletion.
