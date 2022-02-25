#

## `freesurfer:7.1.1` 

```
neurodocker generate docker \
  --base continuumio/miniconda:4.7.12 \
  --pkg-manager apt \
  --freesurfer version=7.1.1 \
  --matlabmcr version=2014b install_path=/opt/MCRv84 \
  --run "ln -s /opt/MCRv84/v84 /opt/freesurfer-7.1.1/MCRv84" \
    | docker build -t pwighton/freesurfer:7.1.1 -
```

## `freesurfer:7.1.1-min` 

### Setup local env vars
```
DATA_DIR=${HOME}/data/fs711-20201127
INPUT_DIR=${HOME}/data/fs711-20201127/input
SUBJECTS_DIR=${DATA_DIR}/fs-subjects-20201128
mkdir -p ${DATA_DIR}
mkdir -p ${SUBJECTS_DIR}
mkdir -p ${INPUT_DIR}
```

### Get data

```
mkdir -p ${INPUT_DIR}/ds000239-sub-01_run-01
mkdir -p ${INPUT_DIR}/ds000239-sub-01_run-02
mkdir -p ${INPUT_DIR}/studyforrest-sub-01
aws s3 cp s3://openneuro.org/ds000239/sub-01/anat/sub-01_run-01_T1w.nii.gz ${INPUT_DIR}/ds000239-sub-01_run-01/
aws s3 cp s3://openneuro.org/ds000239/sub-01/anat/sub-01_run-02_T1w.nii.gz ${INPUT_DIR}/ds000239-sub-01_run-02/
rsync -aL psydata.ovgu.de::studyforrest/structural/sub-01/anat/ ${INPUT_DIR}/studyforrest-sub-01
```

### Delete any FreeSurfer generated files
From previous runs

```
rm -rf ${SUBJECTS_DIR}/*
```
### Minify

Generate the `freesurfer:7.1.1` container, then invoke it:
```
docker run --rm -it --security-opt seccomp:unconfined \
  -v ${INPUT_DIR}:/data/input \
  -v ${SUBJECTS_DIR}:/data/fs-subjects \
  -e SUBJECTS_DIR='/data/fs-subjects' \
  -v ${DATA_DIR}/license.txt:/data/license.txt \
  -e FS_LICENSE='/data/license.txt' \
  -e FS_TIME_ALLOW=0 \
  --name fs711 \
  pwighton/freesurfer:7.1.1 \
    /bin/bash
```

Then, in another terminal, trace the following commands:
```
cmd1="recon-all -s studyforrest-sub-01 -all -i /data/input/studyforrest-sub-01/sub-01_T1w.nii.gz  -T2 /data/input/studyforrest-sub-01/sub-01_T2w.nii.gz -T2pial"
cmd2="recon-all -s ds000239-sub-01_run-01 -all -i /data/input/ds000239-sub-01_run-01/sub-01_run-01_T1w.nii.gz"
cmd3="recon-all -s ds000239-sub-01_run-02 -all -i /data/input/ds000239-sub-01_run-02/sub-01_run-02_T1w.nii.gz"
cmd4="asegstats2table --subjects studyforrest-sub-01 --segno 11 17 18 --meas mean --tablefile aseg.mean-intensity.table"
cmd5="aparcstats2table --hemi lh --subjects studyforrest-sub-01 --tablefile lh.aparc.area.table"
cmd6="recon-all -base long-base-ds000239-sub-01 -tp ds000239-sub-01_run-01 -tp ds000239-sub-01_run-02 -all"
cmd7="recon-all -long ds000239-sub-01_run-01 long-base-ds000239-sub-01 -all"
```

With `neurodocker-minify`:
```
yes | neurodocker-minify --container fs711 \
  --dirs-to-prune /opt \
  --commands "$cmd1" "$cmd2" "$cmd3" "$cmd4" "$cmd5" "$cmd6" "$cmd7" && \
docker export fs711 | docker import - pwighton/freesurfer:7.1.1-min
```

Careful with the `yes` command.  This bypasses user verification of the files `neurodocker-minify` has identified for deletion.
