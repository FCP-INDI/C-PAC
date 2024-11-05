#!/bin/bash

# Human Connectome Project [Pipelines][Pipelines] = THIS SOFTWARE

# Copyright (c) 2011-2014 [The Human Connectome Project][HCP]

# Redistribution and use in source and binary forms, with or without modification,
# is permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions, and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions, and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * The names of Washington University in St. Louis, the University of Minnesota,
#   Oxford University, the Human Connectome Project, or any contributors
#   to this software may *not* be used to endorse or promote products derived
#   from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# <!-- References -->

# [HCP]: http://www.humanconnectome.org
# [Pipelines]: https://github.com/Washington-University/Pipelines

# Source: https://github.com/DCAN-Labs/DCAN-HCP/blob/32254dc/fMRISurface/scripts/SubcorticalProcessing.sh

# Modifications copyright (C) 2021 - 2024  C-PAC Developers
# This file is part of C-PAC.

set -e
script_name="SubcorticalProcessing.sh"
echo "${script_name}: START"

AtlasSpaceFolder="$1"
echo "${script_name}: AtlasSpaceFolder: ${AtlasSpaceFolder}"

ROIFolder="$2"
echo "${script_name}: ROIFolder: ${ROIFolder}"

FinalfMRIResolution="$3"
echo "${script_name}: FinalfMRIResolution: ${FinalfMRIResolution}"

ResultsFolder="$4"
echo "${script_name}: ResultsFolder: ${ResultsFolder}"

NameOffMRI="$5"
echo "${script_name}: NameOffMRI: ${NameOffMRI}"

SmoothingFWHM="$6"
echo "${script_name}: SmoothingFWHM: ${SmoothingFWHM}"

BrainOrdinatesResolution="$7"
echo "${script_name}: BrainOrdinatesResolution: ${BrainOrdinatesResolution}"

VolumefMRI="${ResultsFolder}/${NameOffMRI}"
echo "${script_name}: VolumefMRI: ${VolumefMRI}"

Sigma=`echo "$SmoothingFWHM / ( 2 * ( sqrt ( 2 * l ( 2 ) ) ) )" | bc -l`
echo "${script_name}: Sigma: ${Sigma}"

#NOTE: wmparc has dashes in structure names, which -cifti-create-* won't accept
#ROIs files have acceptable structure names


#deal with fsl_sub being silly when we want to use numeric equality on decimals
unset POSIXLY_CORRECT

## Create scratch space directory to run I/O intensive wb_commands
#TMPDIR=${TMPDIR:-/tmp/$USER}
#if [ ! -d ${TMPDIR} ]; then
#    mkdir -p ${TMDPIR}
#    chmod 770 ${TMPDIR} || true
#fi
#RandomHash=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13 ; echo '')
#TempSubjectDIR="${TMPDIR}/$RandomHash"
#mkdir -p $TempSubjectDIR

#hj edit: make a temp dir
TempSubjectDIR="$ResultsFolder/TempSubjectDIR"
mkdir -p $TempSubjectDIR
chmod 775 $TempSubjectDIR

#function clean_up {
#    echo Exit code caught. Removing temp scratch space directory
#    rm -fR ${TempSubjectDIR}
#}


#generate subject-roi space fMRI cifti for subcortical
if [[ `echo "$BrainOrdinatesResolution == $FinalfMRIResolution" | bc -l | cut -f1 -d.` == "1" ]]
then
    echo "${script_name}: Creating subject-roi subcortical cifti at same resolution as output"
    wb_command -cifti-create-dense-timeseries ${ResultsFolder}/${NameOffMRI}_temp_subject.dtseries.nii -volume "$VolumefMRI".nii.gz "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz
else
    echo "${script_name}: Creating subject-roi subcortical cifti at differing fMRI resolution"
    wb_command -volume-affine-resample "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz $FSLDIR/etc/flirtsch/ident.mat "$VolumefMRI".nii.gz ENCLOSING_VOXEL "$ResultsFolder"/ROIs."$FinalfMRIResolution".nii.gz
    wb_command -cifti-create-dense-timeseries ${ResultsFolder}/${NameOffMRI}_temp_subject.dtseries.nii -volume "$VolumefMRI".nii.gz "$ResultsFolder"/ROIs."$FinalfMRIResolution".nii.gz
    rm -f "$ResultsFolder"/ROIs."$FinalfMRIResolution".nii.gz
fi

echo "${script_name}: Dilating out zeros"
#dilate out any exact zeros in the input data, for instance if the brain mask is wrong
wb_command -cifti-dilate ${ResultsFolder}/${NameOffMRI}_temp_subject.dtseries.nii COLUMN 0 10 ${ResultsFolder}/${NameOffMRI}_temp_subject_dilate.dtseries.nii
rm -f ${ResultsFolder}/${NameOffMRI}_temp_subject.dtseries.nii

echo "${script_name}: Generate atlas subcortical template cifti"
wb_command -cifti-create-label ${ResultsFolder}/${NameOffMRI}_temp_template.dlabel.nii -volume "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz

echo "${script_name}: Running resampling in scratch space to avoid I/O limit"
cp ${ResultsFolder}/${NameOffMRI}_temp_subject_dilate.dtseries.nii ${TempSubjectDIR}/${NameOffMRI}_temp_subject_dilate.dtseries.nii
cp ${ResultsFolder}/${NameOffMRI}_temp_template.dlabel.nii ${TempSubjectDIR}/${NameOffMRI}_temp_template.dlabel.nii

if [[ `echo "${Sigma} > 0" | bc -l | cut -f1 -d.` == "1" ]]
then
    echo "${script_name}: Smoothing and resampling"
    #this is the whole timeseries, so don't overwrite, in order to allow on-disk writing, then delete temporary
    wb_command -cifti-smoothing ${ResultsFolder}/${NameOffMRI}_temp_subject_dilate.dtseries.nii 0 ${Sigma} COLUMN ${TempSubjectDIR}/${NameOffMRI}_temp_subject_smooth.dtseries.nii -fix-zeros-volume
    #resample, delete temporary
    wb_command -cifti-resample ${TempSubjectDIR}/${NameOffMRI}_temp_subject_smooth.dtseries.nii COLUMN ${TempSubjectDIR}/${NameOffMRI}_temp_template.dlabel.nii COLUMN ADAP_BARY_AREA CUBIC ${TempSubjectDIR}/${NameOffMRI}_temp_atlas.dtseries.nii -volume-predilate 10
    #rm -f ${ResultsFolder}/${NameOffMRI}_temp_subject_smooth.dtseries.nii
    cp ${TempSubjectDIR}/${NameOffMRI}_temp_atlas.dtseries.nii ${ResultsFolder}/${NameOffMRI}_temp_atlas.dtseries.nii
else
    echo "${script_name}: Resampling"
    wb_command -cifti-resample ${TempSubjectDIR}/${NameOffMRI}_temp_subject_dilate.dtseries.nii COLUMN ${TempSubjectDIR}/${NameOffMRI}_temp_template.dlabel.nii COLUMN ADAP_BARY_AREA CUBIC ${TempSubjectDIR}/${NameOffMRI}_temp_atlas.dtseries.nii -volume-predilate 10
    cp ${TempSubjectDIR}/${NameOffMRI}_temp_atlas.dtseries.nii ${ResultsFolder}/${NameOffMRI}_temp_atlas.dtseries.nii
fi


#delete the temp space directory
# trap clean_up EXIT SIGTERM SIGHUP SIGINT SIGQUIT
rm -rf ${TempSubjectDIR}

#delete common temporaries
rm -f ${ResultsFolder}/${NameOffMRI}_temp_subject_dilate.dtseries.nii
rm -f ${ResultsFolder}/${NameOffMRI}_temp_template.dlabel.nii

#write output volume, delete temporary
#NOTE: $VolumefMRI contains a path in it, it is not a file in the current directory
wb_command -cifti-separate ${ResultsFolder}/${NameOffMRI}_temp_atlas.dtseries.nii COLUMN -volume-all "$VolumefMRI"_AtlasSubcortical_s"$SmoothingFWHM".nii.gz
rm -f ${ResultsFolder}/${NameOffMRI}_temp_atlas.dtseries.nii

echo "${script_name}: END"
