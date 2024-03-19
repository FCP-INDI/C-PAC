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

# Source: https://github.com/DCAN-Labs/DCAN-HCP/blob/32254dc/PostFreeSurfer/PostFreeSurferPipeline.sh#L53

# Modifications copyright (C) 2021 - 2024  C-PAC Developers
# This file is part of C-PAC.

# Input Variables
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --post_freesurfer_folder) StudyFolder="$2"; shift ;;
        --subject) Subject="$2"; shift ;;
		--surfatlasdir) SurfaceAtlasDIR="$2"; shift ;;
		--grayordinatesdir) GrayordinatesSpaceDIR="$2"; shift ;;
		--grayordinatesres) GrayordinatesResolutions="$2"; shift ;;
		--hiresmesh) HighResMesh="$2"; shift ;;
		--lowresmesh) LowResMeshes="$2"; shift ;;
		--subcortgraylabels) SubcorticalGrayLabels="$2"; shift ;;
		--freesurferlabels) FreeSurferLabels="$2"; shift ;;
		--freesurfer_folder) FreeSurferFolder="$2"; shift ;;
		--atlas_transform) AtlasTransformCPAC="$2"; shift ;;
		--inverse_atlas_transform) InverseAtlasTransformCPAC="$2"; shift ;;
		--atlas_t1w) AtlasSpaceT1wImageCPAC="$2"; shift ;;
		--t1w_restore) T1wRestoreImageCPAC="$2"; shift ;;
    esac
    shift
done

RegName=MSMSulc
# RegName=FS
useT2=false

# default parameters
CorrectionSigma=$(echo "sqrt ( 200 )" | bc -l)
InflateExtraScale=1

#Naming Conventions
T1wImage="T1w_acpc_dc"
T1wFolder="T1w" #Location of T1w images
T2wFolder="T2w" #Location of T1w images
T2wImage="T2w_acpc_dc"
AtlasSpaceFolder="MNINonLinear"
NativeFolder="Native"
FreeSurferInput="T1w_acpc_dc_restore_1mm"
AtlasTransform="acpc_dc2standard"
InverseAtlasTransform="standard2acpc_dc"
AtlasSpaceT1wImage="T1w_restore"
AtlasSpaceT2wImage="T2w_restore"
T1wRestoreImage="T1w_acpc_dc_restore"
T2wRestoreImage="T2w_acpc_dc_restore"
OrginalT1wImage="T1w"
OrginalT2wImage="T2w"
T1wImageBrainMask="brainmask_fs"
InitialT1wTransform="acpc.mat"
dcT1wTransform="T1w_dc.nii.gz"
InitialT2wTransform="acpc.mat"
dcT2wTransform="T2w_reg_dc.nii.gz"
FinalT2wTransform="${Subject}/mri/transforms/T2wtoT1w.mat"
BiasField="BiasField_acpc_dc"
OutputT1wImage="T1w_acpc_dc"
OutputT1wImageRestore="T1w_acpc_dc_restore"
OutputT1wImageRestoreBrain="T1w_acpc_dc_restore_brain"
OutputMNIT1wImage="T1w"
OutputMNIT1wImageRestore="T1w_restore"
OutputMNIT1wImageRestoreBrain="T1w_restore_brain"
OutputT2wImage="T2w_acpc_dc"
OutputT2wImageRestore="T2w_acpc_dc_restore"
OutputT2wImageRestoreBrain="T2w_acpc_dc_restore_brain"
OutputMNIT2wImage="T2w"
OutputMNIT2wImageRestore="T2w_restore"
OutputMNIT2wImageRestoreBrain="T2w_restore_brain"
OutputOrigT1wToT1w="OrigT1w2T1w.nii.gz"
OutputOrigT1wToStandard="OrigT1w2standard.nii.gz" #File was OrigT2w2standard.nii.gz, regnerate and apply matrix
OutputOrigT2wToT1w="OrigT2w2T1w.nii.gz" #mv OrigT1w2T2w.nii.gz OrigT2w2T1w.nii.gz
OutputOrigT2wToStandard="OrigT2w2standard.nii.gz"
BiasFieldOutput="BiasField"
Jacobian="NonlinearRegJacobians.nii.gz"

T1wFolder="$StudyFolder"/"$T1wFolder"
T2wFolder="$StudyFolder"/"$T2wFolder"
AtlasSpaceFolder="$StudyFolder"/"$AtlasSpaceFolder"
AtlasTransform="$AtlasSpaceFolder"/xfms/"$AtlasTransform"
InverseAtlasTransform="$AtlasSpaceFolder"/xfms/"$InverseAtlasTransform"

# Added by XL
cd $StudyFolder

if [ ! -e "$T1wFolder" ] ; then
	mkdir -p "$T1wFolder"
fi

if [ ! -e "$AtlasSpaceFolder" ] ; then
	mkdir -p "$AtlasSpaceFolder"
fi

if [ ! -e "$AtlasSpaceFolder"/xfms ] ; then
	mkdir -p "$AtlasSpaceFolder"/xfms
fi

if [ ! -e "$T1wFolder"/"$T1wRestoreImage".nii.gz ] ; then
	cp ${T1wRestoreImageCPAC} "$T1wFolder"/"$T1wRestoreImage".nii.gz
fi

if [ ! -e "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz ] ; then
	cp ${AtlasSpaceT1wImageCPAC} "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz
fi

if [ ! -e ${AtlasTransform}.nii.gz ] ; then
	cp ${AtlasTransformCPAC} ${AtlasTransform}.nii.gz
fi

if [ ! -e ${InverseAtlasTransform}.nii.gz ] ; then
	cp ${InverseAtlasTransformCPAC} ${InverseAtlasTransform}.nii.gz
fi

echo "$StudyFolder" "$Subject" "$T1wFolder" "$AtlasSpaceFolder" "$NativeFolder" "$FreeSurferFolder" "$FreeSurferInput" "$T1wRestoreImage" "$T2wRestoreImage" "$SurfaceAtlasDIR" "$HighResMesh" "$LowResMeshes" "$AtlasTransform" "$InverseAtlasTransform" "$AtlasSpaceT1wImage" "$AtlasSpaceT2wImage" "$T1wImageBrainMask" "$FreeSurferLabels" "$GrayordinatesSpaceDIR" "$GrayordinatesResolutions" "$SubcorticalGrayLabels" "$RegName" "$InflateExtraScale" "$useT2"

bash /code/CPAC/surface/PostFreeSurfer/FreeSurfer2CaretConvertAndRegisterNonlinear.sh "$StudyFolder" "$Subject" "$T1wFolder" "$AtlasSpaceFolder" "$NativeFolder" "$FreeSurferFolder" "$FreeSurferInput" "$T1wRestoreImage" "$T2wRestoreImage" "$SurfaceAtlasDIR" "$HighResMesh" "$LowResMeshes" "$AtlasTransform" "$InverseAtlasTransform" "$AtlasSpaceT1wImage" "$AtlasSpaceT2wImage" "$T1wImageBrainMask" "$FreeSurferLabels" "$GrayordinatesSpaceDIR" "$GrayordinatesResolutions" "$SubcorticalGrayLabels" "$RegName" "$InflateExtraScale" "$useT2"
