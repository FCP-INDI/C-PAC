
#!/bin/bash

#!/bin/bash

echo "START"

StudyFolder="$1"
echo "StudyFolder: ${StudyFolder}"

FreeSurferFolder="$2"
echo "FreeSurferFolder: ${FreeSurferFolder}"

Subject="$3"
echo "Subject: ${Subject}"

T1wRestoreImageCPAC="$4"
echo "T1wRestoreImageCPAC: ${T1wRestoreImageCPAC}"

AtlasSpaceT1wImageCPAC="$5"
echo "AtlasSpaceT1wImageCPAC: ${AtlasSpaceT1wImageCPAC}"

AtlasTransformCPAC="$6"
echo "AtlasTransformCPAC ${AtlasTransformCPAC}"

InverseAtlasTransformCPAC="$7"
echo "InverseAtlasTransformCPAC: ${InverseAtlasTransformCPAC}"

SurfaceAtlasDIR="$8"
echo "SurfaceAtlasDIR: ${SurfaceAtlasDIR}"

GrayordinatesSpaceDIR="$9"
echo "GrayordinatesSpaceDIR: ${GrayordinatesSpaceDIR}"

GrayordinatesResolutions="${10}"
echo "GrayordinatesResolutions: ${GrayordinatesResolutions}"

HighResMesh="${11}"
echo "HighResMesh: ${HighResMesh}"

LowResMeshes="${12}"
echo "LowResMeshes: ${LowResMeshes}"

SubcorticalGrayLabels="${13}"
echo "SubcorticalGrayLabels: ${SubcorticalGrayLabels}"

FreeSurferLabels="${14}"
echo "FreeSurferLabels: ${FreeSurferLabels}"

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


LowResMeshes=${LowResMeshes//@/ }
echo "LowResMeshes: ${LowResMeshes}"

GrayordinatesResolutions=${GrayordinatesResolutions//@/ }
echo "GrayordinatesResolutions: ${GrayordinatesResolutions}"


# ------------------------------------------------------------------------------
#  Load function libraries
# ------------------------------------------------------------------------------


HCPPIPEDIR=/code/CPAC/resources
source ${HCPPIPEDIR}/global/scripts/log.shlib # Logging related functions
echo "HCPPIPEDIR: ${HCPPIPEDIR}"
MSMCONFIGDIR=${HCPPIPEDIR}/MSMConfig

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


# ------------------------------------------------------------------------------
#  Gather and show positional parameters
# ------------------------------------------------------------------------------

#Make some folders for this and later scripts
if [ ! -e "$T1wFolder"/"$NativeFolder" ] ; then
	mkdir -p "$T1wFolder"/"$NativeFolder"
fi
if [ ! -e "$AtlasSpaceFolder"/ROIs ] ; then
	mkdir -p "$AtlasSpaceFolder"/ROIs
fi
if [ ! -e "$AtlasSpaceFolder"/Results ] ; then
	mkdir "$AtlasSpaceFolder"/Results
fi
if [ ! -e "$AtlasSpaceFolder"/"$NativeFolder" ] ; then
	mkdir "$AtlasSpaceFolder"/"$NativeFolder"
fi
if [ ! -e "$AtlasSpaceFolder"/fsaverage ] ; then
	mkdir "$AtlasSpaceFolder"/fsaverage
fi
for LowResMesh in ${LowResMeshes} ; do
	if [ ! -e "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k ] ; then
		mkdir "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k
	fi
	if [ ! -e "$T1wFolder"/fsaverage_LR"$LowResMesh"k ] ; then
		mkdir "$T1wFolder"/fsaverage_LR"$LowResMesh"k
	fi
done
