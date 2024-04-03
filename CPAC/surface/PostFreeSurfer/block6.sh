
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

HCPPIPEDIR=/code/CPAC/resources
source ${HCPPIPEDIR}/global/scripts/log.shlib # Logging related functions
echo "HCPPIPEDIR: ${HCPPIPEDIR}"
MSMCONFIGDIR=${HCPPIPEDIR}/MSMConfig

for Hemisphere in L R ; do
	#Set a bunch of different ways of saying left and right
	if [ $Hemisphere = "L" ] ; then
		hemisphere="l"
		Structure="CORTEX_LEFT"
	elif [ $Hemisphere = "R" ] ; then
		hemisphere="r"
		Structure="CORTEX_RIGHT"
	fi
	
cd ${StudyFolder}

STRINGII=""
for LowResMesh in ${LowResMeshes} ; do
	STRINGII=$(echo "${STRINGII}${AtlasSpaceFolder}/fsaverage_LR${LowResMesh}k@${AtlasSpaceFolder}/fsaverage_LR${LowResMesh}k@${LowResMesh}k_fs_LR ${T1wFolder}/fsaverage_LR${LowResMesh}k@${AtlasSpaceFolder}/fsaverage_LR${LowResMesh}k@${LowResMesh}k_fs_LR ")
done

#Add CIFTI Maps to Spec Files
for STRING in "$T1wFolder"/"$NativeFolder"@"$AtlasSpaceFolder"/"$NativeFolder"@native "$AtlasSpaceFolder"/"$NativeFolder"@"$AtlasSpaceFolder"/"$NativeFolder"@native "$AtlasSpaceFolder"@"$AtlasSpaceFolder"@"$HighResMesh"k_fs_LR ${STRINGII} ; do
	FolderI=$(echo $STRING | cut -d "@" -f 1)
	FolderII=$(echo $STRING | cut -d "@" -f 2)
	Mesh=$(echo $STRING | cut -d "@" -f 3)
	for STRINGII in sulc@dscalar thickness@dscalar curvature@dscalar aparc@dlabel aparc.a2009s@dlabel ; do #Remove BA@dlabel because it doesn't convert properly
		Map=$(echo $STRINGII | cut -d "@" -f 1)
		Ext=$(echo $STRINGII | cut -d "@" -f 2)
		if [ -e "$FolderII"/"$Subject"."$Map"."$Mesh"."$Ext".nii ] ; then
			wb_command -add-to-spec-file "$FolderI"/"$Subject"."$Mesh".wb.spec INVALID "$FolderII"/"$Subject"."$Map"."$Mesh"."$Ext".nii
		fi
	done
done
done