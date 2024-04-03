
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

HCPPIPEDIR=/code/CPAC/resources
source ${HCPPIPEDIR}/global/scripts/log.shlib # Logging related functions
echo "HCPPIPEDIR: ${HCPPIPEDIR}"
MSMCONFIGDIR=${HCPPIPEDIR}/MSMConfig

cd ${FreeSurferFolder}
#rm mri/c_ras.mat
MatrixX=$(mri_info mri/brain.finalsurfs.mgz | grep "c_r" | cut -d "=" -f 5 | sed s/" "/""/g)
MatrixY=$(mri_info mri/brain.finalsurfs.mgz | grep "c_a" | cut -d "=" -f 5 | sed s/" "/""/g)
MatrixZ=$(mri_info mri/brain.finalsurfs.mgz | grep "c_s" | cut -d "=" -f 5 | sed s/" "/""/g)
echo "1 0 0 ""$MatrixX" >> mri/c_ras.mat
echo "0 1 0 ""$MatrixY" >> mri/c_ras.mat
echo "0 0 1 ""$MatrixZ" >> mri/c_ras.mat
echo "0 0 0 1" >> mri/c_ras.mat

cd ${StudyFolder}

#Convert FreeSurfer Volumes
for Image in wmparc aparc.a2009s+aseg aparc+aseg ; do
	if [ -e "$FreeSurferFolder"/mri/"$Image".mgz ] ; then

		mri_convert -rt nearest -rl "$T1wFolder"/"$T1wRestoreImage".nii.gz "$FreeSurferFolder"/mri/"$Image".mgz "$T1wFolder"/"$Image"_1mm.nii.gz
		applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$Image".nii.gz
		applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$Image".nii.gz
		wb_command -volume-label-import "$T1wFolder"/"$Image".nii.gz "$FreeSurferLabels" "$T1wFolder"/"$Image".nii.gz -drop-unused-labels
		wb_command -volume-label-import "$AtlasSpaceFolder"/"$Image".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/"$Image".nii.gz -drop-unused-labels
	fi
done

# #Create FreeSurfer Brain Mask (Now done in PostFreeSurfer.sh so brainmask_fs.nii.gz exists for ANTs Registration)
fslmaths "$T1wFolder"/wmparc_1mm.nii.gz -bin -dilD -dilD -dilD -ero -ero "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
wb_command -volume-fill-holes "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
fslmaths "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -bin "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
applywarp --rel --interp=nn -i "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$T1wImageBrainMask".nii.gz
applywarp --rel --interp=nn -i "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$T1wImageBrainMask".nii.gz

#Add volume files to spec files
if $useT2; then
wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$T1wFolder"/"$T2wImage".nii.gz
fi
wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$T1wFolder"/"$T1wRestoreImage".nii.gz

if $useT2; then
wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT2wImage".nii.gz
fi
wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz

if $useT2; then
wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT2wImage".nii.gz
fi
wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz

for LowResMesh in ${LowResMeshes} ; do
  if $useT2; then
  wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT2wImage".nii.gz
  fi
  wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz

  if $useT2; then
  wb_command -add-to-spec-file "$T1wFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec INVALID "$T1wFolder"/"$T2wImage".nii.gz
  fi
  wb_command -add-to-spec-file "$T1wFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec INVALID "$T1wFolder"/"$T1wRestoreImage".nii.gz
done

#Import Subcortical ROIs
for GrayordinatesResolution in ${GrayordinatesResolutions} ; do
  cp "$GrayordinatesSpaceDIR"/Atlas_ROIs."$GrayordinatesResolution".nii.gz "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz
  applywarp --interp=nn -i "$AtlasSpaceFolder"/wmparc.nii.gz -r "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz
  wb_command -volume-label-import "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz -drop-unused-labels
  applywarp --interp=nn -i "$SurfaceAtlasDIR"/Avgwmparc.nii.gz -r "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$AtlasSpaceFolder"/ROIs/Atlas_wmparc."$GrayordinatesResolution".nii.gz
  wb_command -volume-label-import "$AtlasSpaceFolder"/ROIs/Atlas_wmparc."$GrayordinatesResolution".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/ROIs/Atlas_wmparc."$GrayordinatesResolution".nii.gz -drop-unused-labels
  wb_command -volume-label-import "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz ${SubcorticalGrayLabels} "$AtlasSpaceFolder"/ROIs/ROIs."$GrayordinatesResolution".nii.gz -discard-others
  if $useT2; then
  applywarp --interp=spline -i "$AtlasSpaceFolder"/"$AtlasSpaceT2wImage".nii.gz -r "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$AtlasSpaceFolder"/"$AtlasSpaceT2wImage"."$GrayordinatesResolution".nii.gz
  fi
  applywarp --interp=spline -i "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz -r "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage"."$GrayordinatesResolution".nii.gz
done 
