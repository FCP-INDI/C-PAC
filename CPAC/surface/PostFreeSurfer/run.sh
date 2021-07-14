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
        # *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# StudyFolder=`opts_GetOpt1 "--post_freesurfer_folder" $@`
# Subject=`opts_GetOpt1 "--subject" $@`
# SurfaceAtlasDIR=`opts_GetOpt1 "--surfatlasdir" $@`
# GrayordinatesSpaceDIR=`opts_GetOpt1 "--grayordinatesdir" $@`
# GrayordinatesResolutions=`opts_GetOpt1 "--grayordinatesres" $@`
# HighResMesh=`opts_GetOpt1 "--hiresmesh" $@`
# LowResMeshes=`opts_GetOpt1 "--lowresmesh" $@`
# SubcorticalGrayLabels=`opts_GetOpt1 "--subcortgraylabels" $@`
# FreeSurferLabels=`opts_GetOpt1 "--freesurferlabels" $@`
# ReferenceMyelinMaps=`opts_GetOpt1 "--refmyelinmaps" $@`
# CorrectionSigma=`opts_GetOpt1 "--mcsigma" $@`
RegName=MSMSulc
# InflateExtraScale=`opts_GetOpt1 "--inflatescale" $@`
useT2=false

# Extra arguments for ANTs based Atlas Registration
# useStudyTemplate=`opts_GetOpt1 "--useStudyTemplate" $@`
# StudyTemplate=`opts_GetOpt1 "--studytemplate" $@`
# StudyTemplateBrain=`opts_GetOpt1 "--studytemplatebrain" $@`
# T1wTemplate=`opts_GetOpt1 "--t1template" $@`
# T1wTemplateBrain=`opts_GetOpt1 "--t1templatebrain" $@`
# T1wTemplate2mm=`opts_GetOpt1 "--t1template2mm" $@`
# T1wTemplate2mmBrain=`opts_GetOpt1 "--t1template2mmbrain" $@`
# T2wTemplate=`opts_GetOpt1 "--t2template" $@`
# T2wTemplateBrain=`opts_GetOpt1 "--t2templatebrain" $@`
# T2wTemplate2mm=`opts_GetOpt1 "--t2template2mm" $@`
# TemplateMask=`opts_GetOpt1 "--templatemask" $@`
# Template2mmMask=`opts_GetOpt1 "--template2mmmask" $@`

#################### ABIDE FIX ########################
# Reference2mm=`opts_GetOpt1 "--reference2mm" $@`
# Reference2mmMask=`opts_GetOpt1 "--reference2mmmask" $@`
# FNIRTConfig=`opts_GetOpt1 "--config" $@`
#######################################################

# log_Msg "RegName: ${RegName}"

# default parameters
CorrectionSigma=$(echo "sqrt ( 200 )" | bc -l)
# RegName=FS
InflateExtraScale=1

# PipelineScripts=${HCPPIPEDIR_PostFS}

#Naming Conventions
T1wImage="T1w_acpc_dc"
T1wFolder="T1w" #Location of T1w images
T2wFolder="T2w" #Location of T1w images
T2wImage="T2w_acpc_dc" 
AtlasSpaceFolder="MNINonLinear"
NativeFolder="Native"
# FreeSurferFolder=`opts_GetOpt1 "--freesurfer_folder" $@`
FreeSurferInput="T1w_acpc_dc_restore_1mm"
# AtlasTransformCPAC=`opts_GetOpt1 "--atlas_transform" $@`
# /data3/cnl/xli/reproducibility/out/cpac_abcd/output/cpac_cpac_abcd-options/sub-0025427a_ses-1/anat/sub-0025427a_ses-1_from-T1w_to-template_mode-image_xfm.nii.gz # from-T1w_to-template_mode-image_xfm
AtlasTransform="acpc_dc2standard"
# InverseAtlasTransformCPAC=`opts_GetOpt1 "--inverse_atlas_transform" $@`
# /data3/cnl/cpac1.8.1/abcd_scout/output/cpac_cpac_abcd-options/sub-0025427_ses-1/anat/sub-0025427_ses-1_from-template_to-T1w_mode-image_xfm.nii.gz
# InverseAtlasTransformABCD=/data3/cnl/fmriprep/Lei_working/CPAC_XCP/ABCD/ABCD_Testing_With_Intermediate/sub-0025427/ses-1/files/MNINonLinear/xfms/standard2acpc_dc.nii.gz
InverseAtlasTransform="standard2acpc_dc"
# AtlasSpaceT1wImageCPAC=`opts_GetOpt1 "--atlas_t1w" $@`
# /data3/cnl/xli/reproducibility/out/cpac_abcd/working/cpac_sub-0025427a_ses-1/FSL-ABCD_T1_to_template_81/sub-0025427a_run-1_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz
AtlasSpaceT1wImage="T1w_restore"
AtlasSpaceT2wImage="T2w_restore"
# T1wRestoreImageCPAC=`opts_GetOpt1 "--t1w_restore" $@`
# /data3/cnl/xli/reproducibility/out/cpac_abcd/working/cpac_sub-0025427a_ses-1/fast_bias_field_correction_52/get_anat_restore/sub-0025427a_run-1_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths.nii.gz # DCAN origin file: "T1w_acpc_dc_restore"; CPAC RP item: "desc-restore-brain_T1w"
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
# FreeSurferFolder="$T1wFolder"/"$FreeSurferFolder"
# FreeSurferFolder=/data3/cnl/xli/reproducibility/out/cpac_abcd/working/cpac_sub-0025427a_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all
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