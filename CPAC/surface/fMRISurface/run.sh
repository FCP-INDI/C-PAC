Path=`opts_GetOpt1 "--post_freesurfer_folder" $@`
Subject=`opts_GetOpt1 "--subject" $@`
NameOffMRI=task-rest01
LowResMesh=`opts_GetOpt1 "--lowresmesh" $@`
FinalfMRIResolution=`opts_GetOpt1 "--fmrires" $@`
SmoothingFWHM=`opts_GetOpt1 "--smoothingFWHM" $@`
GrayordinatesResolution=`opts_GetOpt1 "--grayordinatesres" $@`
RegName=MSMSulc

AtlasSpaceFolder="MNINonLinear"
T1wFolder="T1w"
NativeFolder="Native"
ResultsFolder="Results"
DownSampleFolder="fsaverage_LR${LowResMesh}k"
ROIFolder="ROIs"
OutputAtlasDenseTimeseries="${NameOffMRI}_Atlas"

AtlasSpaceFolder="$Path"/"$AtlasSpaceFolder"
T1wFolder="$Path"/"$T1wFolder"
ResultsFolder="$AtlasSpaceFolder"/"$ResultsFolder"/"$NameOffMRI" # we don't have this folder so far, make one
ROIFolder="$AtlasSpaceFolder"/"$ROIFolder"

# func to standard
AtlasSpacefMRI=`opts_GetOpt1 "--fmri" $@`
Scout=`opts_GetOpt1 "--scout" $@`

if [ ! -e "$ResultsFolder" ] ; then
	mkdir "$ResultsFolder"
fi

if [ ! -e "$ResultsFolder"/"$NameOffMRI".nii.gz ] ; then
	cp ${AtlasSpacefMRI} "$ResultsFolder"/"$NameOffMRI".nii.gz
fi

if [ ! -e "$ResultsFolder"/"$NameOffMRI"_SBRef.nii.gz ] ; then
	cp ${Scout} "$ResultsFolder"/"$NameOffMRI"_SBRef.nii.gz
fi

#Make fMRI Ribbon
#Noisy Voxel Outlier Exclusion
#Ribbon-based Volume to Surface mapping and resampling to standard surface
# log_Msg "Make fMRI Ribbon"
# log_Msg "mkdir -p ${ResultsFolder}/RibbonVolumeToSurfaceMapping"
# mkdir -p "$ResultsFolder"/RibbonVolumeToSurfaceMapping
bash RibbonVolumeToSurfaceMapping.sh "$ResultsFolder"/RibbonVolumeToSurfaceMapping "$ResultsFolder"/"$NameOffMRI" "$Subject" "$AtlasSpaceFolder"/"$DownSampleFolder" "$LowResMesh" "$AtlasSpaceFolder"/"$NativeFolder" "${RegName}"

#Surface Smoothing
# log_Msg "Surface Smoothing"
bash SurfaceSmoothing.sh "$ResultsFolder"/"$NameOffMRI" "$Subject" "$AtlasSpaceFolder"/"$DownSampleFolder" "$LowResMesh" "$SmoothingFWHM"

#Subcortical Processing
# log_Msg "Subcortical Processing"
bash SubcorticalProcessing.sh "$AtlasSpaceFolder" "$ROIFolder" "$FinalfMRIResolution" "$ResultsFolder" "$NameOffMRI" "$SmoothingFWHM" "$GrayordinatesResolution"

#Generation of Dense Timeseries
# log_Msg "Generation of Dense Timeseries"
bash CreateDenseTimeseries.sh "$AtlasSpaceFolder"/"$DownSampleFolder" "$Subject" "$LowResMesh" "$ResultsFolder"/"$NameOffMRI" "$SmoothingFWHM" "$ROIFolder" "$ResultsFolder"/"$OutputAtlasDenseTimeseries" "$GrayordinatesResolution"