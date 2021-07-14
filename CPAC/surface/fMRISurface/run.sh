while [[ "$#" -gt 0 ]]; do
    case $1 in
        --post_freesurfer_folder) Path="$2"; shift ;;
        --subject) Subject="$2"; shift ;;
		--lowresmesh) LowResMeshes="$2"; shift ;;
		--fmrires) FinalfMRIResolution="$2"; shift ;;
		--smoothingFWHM) SmoothingFWHM="$2"; shift ;;
		--grayordinatesres) GrayordinatesResolution="$2"; shift ;;
		--fmri) AtlasSpacefMRI="$2"; shift ;;
		--scout) Scout="$2"; shift ;;
        # *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Path=`opts_GetOpt1 "--post_freesurfer_folder" $@`
# Subject=`opts_GetOpt1 "--subject" $@`
# LowResMesh=`opts_GetOpt1 "--lowresmesh" $@`
# FinalfMRIResolution=`opts_GetOpt1 "--fmrires" $@`
# SmoothingFWHM=`opts_GetOpt1 "--smoothingFWHM" $@`
# GrayordinatesResolution=`opts_GetOpt1 "--grayordinatesres" $@`
# AtlasSpacefMRI=`opts_GetOpt1 "--fmri" $@` # func to standard
# Scout=`opts_GetOpt1 "--scout" $@`

NameOffMRI=task-rest01
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
echo "Make fMRI Ribbon"
echo "mkdir -p ${ResultsFolder}/RibbonVolumeToSurfaceMapping"
mkdir -p "$ResultsFolder"/RibbonVolumeToSurfaceMapping
bash RibbonVolumeToSurfaceMapping.sh "$ResultsFolder"/RibbonVolumeToSurfaceMapping "$ResultsFolder"/"$NameOffMRI" "$Subject" "$AtlasSpaceFolder"/"$DownSampleFolder" "$LowResMesh" "$AtlasSpaceFolder"/"$NativeFolder" "${RegName}"

#Surface Smoothing
echo "Surface Smoothing"
bash SurfaceSmoothing.sh "$ResultsFolder"/"$NameOffMRI" "$Subject" "$AtlasSpaceFolder"/"$DownSampleFolder" "$LowResMesh" "$SmoothingFWHM"

#Subcortical Processing
echo "Subcortical Processing"
bash SubcorticalProcessing.sh "$AtlasSpaceFolder" "$ROIFolder" "$FinalfMRIResolution" "$ResultsFolder" "$NameOffMRI" "$SmoothingFWHM" "$GrayordinatesResolution"

#Generation of Dense Timeseries
echo "Generation of Dense Timeseries"
bash CreateDenseTimeseries.sh "$AtlasSpaceFolder"/"$DownSampleFolder" "$Subject" "$LowResMesh" "$ResultsFolder"/"$NameOffMRI" "$SmoothingFWHM" "$ROIFolder" "$ResultsFolder"/"$OutputAtlasDenseTimeseries" "$GrayordinatesResolution"