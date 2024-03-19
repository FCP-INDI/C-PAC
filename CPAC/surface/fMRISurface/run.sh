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

# Source: https://github.com/DCAN-Labs/DCAN-HCP/blob/32254dc/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh

# Modifications copyright (C) 2021 - 2024  C-PAC Developers
# This file is part of C-PAC.

# Input Variables
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --post_freesurfer_folder) Path="$2"; shift ;;
        --subject) Subject="$2"; shift ;;
		--lowresmesh) LowResMesh="$2"; shift ;;
		--fmrires) FinalfMRIResolution="$2"; shift ;;
		--smoothingFWHM) SmoothingFWHM="$2"; shift ;;
		--grayordinatesres) GrayordinatesResolution="$2"; shift ;;
		--fmri) AtlasSpacefMRI="$2"; shift ;;
		--scout) Scout="$2"; shift ;;
    esac
    shift
done

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
ResultsFolder="$AtlasSpaceFolder"/"$ResultsFolder"/"$NameOffMRI"
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
bash /code/CPAC/surface/fMRISurface/RibbonVolumeToSurfaceMapping.sh "$ResultsFolder"/RibbonVolumeToSurfaceMapping "$ResultsFolder"/"$NameOffMRI" "$Subject" "$AtlasSpaceFolder"/"$DownSampleFolder" "$LowResMesh" "$AtlasSpaceFolder"/"$NativeFolder" "${RegName}"

#Surface Smoothing
echo "Surface Smoothing"
bash /code/CPAC/surface/fMRISurface/SurfaceSmoothing.sh "$ResultsFolder"/"$NameOffMRI" "$Subject" "$AtlasSpaceFolder"/"$DownSampleFolder" "$LowResMesh" "$SmoothingFWHM"

#Subcortical Processing
echo "Subcortical Processing"
bash /code/CPAC/surface/fMRISurface/SubcorticalProcessing.sh "$AtlasSpaceFolder" "$ROIFolder" "$FinalfMRIResolution" "$ResultsFolder" "$NameOffMRI" "$SmoothingFWHM" "$GrayordinatesResolution"

#Generation of Dense Timeseries
echo "Generation of Dense Timeseries"
bash /code/CPAC/surface/fMRISurface/CreateDenseTimeseries.sh "$AtlasSpaceFolder"/"$DownSampleFolder" "$Subject" "$LowResMesh" "$ResultsFolder"/"$NameOffMRI" "$SmoothingFWHM" "$ROIFolder" "$ResultsFolder"/"$OutputAtlasDenseTimeseries" "$GrayordinatesResolution"
