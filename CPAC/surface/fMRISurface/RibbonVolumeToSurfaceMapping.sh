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

# Source: https://github.com/DCAN-Labs/DCAN-HCP/blob/32254dc/fMRISurface/scripts/RibbonVolumeToSurfaceMapping.sh

# Modifications copyright (C) 2021 - 2024  C-PAC Developers
# This file is part of C-PAC.

set -e
echo -e "\n START: RibbonVolumeToSurfaceMapping"

WorkingDirectory="$1"
VolumefMRI="$2"
Subject="$3"
DownsampleFolder="$4"
LowResMesh="$5"
AtlasSpaceNativeFolder="$6"
RegName="$7"

if [ ${RegName} = "FS" ]; then
    RegName="reg.reg_LR"
fi

NeighborhoodSmoothing="5"
Factor="0.5"


LeftGreyRibbonValue="1"
RightGreyRibbonValue="1"

for Hemisphere in L R ; do
  if [ $Hemisphere = "L" ] ; then
    GreyRibbonValue="$LeftGreyRibbonValue"
  elif [ $Hemisphere = "R" ] ; then
    GreyRibbonValue="$RightGreyRibbonValue"
  fi
  wb_command -create-signed-distance-volume "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii "$VolumefMRI"_SBRef.nii.gz "$WorkingDirectory"/"$Subject"."$Hemisphere".white.native.nii.gz
  wb_command -create-signed-distance-volume "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii "$VolumefMRI"_SBRef.nii.gz "$WorkingDirectory"/"$Subject"."$Hemisphere".pial.native.nii.gz
  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".white.native.nii.gz -thr 0 -bin -mul 255 "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz
  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz -bin "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz
  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".pial.native.nii.gz -uthr 0 -abs -bin -mul 255 "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz
  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz -bin "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz
  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz -mas "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz -mul 255 "$WorkingDirectory"/"$Subject"."$Hemisphere".ribbon.nii.gz
  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".ribbon.nii.gz -bin -mul $GreyRibbonValue "$WorkingDirectory"/"$Subject"."$Hemisphere".ribbon.nii.gz
  rm "$WorkingDirectory"/"$Subject"."$Hemisphere".white.native.nii.gz "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz "$WorkingDirectory"/"$Subject"."$Hemisphere".pial.native.nii.gz "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz
done

fslmaths "$WorkingDirectory"/"$Subject".L.ribbon.nii.gz -add "$WorkingDirectory"/"$Subject".R.ribbon.nii.gz "$WorkingDirectory"/ribbon_only.nii.gz
rm "$WorkingDirectory"/"$Subject".L.ribbon.nii.gz "$WorkingDirectory"/"$Subject".R.ribbon.nii.gz

fslmaths "$VolumefMRI" -Tmean "$WorkingDirectory"/mean -odt float
fslmaths "$VolumefMRI" -Tstd "$WorkingDirectory"/std -odt float
fslmaths "$WorkingDirectory"/std -div "$WorkingDirectory"/mean "$WorkingDirectory"/cov

fslmaths "$WorkingDirectory"/cov -mas "$WorkingDirectory"/ribbon_only.nii.gz "$WorkingDirectory"/cov_ribbon

fslmaths "$WorkingDirectory"/cov_ribbon -div `fslstats "$WorkingDirectory"/cov_ribbon -M` "$WorkingDirectory"/cov_ribbon_norm
fslmaths "$WorkingDirectory"/cov_ribbon_norm -bin -s $NeighborhoodSmoothing "$WorkingDirectory"/SmoothNorm
fslmaths "$WorkingDirectory"/cov_ribbon_norm -s $NeighborhoodSmoothing -div "$WorkingDirectory"/SmoothNorm -dilD "$WorkingDirectory"/cov_ribbon_norm_s$NeighborhoodSmoothing
fslmaths "$WorkingDirectory"/cov -div `fslstats "$WorkingDirectory"/cov_ribbon -M` -div "$WorkingDirectory"/cov_ribbon_norm_s$NeighborhoodSmoothing "$WorkingDirectory"/cov_norm_modulate
fslmaths "$WorkingDirectory"/cov_norm_modulate -mas "$WorkingDirectory"/ribbon_only.nii.gz "$WorkingDirectory"/cov_norm_modulate_ribbon

STD=`fslstats "$WorkingDirectory"/cov_norm_modulate_ribbon -S`
echo $STD
MEAN=`fslstats "$WorkingDirectory"/cov_norm_modulate_ribbon -M`
echo $MEAN
Lower=`echo "$MEAN - ($STD * $Factor)" | bc -l`
echo $Lower
Upper=`echo "$MEAN + ($STD * $Factor)" | bc -l`
echo $Upper

fslmaths "$WorkingDirectory"/mean -bin "$WorkingDirectory"/mask
fslmaths "$WorkingDirectory"/cov_norm_modulate -thr $Upper -bin -sub "$WorkingDirectory"/mask -mul -1 "$WorkingDirectory"/goodvoxels

for Hemisphere in L R ; do
  for Map in mean cov ; do
    wb_command -volume-to-surface-mapping "$WorkingDirectory"/"$Map".nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii -ribbon-constrained "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii -volume-roi "$WorkingDirectory"/goodvoxels.nii.gz
    wb_command -metric-dilate "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii -nearest
    wb_command -metric-mask "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii
    wb_command -volume-to-surface-mapping "$WorkingDirectory"/"$Map".nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii -ribbon-constrained "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii
    wb_command -metric-mask "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii
    wb_command -metric-resample "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WorkingDirectory"/"$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.func.gii -area-surfs "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    wb_command -metric-mask "$WorkingDirectory"/"$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.func.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$WorkingDirectory"/"$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.func.gii
    wb_command -metric-resample "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WorkingDirectory"/"$Hemisphere"."$Map"_all."$LowResMesh"k_fs_LR.func.gii -area-surfs "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    wb_command -metric-mask "$WorkingDirectory"/"$Hemisphere"."$Map"_all."$LowResMesh"k_fs_LR.func.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$WorkingDirectory"/"$Hemisphere"."$Map"_all."$LowResMesh"k_fs_LR.func.gii
  done
  wb_command -volume-to-surface-mapping "$WorkingDirectory"/goodvoxels.nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii -ribbon-constrained "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii
  wb_command -metric-mask "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii
  wb_command -metric-resample "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WorkingDirectory"/"$Hemisphere".goodvoxels."$LowResMesh"k_fs_LR.func.gii -area-surfs "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
  wb_command -metric-mask "$WorkingDirectory"/"$Hemisphere".goodvoxels."$LowResMesh"k_fs_LR.func.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$WorkingDirectory"/"$Hemisphere".goodvoxels."$LowResMesh"k_fs_LR.func.gii

  wb_command -volume-to-surface-mapping "$VolumefMRI".nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$VolumefMRI"."$Hemisphere".native.func.gii -ribbon-constrained "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii -volume-roi "$WorkingDirectory"/goodvoxels.nii.gz
  wb_command -metric-dilate "$VolumefMRI"."$Hemisphere".native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$VolumefMRI"."$Hemisphere".native.func.gii -nearest
  wb_command -metric-mask  "$VolumefMRI"."$Hemisphere".native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii  "$VolumefMRI"."$Hemisphere".native.func.gii
  wb_command -metric-resample "$VolumefMRI"."$Hemisphere".native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$VolumefMRI"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii -area-surfs "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
  wb_command -metric-mask "$VolumefMRI"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$VolumefMRI"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii
done

echo " END: RibbonVolumeToSurfaceMapping"
