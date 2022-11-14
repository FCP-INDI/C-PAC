<!-- Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with C-PAC. If not, see <https://www.gnu.org/licenses/>. -->
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Added
- Added the ability to downsample to 10K or 2K resolution for freesurfer runs
- Added the ability to run AFNI 3dDespike on template-space BOLD data.
- Added the ability to ingress TotalReadoutTime from epi field map meta-data from the JSON sidecars.
- Added the ability to use TotalReadoutTime of epi field maps in the calculation of FSL topup distortion correction.
- Difference method (``-``) for ``CPAC.utils.configuration.Configuration`` instances
- Calculate reho and alff when timeseries in template space
- Added new default pipeline that uses FSL-BET for brain extraction. Previous default pipleine is now called default-deprecated.
- Added ``fail_fast`` configuration setting and CLI flag

### Changed
- Space labels in output filenames now contain specific template labels for MNI templates
- Added a level of depth to `working` directories to match `log` and `output` directory structure
- Renamed participant-pipeline-level `output` directory prefix to `pipeline_` to match `log` and `working` paths
- Changed the 1mm atlases chosen in the rbc-options preconfig to the 2mm versions
- For Nilearn-generated correlation matrices, diagonals are now set to all `1`s (were all `0`s)
- Added ability to apply nusiance correction to template-space BOLD images
- Removed ability to run single-step-resampling on motion-corrected BOLD data
- Moved default pipeline config into directory with other preconfigs
- Added crash messages from during and before graph building to logs
- Added data-config-specific hash string to C-PAC-generated config files
- Updated `rbc-options` preconfig to use `fmriprep-options` preprocessing
- Changed minimized pipeline base from `default` to `blank`
- Removed deprecated `--disable_file_logging` CLI flag
- Improved flexibility of some pipeline options regarding the application of distortion correction transforms
- Pinned AFNI to AFNI_21.1.00

### Fixed
- Fixed an issue where the distortion correction un-warps were not being applied to the final template-space BOLD time series data depending on pipeline configuration decisions.
- Fixed [a bug](https://github.com/FCP-INDI/C-PAC/issues/1779) in which generated pipeline configs were not 100% accurate. The only affected configurable option discovered in testing was seed-based correlation analysis always reverting to the default configuration.
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1795) that was causing `cpac run` to fail when passing a manual random seed via `--random_seed`.
- Replaces ``DwellTime`` with ``EffectiveEchoSpacing`` for FSL usage of the term
- Fixed an issue that was causing some epi field maps to not be ingressed if the BIDS tags were not in the correct order.
- Fixed an issue where some phase-difference GRE field map files were not properly being ingressed if the filenames were not expected.
- Fixed a bug where ALFF & f/ALFF would not run if frequency filtering was disabled earlier in the pipeline.

## [v1.8.4] - 2022-06-27

### Added
- Added the ability to follow symlinks for BIDS directories
- Added log of expected outputs, generated at the beginning of the run
- Added additional surface derivatives to outputs directory
- Added additional time series outputs from ABCD-options related processes to outputs directory
- Added the ability to disable the exception raised if the initial resource check estimates more memory is needed
- Added `--runtime_usage` and `--runtime_buffer` flags and related pipeline config entries and functionality
- Added additional error checks and a more informative message for node block connecting in the engine
- Expanded some surface post-processing workflows to be more flexible with other pipeline configurations
- Added [lint definition file](./.pylintrc) (for developers)
- Added list of preconfigured pipelines to the usage string
- Added CI smoke tests for preconfigs

### Changed
- Made ndmg correlation matrices a configurable option
- Made surface output filenames BIDSier
- Uses max instead of sum for intial memory estimation
- Updated rbc-options configuration
- Updated XCP-QC files to better adhere to XCP
- Updated CI to only rebuild software dependencies on change
- Replaced deprecated `optparse.OptionError` with to `click.BadParameter`
- Relicensed C-PAC from BSD-3-Clause to LGPL-3.0-or-later

### Fixed
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1741) that was causing `single_step_resampling` to inadvertently cause unexpected forks in the pipeline past transform application.
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1702) that was causing `single_step_resampling` to crash with `3dVolReg`
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1686) that was causing datasets containing phase-difference field maps to crash.
- Fixed merge error preventing QC files and surface derivatives copying to output directory and renaming connectome → connectivity matrix files
- Fixed some incorrect connections in XCP-QC file generation
- Fixed a bug where subsequent subjects' logs were being appended to prior subjects' logs
- Fixed templates used for rodent pipeline and added outputs that were missing
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1556) in which ITK header imprecision was causing N4 bias correction to crash
- Fixed a bug where a node with zero nanoseconds in timing information causes report generation to fail.
- Fixed a bug in which `--tracking_opt-out` was not applied to calls to `utils`

## [1.8.3] - 2022-02-11

### Added
- Added XCP-style quality control file
- Added RBC-options pipeline preconfiguration
- Added `engine.log` (when verbose debugging is on)
- Added ability to fix random seed for
    - `antsAI`
    - `antsRegistration`
    - `Atropos` (fixed but not specified)
    - `fslmaths`
    - `mri_vol2vol`
    - `recon-all`
- Added ability to use lateral ventricles mask in place of cerebrospinal fluid mask when when segmentation is Off, specifically for the rodent pipeline, but works on any dataset when segmentation is off

### Changed
- In a given pipeline configuration, segmentation probability maps and binary tissue masks are warped to template space, and those warped masks are included in the output directory
    - if `registration_workflows['functional_registration']['EPI_registration']['run segmentation']` is `On` and `segmentation['tissue_segmentation']['Template_Based']['template_for_segmentation']` includes `EPI_Template`
    
      and/or
    - if `registration_workflows['anatomical_registration']['run']` is `On` and `segmentation['tissue_segmentation']['Template_Based']['template_for_segmentation']` includes `T1_Template`
- Renamed connectivity matrices from `*_connectome.tsv` to `*_correlations.tsv`
- Moved some ephemeral logging statements into `pypeline.log`

### Fixed
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1638) in which working connectivity matrix filepaths were generated incorrectly, preventing generating matrices depending on container bindings
- Fixed broken links in README
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1575) in which anatomical-only configurations required functional data directories
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1532) in which nuisance regressors would crash when segmentation is off and no CSF mask is provided

## [1.8.2] - 2021-12-02

### Added 

- Added FSL-TOPUP as an option for distortion correction.
- Added changelog
- Added CHD8 mouse template (`/cpac_templates/chd8_functional_template_noise_mask_ag.nii.gz`) 
- Added commandline flags `--T1w_label` and `--bold_label`
- Added the ability to ingress an entire FreeSurfer output directory to bypass surface analysis if already completed elsewhere
- Added AFNI and Nilearn implementations of Pearson and partial correlation matrices

### Changed

- Expanded meta-data ingress for EPI field maps to include more fields when parsing BIDS sidecar JSONs.
- Updated possible inputs for T2w processing and ACPC-alignment blocks to increase the modularity of these pipeline options.
- `master` branch renamed `main`
- Packaged templates in https://github.com/FCP-INDI/C-PAC_templates

### Deprecated

- `master` branch name (renamed `main`)

### Fixed

- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1620) in which the preprocessed T2w data would not be found in the resource pool when expected.
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1582) in which distortion correction-related field map ingress would raise `IndexError: list index out of range` when ingressing the field maps.
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1572) in which some nodes would raise `KeyError: 'in_file'` when estimating memory allocation
- Improved memory management for multi-core node allocation.
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1548) where `--participant_label [A B C]` would skip first and last labels (`A` and `C`).
- Stripped the ABI tag note (which was preventing the library from loading dynamically on some host operating systems) from `libQt5Core.so.5` in the ABCD-HCP variant image.
- Fixed an issue blocking non-C-PAC output data from being read in without sidecar meta-data.

## [1.8.1] - 2021-09-17

See [Version 1.8.1 Beta](https://fcp-indi.github.io/docs/user/release_notes/v1.8.1) for release notes for v1.8.1 and [Release Notes](https://fcp-indi.github.io/docs/user/release_notes) for all release notes back to v0.1.1.

[unreleased]: https://github.com/FCP-INDI/C-PAC/compare/v1.8.4...develop
[1.8.4]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.4
[1.8.3]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.3
[1.8.2]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.2
[1.8.1]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.1
