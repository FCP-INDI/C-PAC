# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

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

### Fixed
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1741) that was causing `single_step_resampling` to inadvertently cause unexpected forks in the pipeline past transform application.
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1702) that was causing `single_step_resampling` to crash with `3dVolReg`
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1686) that was causing datasets containing phase-difference field maps to crash.
- Fixed merge error preventing QC files and surface derivatives copying to output directory and renaming connectome → connectivity matrix files
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

[unreleased]: https://github.com/FCP-INDI/C-PAC/compare/v1.8.3...develop
[1.8.3]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.3
[1.8.2]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.2
[1.8.1]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.1
