<!-- Copyright (C) 2022-2023  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with C-PAC. If not, see <https://www.gnu.org/licenses/>. -->
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.8.7.post1] - unreleased

### Changed

- Disabled variant image builds.

### Fixed

- Graph-building bugs that prevented longitudinal workflows from running.

## [1.8.7] - 2024-05-03

### Added

- `Robustfov` feature in `FSL-BET` to crop images ensuring removal of neck regions that may appear in the skull-stripped images. 
- Ability to throttle nodes, estimating all available memory when threading.
- Ability to configure FreeSurfer ingress from the command line.

### Changed

- The ABCD-pipeline based surface post-processing workflows have been modularized to be more robust, resolving a running issue with this part of the pipeline stalling or crashing in some runs.
- Moved autoversioning from CI to pre-commit
- Updated `FSL-BET` config to default `-mask-boolean` flag as on, and removed all removed `mask-boolean` keys from configs.
- Added `dvars` as optional output in `cpac_outputs`.

### Fixed

- Fixed a bug where ingressing fmriprep outputs into C-PAC with a blank nuisance confounds field in the C-PAC pipeline configuration file would cause a crash.
- Fixed a bug where spatial smoothing and z-scoring of final outputs would sometimes fail to run when running a C-PAC pipeline that would ingress fmriprep outputs.
- Fixed a bug where ingress of distortion correction-related field map metadata would sometimes fail to recognize both echo times, when there were two present, leading to an error message claiming an echo time is missing.
- Changed an extraneous default pipeline configuration setting - `surface_connectivity` is now disabled in the default configuration as intended.

## [1.8.6] - 2024-01-15

### Added

- Some automatic handling of user-provided BIDSy atlas names.
- `sig_imports` static method decorator for `Function` nodes, to accommodate type hinting in signatures of `Function` node functions.
- Ability to ingress an fmriprep output directory and run derivatives with C-PAC
- fmriprep-ingress preconfig that runs derivatives
- String representations for NodeBlock and ResourcePool class instances
- `switch_is_off`, `switch_is_on` and `switch_is_on_off` methods to `Configuration` class
- `__repr__` and `__str__` methods to `ResourcePool`s and `NodeBlockFunction`s

### Fixed

- Fixed a bug where some connectivity matrices wouldn't generate if anatomical and functional outputs were in different resolutions.
- Handling of `3dECM` outputs for AFNI ≥ 21.1.1.
- Fixed a bug where sparsity thresholds were not being scaled for network centrality.
- Fixed a bug where `calculate_motion_first` would not calculate motion at all.
- Fixed a bug in parsing `FROM: /file/path` syntax

### Changed

- Updates develop version numbers to adhere to [PEP440](https://peps.python.org/pep-0440) by changing `{major}.{minor}.{patch}.{SHA}-{dev}` to `{major}.{minor}.{patch}.dev{integer}+{SHA}`
- Adds checksum steps to `curl`d steps in Docker build process (for standard and `lite` images)
- Makes in-container root directory writable by all
- Updates multiword CLI commands and options to accept either standard `-`s or backwards-compatible `_`s interchangeably
- Disabled `--use-estimate-learning-rate-once` in `antsRegistration` (ANTsX/ANTs#1405; ANTsX/ANTs#1411)
- Removes `torch` from preinstalled dependencies and only installs if we're running `unet`
- Uses highest resolution available locally as reference when resampling a template to a non-packaged resolution (was always using 1mm reference before)
- Updates config boolean validation from anything-truthy-is-True (e.g., `[True, False]`, or `[False]`, or a typo like `Offf`) to only accepting bools, ints, and YAML boolean strings like "On" and "Off" as boolean
- When applying a filter to motion parameters, now C-PAC reports both the original and the filtered motion parameters and uses the original parameters for qc. Previous versions only reported the filtered parameters and used the filtered parameters for qc.
- Makes nuisance regression space non-forkable. In v1.8.5, nuisance regression forking was broken, so this change should not cause backwards-compatibility issues.

### Added dependencies

- `click-aliases`
- `dc`
- `semver`

### Removed dependencies

- `bids-validator`
- `MSM` (now packaged with `FSL`)
- `Node`
- `NVM`
- `simplejson`
- `wxpython`
- `yamlordereddictloader`

### Upgraded dependencies

- `AFNI` 21.1.00 'Domitian' → 23.3.09 'Septimius Severus'
- `ANTs` 2.3.3 'Leptomyrmex' → 2.4.3 'Emplastus'
- `boto3` 1.7.37 → 1.28.4
- `click` 6.7 → 8.1.5
- `configparser` 3.7.4 → 5.3.0
- `FSL` 5.0.10 (5.0.9) → 6.0.6.5
- `future` 0.16.0 → 0.18.3
- `ICA-AROMA` 0.4.3-beta → 0.4.4-beta
- `INDI-Tools` 0.0.7@`main` → 0.0.7@`998c246`
- `joblib` 1.0.1 → 1.2.0
- `matplotlib` 3.1.3 → 3.7.1
- `networkx` 2.4 → 3.1
- `nibabel` 3.0.1 → 5.1.0
- `nilearn` 0.4.1 → 0.10.0
- `nipype` 1.5.1 → 1.8.6
- `numpy` 1.21.0 → 1.25.1
- `pandas` 1.0.5 → 2.0.3
- `patsy` 0.5.0 → 0.5.3
- `prov` 1.5.2 → 2.0.0
- `psutil` 5.6.6 → 5.9.5
- `pybids` 0.15.1 → 0.15.6
- `PyPEER` 1.0@`main` → 1.1@`6965d2b`
- `Python` 3.7.13 → 3.10.6 (holding back from 3.11 for `sdcflows`, `torch`, and `torchvision`)
- `python-dateutil` 2.7.3 → 2.8.2
- `PyYAML` 5.4 → 6.0
- `requests` 2.21.0 → 2.28.2
- `scipy` 1.6.3 → 1.11.1
- `sdcflows` 2.0.5 → 2.4.0
- `torch` 1.2.0+cu92 → 1.13.1
- `torchvision` 0.4.0+cu92 → 0.14.1
- `traits` 4.6.0 → 6.3.2
- `Ubuntu` 18.04 'Bionic Beaver' → 22.04 'Jammy Jellyfish'
- `voluptuous` 0.12.0 → 0.13.1
- `wb_command` neurodebian latest → 1.5.0

## [1.8.5] - 2023-05-24

### Added

- Added the ability to downsample to 10K or 2K resolution for freesurfer runs
- Added the ability to run AFNI 3dDespike on template-space BOLD data.
- Added the ability to ingress TotalReadoutTime from epi field map meta-data from the JSON sidecars.
- Added the ability to use TotalReadoutTime of epi field maps in the calculation of FSL topup distortion correction.
- Difference method (`-`) for `CPAC.utils.configuration.Configuration` instances
- Calculate reho and alff when timeseries in template space
- Added new default pipeline that uses FSL-BET for brain extraction. Previous default pipleine is now called default-deprecated.
- Added `fail_fast` configuration setting and CLI flag
- Added abililty to fork on motion filter
- Added [`sdcflows`](https://www.nipreps.org/sdcflows/2.0/) to CPAC requirements
- Added NodeBlock information to `pypeline.log` when verbose debugging is on
- Added the ability to ingress FreeSurfer data into CPAC
- Added the ability to toggle FreeSurfer derived masks for brain extraction
- Added an optional volume center to FD-J calculation
- Added new preconfig `abcd-prep`, which performs minimal preprocessing on the T1w data in preparation for Freesurfer Recon-All

### Changed

- Freesurfer output directory ingress moved to the data configuration YAML
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
- Updated some output filenaming conventions for human-readability and to move closer to BIDS-derivatives compliance
- Changed motion filter from single dictionary to list of dictionaries
- Changed CI logic to allow non-release tags

### Upgraded dependencies

- `nibabel` 2.3.3 → 3.0.1
- `numpy` 1.16.4 → 1.21.0
- `pandas` 0.23.4 → 1.0.5
- `pybids` 0.13.2 → 0.15.1
- `scipy` 1.4.1 → 1.6.0

### Fixed

- Fixed an issue where the distortion correction un-warps were not being applied to the final template-space BOLD time series data depending on pipeline configuration decisions.
- Fixed [a bug](https://github.com/FCP-INDI/C-PAC/issues/1779) in which generated pipeline configs were not 100% accurate. The only affected configurable option discovered in testing was seed-based correlation analysis always reverting to the default configuration.
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1795) that was causing `cpac run` to fail when passing a manual random seed via `--random_seed`.
- Replaces `DwellTime` with `EffectiveEchoSpacing` for FSL usage of the term
- Fixed an issue that was causing some epi field maps to not be ingressed if the BIDS tags were not in the correct order.
- Fixed an issue where some phase-difference GRE field map files were not properly being ingressed if the filenames were not expected.
- Fixed a bug where ALFF & f/ALFF would not run if frequency filtering was disabled earlier in the pipeline.
- Fixed a bug where `surface_analysis.freesurfer.freesurfer_dir` in the pipeline config was not ingressed at runtime.
- Added public read access to some overly restricted packaged templates
- Fixed a bug where notch filter was always assuming the sampling frequency was `2.0`.

## [1.8.4] - 2022-06-27

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
- Packaged templates in [:octocat:/FCP-INDI/C-PAC_templates](https://github.com/FCP-INDI/C-PAC_templates)

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

[1.8.7.post1]: https://github.com/FCP-INDI/C-PAC/compare/v1.8.7...v1.8.7.post1.dev3
[1.8.7]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.7
[1.8.6]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.6
[1.8.5]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.5
[1.8.4]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.4
[1.8.3]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.3
[1.8.2]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.2
[1.8.1]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.1
