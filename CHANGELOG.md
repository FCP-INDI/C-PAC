# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Added
- Added XCP-style quality control file
- Added `engine.log` (when verbose debugging is on)

### Changed
- Renamed connectivity matrices from `*_connectome.tsv` to `*_correlations.tsv`

### Fixed
- Fixed [bug](https://github.com/FCP-INDI/C-PAC/issues/1638) in which working connectivity matrix filepaths were generated incorrectly, preventing generating matrices depending on container bindings

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

[unreleased]: https://github.com/FCP-INDI/C-PAC/compare/v1.8.2...develop
[1.8.2]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.2
[1.8.1]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.1
