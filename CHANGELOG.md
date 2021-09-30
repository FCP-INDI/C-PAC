# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

- Added changelog

### Changed

- `master` branch renamed `main`

### Deprecated

- `master` branch name (renamed `main`)

### Fixed

- Bug (https://github.com/FCP-INDI/C-PAC/issues/1582) in which blip-up/down distortion correction and EPI field map ingress would raise `IndexError: list index out of range` when ingressing the field maps.
- Bug (https://github.com/FCP-INDI/C-PAC/issues/1572) in which some nodes would raise `KeyError: 'in_file'` when estimating memory allocation

## [1.8.1] - 2021-09-17

See [Version 1.8.1 Beta](https://fcp-indi.github.io/docs/user/release_notes/v1.8.1) for release notes for v1.8.1 and [Release Notes](https://fcp-indi.github.io/docs/user/release_notes) for all release notes back to v0.1.1.

[unreleased]: https://github.com/FCP-INDI/C-PAC/compare/v1.8.1...develop
[1.8.1]: https://github.com/FCP-INDI/C-PAC/releases/tag/v1.8.1
