<!-- Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with C-PAC. If not, see <https://www.gnu.org/licenses/>. -->
## Software dependencies and variant images
We currently have one main and 3 variant images:
* `ABCD-HCP`: dependency versions matched to [ABCD-HCP BIDS fMRI Pipeline](https://github.com/DCAN-Labs/abcd-hcp-pipeline/releases/tag/v0.1.1) versions
* `fMRIPrep-LTS`: dependency versions matched to [fMRIPrep Long-term support](https://reproducibility.stanford.edu/fmriprep-lts/) versions
* `lite`: same dependency versions as main image without FreeSurfer (smaller image)

To save time building Docker images, our continuous integration is set up to build and [store staging images](https://github.com/FCP-INDI?tab=packages&repo_name=C-PAC) independent of changes to C-PAC itself. Images are rebuilt when their Dockerfile changes on a Git branch.

All Dockerfiles are stored in [`.github/Dockerfiles`](./.github/Dockerfiles).

### Naming conventions

We have 3 types of staging Dockerfiles: operating system, software dependency, and C-PAC variant.

#### operating system

`{OS name}.{version}.Dockerfile`

#### software dependency

`{software name}.{version}-{OS version}.Dockerfile`

#### C-PAC variant

`C-PAC.{version}[-{variant}]-{OS version}.Dockerfile` (`-{variant}` is omitted for the main image)

### Changing dependencies

* To change a dependency in a C-PAC image, update the stage images at the top of the relevant `.github/Dockerfiles/C-PAC.develop-*.Dockerfile`.
* If a Dockerfile does not yet exist for the added dependency, create a Dockerfile for the new dependency and add the filename (without extension) to [`jobs.stages.strategy.matrix.Dockerfile` in `.github/workflows/build_stages.yml`](https://github.com/FCP-INDI/C-PAC/blob/4e18916384e52c3dc9610aea3eed537c19d480e3/.github/workflows/build_stages.yml#L77-L97)
* If no Dockerfiles use the removed dependency, remove the Dockerfile for the dependency and remove the filename from [`jobs.stages.strategy.matrix.Dockerfile` in `.github/workflows/build_stages.yml`](https://github.com/FCP-INDI/C-PAC/blob/4e18916384e52c3dc9610aea3eed537c19d480e3/.github/workflows/build_stages.yml#L77-L97)
