<!-- Copyright (C) 2023  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with C-PAC. If not, see <https://www.gnu.org/licenses/>. -->

<!-- Don't rename to '.github/README.md' or it will override the root readme. -->

# C-PAC/.github/README/README.md

The `.github` directory contains Dockerfiles, scripts, requirements files, and workflow configurations for [this repository's GitHub Actions](https://github.com/FCP-INDI/C-PAC/actions).

When updating these files, please update this README as necessary.

[The `workflows` directory](../workflows) contains the configurations for the Actions themselves. The other directories support these configs.

```mermaid
flowchart LR
    subgraph Dockerfiles
      base[base-.*]

      cpacdockerfiles["C-PAC.develop-(?!lite).*"]

      cpaclitedockerfile[C-PAC.develop-lite.*]

      stagedockerfiles[AFNI.*\nANTs.*\nc3d.*\nconnectome-workbench.*\nFSL.*\nICA-AROMA.*\nmsm.*]

      ubuntudockerfiles[Ubuntu.*]
    end
    subgraph scripts
      get_package_id.py
      get_pr_base_shas
      local_ghcr
    end
    subgraph smoke_test_participant.yml
      smoke_test_human
      smoke_test_nhp
      smoke_test_rodent
    end
    subgraph stage_requirements
      .*.txt
    end
    subgraph workflows
      subgraph build_C-PAC.yml
        bCPAC[[C-PAC]]
      end
      subgraph build_and_test.yml
        ubuntu[[Ubnutu]]-->stages[[stages]]-->build-base[[build-base]]-->build-base-standard[[build-base-standard]]
        
        Circle_tests[[Circle_tests]]
      
        build-base-standard-->C-PAC
        C-PAC[[C-PAC]]-->bCPAC
        C-PAC-->Circle_tests
        C-PAC-->smoke-tests-participant
        C-PAC-->C-PAC-lite

        build-base-->C-PAC-lite
        C-PAC-lite[[C-PAC-lite]]-->bCPAC
        C-PAC-lite-->Circle_tests
        C-PAC-lite-->smoke-tests-participant

        build-base-->C-PAC-ABCD-HCP
        C-PAC-ABCD-HCP[[C-PAC-ABCD-HCP]]-->bCPAC
        C-PAC-ABCD-HCP-->Circle_tests
        C-PAC-ABCD-HCP-->smoke-tests-participant

        build-base-->C-PAC-fMRIPREP
        C-PAC-fMRIPrep[[C-PAC-fMRIPrep-LTS]]-->bCPAC
        C-PAC-fMRIPrep-->Circle_tests
        C-PAC-fMRIPrep-->smoke-tests-participant

        smoke-tests-participant[[smoke-tests-participant]]
      end

      on_push.yml-->build_and_test.yml

      delete_images.yml
    end

    subgraph yaml_template[CPAC/utils/configuration/yaml_template.py]

      update_all_preconfigs[[update_all_preconfigs]]
    end

    build-base-standard<-->base<-->build-base

    Circle_tests-->CircleCI((Run tests on Circle CI))

    on_push.yml<-->get_pr_base_shas
    on_push.yml-->update_all_preconfigs

    cpacdockerfiles<-->C-PAC
    cpacdockerfiles<-->C-PAC-ABCD-HCP
    cpacdockerfiles<-->C-PAC-fMRIPrep

    cpaclitedockerfile<-->C-PAC-lite

    delete>delete branch]-->delete_images.yml

    delete_images.yml<-->get_package_id.py

    build-base<-->local_ghcr
    bCPAC<-->local_ghcr
    stages<-->local_ghcr

    push>git push]-->on_push.yml

    smoke-tests-participant-->smoke_test_human
    smoke-tests-participant-->smoke_test_nhp
    smoke-tests-participant-->smoke_test_rodent

    stagedockerfiles<-->stages

    build-base<-->.*.txt

    ubuntudockerfiles<-->ubuntu
```

In this [Mermaid flowchart](https://mermaid.js.org/syntax/flowchart.html), these shapes are used:

concept | shape
---|---
directory | `subgraph`
file | rectangle (`[]`) or `subgraph`
workflow job | subprocess (`[[]]`)
trigger action | asymetric (`>]`)
external API | circle (`(())`)

## Multistage builds

C-PAC images are built in stages like

```mermaid
flowchart LR

OS["operating system"] --> deps["versioned dependencies"]
deps --> base["base image"]
base --> CPAC["C-PAC version"]
```

See [fcp-indi.github.io/docs/nightly/developer/installation](https://fcp-indi.github.io/docs/nightly/developer/installation) for more information about C-PAC's use of multistage builds.

GitHub Actions will dynamically skip building stages based on when a stage's Dockerfile or a workflow file was last modified. To force a stage to rebuild, include

```BASH
[rebuild ${DOCKERFILE_BASENAME}]
```

in a commit message and push. For example, to rebuild `ghcr.io/fcp-indi/c-pac/ubuntu:bionic-non-free`, include

```BASH
[rebuild Ubuntu.bionic-non-free]
```

in the commit message. For this to work, all of these must be true:

1. The Dockerfile is listed in a

   ```YAML
   strategy:
     matrix:
       Dockerfile: 
   ```

   in a job in a [workflow](../workflows) file.
2. The Dockerfile exists.
3. The commit message including

   ```BASH
   [rebuild ${DOCKERFILE_BASENAME}]
   ```

   is for the most recent commit when pushing to GitHub.

## Dockerfile tips

### COPY

For multistage builds, a trailing slash is necessary for Docker to treat a path as a directory, so `COPY` commands for directories should look like

```Dockerfile
COPY --from=STAGE /src/path/ /dest/path/
```

and for files should look like

```Dockerfile
COPY --from=STAGE /src/path /dest/path
```
