This directory contains Dockerfiles, scripts, requirements files, and workflow configurations for [this repository's GitHub Actions](https://github.com/FCP-INDI/C-PAC/actions).

When updating these files, please update this README as necessary.

[The `workflows` directory](./workflows) contains the configurations for the Actions themselves. The other directories support these configs.

```mermaid
flowchart TD
    push([git push])-->check_updated_preconfigs.yml
    delete([delete])-->delete_images.yml
    get_package_id.py-->delete_images.yml
    ubuntudockerfiles-->ubuntu
    stagedockerfiles-->stages
    base-->build-base
    local_ghcr-->stages
    local_ghcr-->build-base
    local_ghcr-->build_C-PAC.yml
    cpacdockerfiles-->C-PAC
    cpaclitedockerfile-->C-PAC-lite
    build_stages.yml-->build_and_test.yml
    get_pr_base_shas-->check_updated_preconfigs.yml
    subgraph Dockerfiles
      ubuntudockerfiles[Ubuntu.*]
      stagedockerfiles[AFNI.*\nANTs.*\nc3d.*\nconnectome-workbench.*\nFSL.*\nICA-AROMA.*\nmsm.*]
      base[base-.*]
      cpacdockerfiles["C-PAC.develop-(?!lite).*"]
      cpaclitedockerfile[C-PAC.develop-lite.*]
    end
    subgraph scripts
      get_package_id.py
      get_pr_base_shas
      local_ghcr
    end
    stage_requirements-->build-base
    subgraph workflows
      subgraph build_stages.yml
        ubuntu([Ubnutu])-->stages([stages])-->build-base([build-base])
      end
      subgraph build_C-PAC.yml
      end
      subgraph build_and_test.yml
        C-PAC-->build_C-PAC.yml
        C-PAC-lite-->build_C-PAC.yml
        C-PAC-ABCD-HCP-->build_C-PAC.yml
        C-PAC-fMRIPrep-->build_C-PAC.yml
        smoke-tests-participant
      end
    build_C-PAC.yml-->smoke-tests-participant
    smoke-tests-participant-->smoke_test_participant.yml-->Circle_tests((Run tests on Circle CI))
    check_updated_preconfigs.yml-->build_stages.yml
    subgraph smoke_test_participant.yml
      smoke_test_human
      smoke_test_nhp
      smoke_test_rodent
    end
    delete_images.yml
    end
```