
# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
import ast
import copy
import hashlib
from itertools import chain
import json
import logging
import os
import re
from typing import Optional
import warnings
import pandas as pd


from nipype import config
from nipype.interfaces.utility import Rename

from CPAC.image_utils.spatial_smoothing import spatial_smoothing
from CPAC.image_utils.statistical_transforms import (
    fisher_z_score_standardize,
    z_score_standardize,
)
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.check_outputs import ExpectedOutputs
from CPAC.pipeline.nodeblock import NodeBlockFunction
from CPAC.pipeline.utils import MOVEMENT_FILTER_KEYS, name_fork, source_set
from CPAC.registration.registration import transform_derivative
from CPAC.resources.templates.lookup_table import lookup_identifier
from CPAC.utils.bids_utils import res_in_filename
from CPAC.utils.configuration import Configuration
from CPAC.utils.datasource import (
    create_anat_datasource,
    create_func_datasource,
    create_general_datasource,
    ingress_func_metadata,
    resolve_resolution,
)
from CPAC.utils.interfaces.datasink import DataSink
from CPAC.utils.interfaces.function import Function
from CPAC.utils.monitoring import (
    getLogger,
    LOGTAIL,
    WARNING_FREESURFER_OFF_WITH_DATA,
    WFLOGGER,
)
from CPAC.utils.outputs import Outputs
from CPAC.utils.utils import (
    check_prov_for_regtool,
    create_id_string,
    get_last_prov_entry,
    read_json,
    write_output_json,
)

from typing import Optional
from .resource import Resource

class ResourcePool:
    def __init__(self, rpool=None, name=None, cfg=None, pipe_list=None):
        if not rpool:
            self.rpool = {}
        else:
            self.rpool = rpool

        if not pipe_list:
            self.pipe_list = []
        else:
            self.pipe_list = pipe_list

        self.name = name
        self.info = {}

        if cfg:
            self.cfg = cfg
            self.logdir = cfg.pipeline_setup["log_directory"]["path"]

            self.num_cpus = cfg.pipeline_setup["system_config"][
                "max_cores_per_participant"
            ]
            self.num_ants_cores = cfg.pipeline_setup["system_config"][
                "num_ants_threads"
            ]

            self.ants_interp = cfg.registration_workflows["functional_registration"][
                "func_registration_to_template"
            ]["ANTs_pipelines"]["interpolation"]
            self.fsl_interp = cfg.registration_workflows["functional_registration"][
                "func_registration_to_template"
            ]["FNIRT_pipelines"]["interpolation"]

            self.func_reg = cfg.registration_workflows["functional_registration"][
                "func_registration_to_template"
            ]["run"]

            self.run_smoothing = (
                "smoothed" in cfg.post_processing["spatial_smoothing"]["output"]
            )
            self.smoothing_bool = cfg.post_processing["spatial_smoothing"]["run"]
            self.run_zscoring = "z-scored" in cfg.post_processing["z-scoring"]["output"]
            self.zscoring_bool = cfg.post_processing["z-scoring"]["run"]
            self.fwhm = cfg.post_processing["spatial_smoothing"]["fwhm"]
            self.smooth_opts = cfg.post_processing["spatial_smoothing"][
                "smoothing_method"
            ]

        self.xfm = [
            "alff",
            "desc-sm_alff",
            "desc-zstd_alff",
            "desc-sm-zstd_alff",
            "falff",
            "desc-sm_falff",
            "desc-zstd_falff",
            "desc-sm-zstd_falff",
            "reho",
            "desc-sm_reho",
            "desc-zstd_reho",
            "desc-sm-zstd_reho",
        ]

    def __repr__(self) -> str:
        params = [
            f"{param}={getattr(self, param)}"
            for param in ["rpool", "name", "cfg", "pipe_list"]
            if getattr(self, param, None) is not None
        ]
        return f'ResourcePool({", ".join(params)})'

    def __str__(self) -> str:
        if self.name:
            return f"ResourcePool({self.name}): {list(self.rpool)}"
        return f"ResourcePool: {list(self.rpool)}"

    def append_name(self, name):
        self.name.append(name)

    def back_propogate_template_name(
        self, wf, resource_idx: str, json_info: dict, id_string: "pe.Node"
    ) -> None:
        """Find and apply the template name from a resource's provenance.

        Parameters
        ----------
        resource_idx : str

        json_info : dict

        id_string : pe.Node

        Returns
        -------
        None
        """
        if "template" in resource_idx and self.check_rpool("derivatives-dir"):
            if self.check_rpool("template"):
                node, out = self.get_data("template")
                wf.connect(node, out, id_string, "template_desc")
        elif "Template" in json_info:
            id_string.inputs.template_desc = json_info["Template"]
        elif (
            "template" in resource_idx and len(json_info.get("CpacProvenance", [])) > 1
        ):
            for resource in source_set(json_info["CpacProvenance"]):
                source, value = resource.split(":", 1)
                if value.startswith("template_") and source != "FSL-AFNI-bold-ref":
                    # 'FSL-AFNI-bold-ref' is currently allowed to be in
                    # a different space, so don't use it as the space for
                    # descendents
                    try:
                        anscestor_json = next(iter(self.rpool.get(source).items()))[
                            1
                        ].get("json", {})
                        if "Description" in anscestor_json:
                            id_string.inputs.template_desc = anscestor_json[
                                "Description"
                            ]
                            return
                    except (IndexError, KeyError):
                        pass
        return

    def get_name(self):
        return self.name

    def check_rpool(self, resource):
        if not isinstance(resource, list):
            resource = [resource]
        for name in resource:
            if name in self.rpool:
                return True
        return False

    def get_pipe_number(self, pipe_idx):
        return self.pipe_list.index(pipe_idx)

    def get_pool_info(self):
        return self.info

    def set_pool_info(self, info_dct):
        self.info.update(info_dct)

    def get_entire_rpool(self):
        return self.rpool

    def get_resources(self):
        return self.rpool.keys()

    def copy_rpool(self):
        return ResourcePool(
            rpool=copy.deepcopy(self.get_entire_rpool()),
            name=self.name,
            cfg=self.cfg,
            pipe_list=copy.deepcopy(self.pipe_list),
        )

    @staticmethod
    def get_raw_label(resource: str) -> str:
        """Remove ``desc-*`` label."""
        for tag in resource.split("_"):
            if "desc-" in tag:
                resource = resource.replace(f"{tag}_", "")
                break
        return resource

    def get_strat_info(self, prov, label=None, logdir=None):
        strat_info = {}
        for entry in prov:
            if isinstance(entry, list):
                strat_info[entry[-1].split(":")[0]] = entry
            elif isinstance(entry, str):
                strat_info[entry.split(":")[0]] = entry.split(":")[1]
        if label:
            if not logdir:
                logdir = self.logdir
            WFLOGGER.info(
                "\n\nPrinting out strategy info for %s in %s\n", label, logdir
            )
            write_output_json(
                strat_info, f"{label}_strat_info", indent=4, basedir=logdir
            )

    def set_json_info(self, resource, pipe_idx, key, val):
        # TODO: actually should probably be able to inititialize resource/pipe_idx
        if pipe_idx not in self.rpool[resource]:
            msg = (
                "\n[!] DEV: The pipeline/strat ID does not exist "
                f"in the resource pool.\nResource: {resource}"
                f"Pipe idx: {pipe_idx}\nKey: {key}\nVal: {val}\n"
            )
            raise Exception(msg)
        if "json" not in self.rpool[resource][pipe_idx]:
            self.rpool[resource][pipe_idx]["json"] = {}
        self.rpool[resource][pipe_idx]["json"][key] = val

    def get_json_info(self, resource, pipe_idx, key):
        # TODO: key checks
        if not pipe_idx:
            for pipe_idx, val in self.rpool[resource].items():
                return val["json"][key]
        return self.rpool[resource][pipe_idx][key]

    @staticmethod
    def get_resource_from_prov(prov):
        # each resource (i.e. "desc-cleaned_bold" AKA nuisance-regressed BOLD
        # data) has its own provenance list. the name of the resource, and
        # the node that produced it, is always the last item in the provenance
        # list, with the two separated by a colon :
        if not len(prov):
            return None
        if isinstance(prov[-1], list):
            return prov[-1][-1].split(":")[0]
        if isinstance(prov[-1], str):
            return prov[-1].split(":")[0]
        return None

    def regressor_dct(self, cfg) -> dict:
        """Return the regressor dictionary for the current strategy if one exists.

        Raises KeyError otherwise.
        """
        # pylint: disable=attribute-defined-outside-init
        if hasattr(self, "_regressor_dct"):  # memoized
            # pylint: disable=access-member-before-definition
            return self._regressor_dct
        key_error = KeyError(
            "[!] No regressors in resource pool. \n\n"
            "Try turning on create_regressors or "
            "ingress_regressors."
        )
        _nr = cfg["nuisance_corrections", "2-nuisance_regression"]
        if not hasattr(self, "timeseries"):
            if _nr["Regressors"]:
                self.regressors = {reg["Name"]: reg for reg in _nr["Regressors"]}
            else:
                self.regressors = []
        if self.check_rpool("parsed_regressors"):  # ingressed regressor
            # name regressor workflow without regressor_prov
            strat_name = _nr["ingress_regressors"]["Regressors"]["Name"]
            if strat_name in self.regressors:
                self._regressor_dct = self.regressors[strat_name]
                return self._regressor_dct
            self.regressor_dct = _nr["ingress_regressors"]["Regressors"]
            return self.regressor_dct
        prov = self.get_cpac_provenance("desc-confounds_timeseries")
        strat_name_components = prov[-1].split("_")
        for _ in list(range(prov[-1].count("_"))):
            reg_name = "_".join(strat_name_components[-_:])
            if reg_name in self.regressors:
                self._regressor_dct = self.regressors[reg_name]
                return self._regressor_dct
        raise key_error

    def set_data(
        self,
        resource,
        node,
        output,
        json_info,
        pipe_idx,
        node_name,
        fork=False,
        inject=False,
    ):
        json_info = json_info.copy()
        cpac_prov = []
        if "CpacProvenance" in json_info:
            cpac_prov = json_info["CpacProvenance"]
        current_prov_list = list(cpac_prov)
        new_prov_list = list(cpac_prov)  # <---- making a copy, it was already a list
        if not inject:
            new_prov_list.append(f"{resource}:{node_name}")
        try:
            res, new_pipe_idx = self.generate_prov_string(new_prov_list)
        except IndexError:
            msg = (
                f"\n\nThe set_data() call for {resource} has no "
                "provenance information and should not be an "
                "injection."
            )
            raise IndexError(msg)
        if not json_info:
            json_info = {
                "RawSources": [
                    resource  # <---- this will be repopulated to the full file path at the end of the pipeline building, in gather_pipes()
                ]
            }
        json_info["CpacProvenance"] = new_prov_list

        if resource not in self.rpool.keys():
            self.rpool[resource] = {}
        elif not fork:  # <--- in the event of multiple strategies/options, this will run for every option; just keep in mind
            search = False
            if self.get_resource_from_prov(current_prov_list) == resource:
                # CHANGING PIPE_IDX, BE CAREFUL DOWNSTREAM IN THIS FUNCTION
                pipe_idx = self.generate_prov_string(current_prov_list)[1]
                if pipe_idx not in self.rpool[resource].keys():
                    search = True
            else:
                search = True
            if search:
                for idx in current_prov_list:
                    if self.get_resource_from_prov(idx) == resource:
                        if isinstance(idx, list):
                            # CHANGING PIPE_IDX, BE CAREFUL DOWNSTREAM IN THIS FUNCTION
                            pipe_idx = self.generate_prov_string(idx)[1]
                        elif isinstance(idx, str):
                            pipe_idx = idx
                        break
            if pipe_idx in self.rpool[resource].keys():
                # in case the resource name is now new, and not the original
                # remove old keys so we don't end up with a new strat for every new node unit (unless we fork)
                del self.rpool[resource][pipe_idx]
        if new_pipe_idx not in self.rpool[resource]:
            self.rpool[resource][new_pipe_idx] = {}
        if new_pipe_idx not in self.pipe_list:
            self.pipe_list.append(new_pipe_idx)

        self.rpool[resource][new_pipe_idx]["data"] = (node, output)
        self.rpool[resource][new_pipe_idx]["json"] = json_info

    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[str] = None,
        report_fetched: Optional[bool] = False,
        optional: Optional[bool] = False,
    ) -> tuple[Optional[dict], Optional[str]] | Optional[dict]:
        # NOTE!!!
        # if this is the main rpool, this will return a dictionary of strats, and inside those, are dictionaries like {'data': (node, out), 'json': info}
        # BUT, if this is a sub rpool (i.e. a strat_pool), this will return a one-level dictionary of {'data': (node, out), 'json': info} WITHOUT THE LEVEL OF STRAT KEYS ABOVE IT
        if not isinstance(resource, list):
            resource = [resource]
        # if a list of potential inputs are given, pick the first one found
        for label in resource:
            if label in self.rpool.keys():
                _found = self.rpool[label]
                if pipe_idx:
                    _found = _found[pipe_idx]
                if report_fetched:
                    return _found, label
                return _found
        if optional:
            if report_fetched:
                return (None, None)
            return None
        msg = (
            "\n\n[!] C-PAC says: None of the listed resources are in "
            f"the resource pool:\n\n  {resource}\n\nOptions:\n- You "
            "can enable a node block earlier in the pipeline which "
            "produces these resources. Check the 'outputs:' field in "
            "a node block's documentation.\n- You can directly "
            "provide this required data by pulling it from another "
            "BIDS directory using 'source_outputs_dir:' in the "
            "pipeline configuration, or by placing it directly in "
            "your C-PAC output directory.\n- If you have done these, "
            "and you still get this message, please let us know "
            "through any of our support channels at: "
            "https://fcp-indi.github.io/\n"
        )
        raise LookupError(msg)

    def get_data(
        self, resource, pipe_idx=None, report_fetched=False, quick_single=False
    ):
        if report_fetched:
            if pipe_idx:
                connect, fetched = self.get(
                    resource, pipe_idx=pipe_idx, report_fetched=report_fetched
                )
                return (connect["data"], fetched)
            connect, fetched = self.get(resource, report_fetched=report_fetched)
            return (connect["data"], fetched)
        if pipe_idx:
            return self.get(resource, pipe_idx=pipe_idx)["data"]
        if quick_single or len(self.get(resource)) == 1:
            for _key, val in self.get(resource).items():
                return val["data"]
        return self.get(resource)["data"]

    def copy_resource(self, resource, new_name):
        try:
            self.rpool[new_name] = self.rpool[resource]
        except KeyError:
            msg = f"[!] {resource} not in the resource pool."
            raise Exception(msg)

    def update_resource(self, resource, new_name):
        # move over any new pipe_idx's
        self.rpool[new_name].update(self.rpool[resource])

    def get_pipe_idxs(self, resource):
        return self.rpool[resource].keys()

    def get_json(self, resource, strat=None):
        # NOTE: resource_strat_dct has to be entered properly by the developer
        # it has to either be rpool[resource][strat] or strat_pool[resource]
        if strat:
            resource_strat_dct = self.rpool[resource][strat]
        else:
            # for strat_pools mainly, where there is no 'strat' key level
            resource_strat_dct = self.rpool[resource]

        # TODO: the below hits the exception if you use get_cpac_provenance on
        # TODO: the main rpool (i.e. if strat=None)
        if "json" in resource_strat_dct:
            strat_json = resource_strat_dct["json"]
        else:
            msg = (
                "\n[!] Developer info: the JSON "
                f"information for {resource} and {strat} "
                f"is incomplete.\n"
            )
            raise Exception(msg)
        return strat_json

    def get_cpac_provenance(self, resource, strat=None):
        # NOTE: resource_strat_dct has to be entered properly by the developer
        # it has to either be rpool[resource][strat] or strat_pool[resource]
        if isinstance(resource, list):
            for _resource in resource:
                try:
                    return self.get_cpac_provenance(_resource, strat)
                except KeyError:
                    continue
        json_data = self.get_json(resource, strat)
        return json_data["CpacProvenance"]

    @staticmethod
    def generate_prov_string(prov):
        # this will generate a string from a SINGLE RESOURCE'S dictionary of
        # MULTIPLE PRECEDING RESOURCES (or single, if just one)
        # NOTE: this DOES NOT merge multiple resources!!! (i.e. for merging-strat pipe_idx generation)
        if not isinstance(prov, list):
            msg = (
                "\n[!] Developer info: the CpacProvenance "
                f"entry for {prov} has to be a list.\n"
            )
            raise TypeError(msg)
        last_entry = get_last_prov_entry(prov)
        resource = last_entry.split(":")[0]
        return (resource, str(prov))

    @staticmethod
    def generate_prov_list(prov_str):
        if not isinstance(prov_str, str):
            msg = (
                "\n[!] Developer info: the CpacProvenance "
                f"entry for {prov_str!s} has to be a string.\n"
            )
            raise TypeError(msg)
        return ast.literal_eval(prov_str)

    @staticmethod
    def get_resource_strats_from_prov(prov):
        # if you provide the provenance of a resource pool output, this will
        # return a dictionary of all the preceding resource pool entries that
        # led to that one specific output:
        #   {rpool entry}: {that entry's provenance}
        #   {rpool entry}: {that entry's provenance}
        resource_strat_dct = {}
        if isinstance(prov, str):
            resource = prov.split(":")[0]
            resource_strat_dct[resource] = prov
        else:
            for spot, entry in enumerate(prov):
                if isinstance(entry, list):
                    resource = entry[-1].split(":")[0]
                    resource_strat_dct[resource] = entry
                elif isinstance(entry, str):
                    resource = entry.split(":")[0]
                    resource_strat_dct[resource] = entry
        return resource_strat_dct

    def flatten_prov(self, prov):
        if isinstance(prov, str):
            return [prov]
        if isinstance(prov, list):
            flat_prov = []
            for entry in prov:
                if isinstance(entry, list):
                    flat_prov += self.flatten_prov(entry)
                else:
                    flat_prov.append(entry)
            return flat_prov
        return None

    def get_strats(self, resources, debug=False):
        # TODO: NOTE: NOT COMPATIBLE WITH SUB-RPOOL/STRAT_POOLS
        # TODO: (and it doesn't have to be)

        import itertools

        linked_resources = []
        resource_list = []
        if debug:
            verbose_logger = getLogger("CPAC.engine")
            verbose_logger.debug("\nresources: %s", resources)
        for resource in resources:
            # grab the linked-input tuples
            if isinstance(resource, tuple):
                linked = []
                for label in list(resource):
                    rp_dct, fetched_resource = self.get(
                        label, report_fetched=True, optional=True
                    )
                    if not rp_dct:
                        continue
                    linked.append(fetched_resource)
                resource_list += linked
                if len(linked) < 2:  # noqa: PLR2004
                    continue
                linked_resources.append(linked)
            else:
                resource_list.append(resource)

        total_pool = []
        variant_pool = {}
        len_inputs = len(resource_list)
        if debug:
            verbose_logger = getLogger("CPAC.engine")
            verbose_logger.debug("linked_resources: %s", linked_resources)
            verbose_logger.debug("resource_list: %s", resource_list)
        for resource in resource_list:
            (
                rp_dct,  # <---- rp_dct has the strats/pipe_idxs as the keys on first level, then 'data' and 'json' on each strat level underneath
                fetched_resource,
            ) = self.get(
                resource,
                report_fetched=True,
                optional=True,  # oh, and we make the resource fetching in get_strats optional so we can have optional inputs, but they won't be optional in the node block unless we want them to be
            )
            if not rp_dct:
                len_inputs -= 1
                continue
            sub_pool = []
            if debug:
                verbose_logger.debug("len(rp_dct): %s\n", len(rp_dct))
            for strat in rp_dct.keys():
                json_info = self.get_json(fetched_resource, strat)
                cpac_prov = json_info["CpacProvenance"]
                sub_pool.append(cpac_prov)
                if fetched_resource not in variant_pool:
                    variant_pool[fetched_resource] = []
                if "CpacVariant" in json_info:
                    for key, val in json_info["CpacVariant"].items():
                        if val not in variant_pool[fetched_resource]:
                            variant_pool[fetched_resource] += val
                            variant_pool[fetched_resource].append(f"NO-{val[0]}")

            if debug:
                verbose_logger = getLogger("CPAC.engine")
                verbose_logger.debug("%s sub_pool: %s\n", resource, sub_pool)
            total_pool.append(sub_pool)

        if not total_pool:
            raise LookupError(
                "\n\n[!] C-PAC says: None of the listed "
                "resources in the node block being connected "
                "exist in the resource pool.\n\nResources:\n"
                "%s\n\n" % resource_list
            )

        # TODO: right now total_pool is:
        # TODO:    [[[T1w:anat_ingress, desc-preproc_T1w:anatomical_init, desc-preproc_T1w:acpc_alignment], [T1w:anat_ingress,desc-preproc_T1w:anatomical_init]],
        # TODO:     [[T1w:anat_ingress, desc-preproc_T1w:anatomical_init, desc-preproc_T1w:acpc_alignment, desc-brain_mask:brain_mask_afni], [T1w:anat_ingress, desc-preproc_T1w:anatomical_init, desc-brain_mask:brain_mask_afni]]]

        # TODO: and the code below thinks total_pool is a list of lists, like [[pipe_idx, pipe_idx], [pipe_idx, pipe_idx, pipe_idx], etc.]
        # TODO: and the actual resource is encoded in the tag: of the last item, every time!
        # keying the strategies to the resources, inverting it
        if len_inputs > 1:
            strats = itertools.product(*total_pool)

            # we now currently have "strats", the combined permutations of all the strategies, as a list of tuples, each tuple combining one version of input each, being one of the permutations.
            # OF ALL THE DIFFERENT INPUTS. and they are tagged by their fetched inputs with {name}:{strat}.
            # so, each tuple has ONE STRAT FOR EACH INPUT, so if there are three inputs, each tuple will have 3 items.
            new_strats = {}

            # get rid of duplicates - TODO: refactor .product
            strat_str_list = []
            strat_list_list = []
            for strat_tuple in strats:
                strat_list = list(copy.deepcopy(strat_tuple))
                strat_str = str(strat_list)
                if strat_str not in strat_str_list:
                    strat_str_list.append(strat_str)
                    strat_list_list.append(strat_list)

            if debug:
                verbose_logger = getLogger("CPAC.engine")
                verbose_logger.debug("len(strat_list_list): %s\n", len(strat_list_list))
            for strat_list in strat_list_list:
                json_dct = {}
                for strat in strat_list:
                    # strat is a prov list for a single resource/input
                    strat_resource, strat_idx = self.generate_prov_string(strat)
                    strat_json = self.get_json(strat_resource, strat=strat_idx)
                    json_dct[strat_resource] = strat_json

                drop = False
                if linked_resources:
                    for linked in linked_resources:  # <--- 'linked' is each tuple
                        if drop:
                            break
                        for xlabel in linked:
                            if drop:
                                break
                            xjson = copy.deepcopy(json_dct[xlabel])
                            for ylabel in linked:
                                if xlabel == ylabel:
                                    continue
                                yjson = copy.deepcopy(json_dct[ylabel])

                                if "CpacVariant" not in xjson:
                                    xjson["CpacVariant"] = {}
                                if "CpacVariant" not in yjson:
                                    yjson["CpacVariant"] = {}

                                current_strat = []
                                for key, val in xjson["CpacVariant"].items():
                                    if isinstance(val, list):
                                        current_strat.append(val[0])
                                    else:
                                        current_strat.append(val)
                                current_spread = list(set(variant_pool[xlabel]))
                                for spread_label in current_spread:
                                    if "NO-" in spread_label:
                                        continue
                                    if spread_label not in current_strat:
                                        current_strat.append(f"NO-{spread_label}")

                                other_strat = []
                                for key, val in yjson["CpacVariant"].items():
                                    if isinstance(val, list):
                                        other_strat.append(val[0])
                                    else:
                                        other_strat.append(val)
                                other_spread = list(set(variant_pool[ylabel]))
                                for spread_label in other_spread:
                                    if "NO-" in spread_label:
                                        continue
                                    if spread_label not in other_strat:
                                        other_strat.append(f"NO-{spread_label}")

                                for variant in current_spread:
                                    in_current_strat = False
                                    in_other_strat = False
                                    in_other_spread = False

                                    if variant is None:
                                        in_current_strat = True
                                        if None in other_spread:
                                            in_other_strat = True
                                    if variant in current_strat:
                                        in_current_strat = True
                                    if variant in other_strat:
                                        in_other_strat = True
                                    if variant in other_spread:
                                        in_other_spread = True

                                    if not in_other_strat:
                                        if in_other_spread:
                                            if in_current_strat:
                                                drop = True
                                                break

                                    if in_other_strat:
                                        if in_other_spread:
                                            if not in_current_strat:
                                                drop = True
                                                break
                                if drop:
                                    break
                if drop:
                    continue

                # make the merged strat label from the multiple inputs
                # strat_list is actually the merged CpacProvenance lists
                pipe_idx = str(strat_list)
                new_strats[pipe_idx] = ResourcePool()
                # new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                # placing JSON info at one level higher only for copy convenience
                new_strats[pipe_idx].rpool["json"] = {}
                new_strats[pipe_idx].rpool["json"]["subjson"] = {}
                new_strats[pipe_idx].rpool["json"]["CpacProvenance"] = strat_list

                # now just invert resource:strat to strat:resource for each resource:strat
                for cpac_prov in strat_list:
                    resource, strat = self.generate_prov_string(cpac_prov)
                    resource_strat_dct = self.rpool[resource][strat]
                    # remember, `resource_strat_dct` is the dct of 'data' and 'json'.
                    new_strats[pipe_idx].rpool[resource] = resource_strat_dct
                    # `new_strats` is A DICTIONARY OF RESOURCEPOOL OBJECTS! each one is a new slice of the resource pool combined together.
                    self.pipe_list.append(pipe_idx)
                    if "CpacVariant" in resource_strat_dct["json"]:
                        if "CpacVariant" not in new_strats[pipe_idx].rpool["json"]:
                            new_strats[pipe_idx].rpool["json"]["CpacVariant"] = {}
                        for younger_resource, variant_list in resource_strat_dct[
                            "json"
                        ]["CpacVariant"].items():
                            if (
                                younger_resource
                                not in new_strats[pipe_idx].rpool["json"]["CpacVariant"]
                            ):
                                new_strats[pipe_idx].rpool["json"]["CpacVariant"][
                                    younger_resource
                                ] = variant_list
                    # preserve each input's JSON info also
                    data_type = resource.split("_")[-1]
                    if data_type not in new_strats[pipe_idx].rpool["json"]["subjson"]:
                        new_strats[pipe_idx].rpool["json"]["subjson"][data_type] = {}
                    new_strats[pipe_idx].rpool["json"]["subjson"][data_type].update(
                        copy.deepcopy(resource_strat_dct["json"])
                    )
        else:
            new_strats = {}
            for resource_strat_list in total_pool:
                # total_pool will have only one list of strats, for the one input
                for cpac_prov in resource_strat_list:  # <------- cpac_prov here doesn't need to be modified, because it's not merging with other inputs
                    resource, pipe_idx = self.generate_prov_string(cpac_prov)
                    resource_strat_dct = self.rpool[resource][pipe_idx]
                    # remember, `resource_strat_dct` is the dct of 'data' and 'json'.
                    new_strats[pipe_idx] = ResourcePool(
                        rpool={resource: resource_strat_dct}
                    )  # <----- again, new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                    # placing JSON info at one level higher only for copy convenience
                    new_strats[pipe_idx].rpool["json"] = resource_strat_dct["json"]
                    # TODO: WARNING- THIS IS A LEVEL HIGHER THAN THE ORIGINAL 'JSON' FOR EASE OF ACCESS IN CONNECT_BLOCK WITH THE .GET(JSON)
                    new_strats[pipe_idx].rpool["json"]["subjson"] = {}
                    new_strats[pipe_idx].rpool["json"]["CpacProvenance"] = cpac_prov
                    # preserve each input's JSON info also
                    data_type = resource.split("_")[-1]
                    if data_type not in new_strats[pipe_idx].rpool["json"]["subjson"]:
                        new_strats[pipe_idx].rpool["json"]["subjson"][data_type] = {}
                    new_strats[pipe_idx].rpool["json"]["subjson"][data_type].update(
                        copy.deepcopy(resource_strat_dct["json"])
                    )
        return new_strats

    def derivative_xfm(self, wf, label, connection, json_info, pipe_idx, pipe_x):
        if label in self.xfm:
            json_info = dict(json_info)

            # get the bold-to-template transform from the current strat_pool info
            xfm_idx = None
            xfm_label = "from-bold_to-template_mode-image_xfm"
            for entry in json_info["CpacProvenance"]:
                if isinstance(entry, list):
                    if entry[-1].split(":")[0] == xfm_label:
                        xfm_prov = entry
                        xfm_idx = self.generate_prov_string(xfm_prov)[1]
                        break

            # but if the resource doesn't have the bold-to-template transform
            # in its provenance/strategy, find the appropriate one for this
            # current pipe_idx/strat
            if not xfm_idx:
                xfm_info = []
                for pipe_idx, entry in self.get(xfm_label).items():
                    xfm_info.append((pipe_idx, entry["json"]["CpacProvenance"]))
            else:
                xfm_info = [(xfm_idx, xfm_prov)]

            for num, xfm_entry in enumerate(xfm_info):
                xfm_idx, xfm_prov = xfm_entry
                reg_tool = check_prov_for_regtool(xfm_prov)

                xfm = transform_derivative(
                    f"{label}_xfm_{pipe_x}_{num}",
                    label,
                    reg_tool,
                    self.num_cpus,
                    self.num_ants_cores,
                    ants_interp=self.ants_interp,
                    fsl_interp=self.fsl_interp,
                    opt=None,
                )
                wf.connect(connection[0], connection[1], xfm, "inputspec.in_file")

                node, out = self.get_data("T1w-brain-template-deriv", quick_single=True)
                wf.connect(node, out, xfm, "inputspec.reference")

                node, out = self.get_data(
                    "from-bold_to-template_mode-image_xfm", pipe_idx=xfm_idx
                )
                wf.connect(node, out, xfm, "inputspec.transform")

                label = f"space-template_{label}"
                json_info["Template"] = self.get_json_info(
                    "T1w-brain-template-deriv", None, "Description"
                )
                new_prov = json_info["CpacProvenance"] + xfm_prov
                json_info["CpacProvenance"] = new_prov
                new_pipe_idx = self.generate_prov_string(new_prov)
                self.set_data(
                    label,
                    xfm,
                    "outputspec.out_file",
                    json_info,
                    new_pipe_idx,
                    f"{label}_xfm_{num}",
                    fork=True,
                )

        return wf

    
    @property
    def filtered_movement(self) -> bool:
        """
        Check if the movement parameters have been filtered in this strat_pool.

        Returns
        -------
        bool
        """
        try:
            return "motion_estimate_filter" in str(
                self.get_cpac_provenance("desc-movementParameters_motion")
            )
        except KeyError:
            # not a strat_pool or no movement parameters in strat_pool
            return False

    def filter_name(self, cfg: Configuration) -> str:
        """
        Return the name of the filter for this strategy.

        In a strat_pool with filtered movement parameters.
        """
        motion_filters = cfg[
            "functional_preproc",
            "motion_estimates_and_correction",
            "motion_estimate_filter",
            "filters",
        ]
        if len(motion_filters) == 1 and cfg.switch_is_on(
            [
                "functional_preproc",
                "motion_estimates_and_correction",
                "motion_estimate_filter",
                "run",
            ],
            exclusive=True,
        ):
            return motion_filters[0]["Name"]
        try:
            key = "motion"
            sidecar = self.get_json("desc-movementParameters_motion")
        except KeyError:
            sidecar = None
        if sidecar is not None and "CpacVariant" in sidecar:
            if sidecar["CpacVariant"][key]:
                return sidecar["CpacVariant"][key][0][::-1].split("_", 1)[0][::-1]
        return "none"

    def post_process(self, wf, label, connection, json_info, pipe_idx, pipe_x, outs):
        input_type = "func_derivative"

        post_labels = [(label, connection[0], connection[1])]

        if re.match(r"(.*_)?[ed]c[bw]$", label) or re.match(r"(.*_)?lfcd[bw]$", label):
            # suffix: [eigenvector or degree] centrality [binarized or weighted]
            # or lfcd [binarized or weighted]
            mask = "template-specification-file"
        elif "space-template" in label:
            if "space-template_res-derivative_desc-bold_mask" in self.rpool.keys():
                mask = "space-template_res-derivative_desc-bold_mask"
            else:
                mask = "space-template_desc-bold_mask"
        else:
            mask = "space-bold_desc-brain_mask"

        mask_idx = None
        for entry in json_info["CpacProvenance"]:
            if isinstance(entry, list):
                if entry[-1].split(":")[0] == mask:
                    mask_prov = entry
                    mask_idx = self.generate_prov_string(mask_prov)[1]
                    break

        if self.smoothing_bool:
            if label in Outputs.to_smooth:
                for smooth_opt in self.smooth_opts:
                    sm = spatial_smoothing(
                        f"{label}_smooth_{smooth_opt}_{pipe_x}",
                        self.fwhm,
                        input_type,
                        smooth_opt,
                    )
                    wf.connect(connection[0], connection[1], sm, "inputspec.in_file")
                    node, out = self.get_data(
                        mask, pipe_idx=mask_idx, quick_single=mask_idx is None
                    )
                    wf.connect(node, out, sm, "inputspec.mask")

                    if "desc-" not in label:
                        if "space-" in label:
                            for tag in label.split("_"):
                                if "space-" in tag:
                                    smlabel = label.replace(tag, f"{tag}_desc-sm")
                                    break
                        else:
                            smlabel = f"desc-sm_{label}"
                    else:
                        for tag in label.split("_"):
                            if "desc-" in tag:
                                newtag = f"{tag}-sm"
                                smlabel = label.replace(tag, newtag)
                                break

                    post_labels.append((smlabel, sm, "outputspec.out_file"))

                    self.set_data(
                        smlabel,
                        sm,
                        "outputspec.out_file",
                        json_info,
                        pipe_idx,
                        f"spatial_smoothing_{smooth_opt}",
                        fork=True,
                    )
                    self.set_data(
                        "fwhm",
                        sm,
                        "outputspec.fwhm",
                        json_info,
                        pipe_idx,
                        f"spatial_smoothing_{smooth_opt}",
                        fork=True,
                    )

        if self.zscoring_bool:
            for label_con_tpl in post_labels:
                label = label_con_tpl[0]
                connection = (label_con_tpl[1], label_con_tpl[2])
                if label in Outputs.to_zstd:
                    zstd = z_score_standardize(f"{label}_zstd_{pipe_x}", input_type)

                    wf.connect(connection[0], connection[1], zstd, "inputspec.in_file")

                    node, out = self.get_data(mask, pipe_idx=mask_idx)
                    wf.connect(node, out, zstd, "inputspec.mask")

                    if "desc-" not in label:
                        if "space-template" in label:
                            new_label = label.replace(
                                "space-template", "space-template_desc-zstd"
                            )
                        else:
                            new_label = f"desc-zstd_{label}"
                    else:
                        for tag in label.split("_"):
                            if "desc-" in tag:
                                newtag = f"{tag}-zstd"
                                new_label = label.replace(tag, newtag)
                                break

                    post_labels.append((new_label, zstd, "outputspec.out_file"))

                    self.set_data(
                        new_label,
                        zstd,
                        "outputspec.out_file",
                        json_info,
                        pipe_idx,
                        "zscore_standardize",
                        fork=True,
                    )

                elif label in Outputs.to_fisherz:
                    zstd = fisher_z_score_standardize(
                        f"{label}_zstd_{pipe_x}", label, input_type
                    )

                    wf.connect(
                        connection[0], connection[1], zstd, "inputspec.correlation_file"
                    )

                    # if the output is 'space-template_desc-MeanSCA_correlations', we want 'desc-MeanSCA_timeseries'
                    oned = label.replace("correlations", "timeseries")

                    node, out = outs[oned]
                    wf.connect(node, out, zstd, "inputspec.timeseries_oned")

                    post_labels.append((new_label, zstd, "outputspec.out_file"))

                    self.set_data(
                        new_label,
                        zstd,
                        "outputspec.out_file",
                        json_info,
                        pipe_idx,
                        "fisher_zscore_standardize",
                        fork=True,
                    )

        return (wf, post_labels)

    def gather_pipes(self, wf, cfg, all=False, add_incl=None, add_excl=None):
        excl = []
        substring_excl = []
        outputs_logger = getLogger(f'{cfg["subject_id"]}_expectedOutputs')
        expected_outputs = ExpectedOutputs()

        if add_excl:
            excl += add_excl

        if "nonsmoothed" not in cfg.post_processing["spatial_smoothing"]["output"]:
            excl += Outputs.native_nonsmooth
            excl += Outputs.template_nonsmooth

        if "raw" not in cfg.post_processing["z-scoring"]["output"]:
            excl += Outputs.native_raw
            excl += Outputs.template_raw

        if not cfg.pipeline_setup["output_directory"]["write_debugging_outputs"]:
            # substring_excl.append(['bold'])
            excl += Outputs.debugging

        for resource in self.rpool.keys():
            if resource not in Outputs.any:
                continue

            if resource in excl:
                continue

            drop = False
            for substring_list in substring_excl:
                bool_list = []
                for substring in substring_list:
                    if substring in resource:
                        bool_list.append(True)
                    else:
                        bool_list.append(False)
                for item in bool_list:
                    if not item:
                        break
                else:
                    drop = True
                if drop:
                    break
            if drop:
                continue

            subdir = "other"
            if resource in Outputs.anat:
                subdir = "anat"
                # TODO: get acq- etc.
            elif resource in Outputs.func:
                subdir = "func"
                # TODO: other stuff like acq- etc.

            for pipe_idx in self.rpool[resource]:
                unique_id = self.get_name()
                part_id = unique_id.split("_")[0]
                ses_id = unique_id.split("_")[1]

                if "ses-" not in ses_id:
                    ses_id = f"ses-{ses_id}"

                out_dir = cfg.pipeline_setup["output_directory"]["path"]
                pipe_name = cfg.pipeline_setup["pipeline_name"]
                container = os.path.join(f"pipeline_{pipe_name}", part_id, ses_id)
                filename = f"{unique_id}_{res_in_filename(self.cfg, resource)}"

                out_path = os.path.join(out_dir, container, subdir, filename)

                out_dct = {
                    "unique_id": unique_id,
                    "out_dir": out_dir,
                    "container": container,
                    "subdir": subdir,
                    "filename": filename,
                    "out_path": out_path,
                }
                self.rpool[resource][pipe_idx]["out"] = out_dct

                # TODO: have to link the pipe_idx's here. and call up 'desc-preproc_T1w' from a Sources in a json and replace. here.
                # TODO: can do the pipeline_description.json variants here too!

        for resource in self.rpool.keys():
            if resource not in Outputs.any:
                continue

            if resource in excl:
                continue

            drop = False
            for substring_list in substring_excl:
                bool_list = []
                for substring in substring_list:
                    if substring in resource:
                        bool_list.append(True)
                    else:
                        bool_list.append(False)
                for item in bool_list:
                    if not item:
                        break
                else:
                    drop = True
                if drop:
                    break
            if drop:
                continue

            num_variant = 0
            if len(self.rpool[resource]) == 1:
                num_variant = ""
            all_jsons = [
                self.rpool[resource][pipe_idx]["json"]
                for pipe_idx in self.rpool[resource]
            ]
            unlabelled = {
                key
                for json_info in all_jsons
                for key in json_info.get("CpacVariant", {}).keys()
                if key not in (*MOVEMENT_FILTER_KEYS, "regressors")
            }
            if "bold" in unlabelled:
                all_bolds = list(
                    chain.from_iterable(
                        json_info["CpacVariant"]["bold"]
                        for json_info in all_jsons
                        if "CpacVariant" in json_info
                        and "bold" in json_info["CpacVariant"]
                    )
                )
                # not any(not) because all is overloaded as a parameter here
                if not any(
                    not re.match(
                        r"apply_(phasediff|blip)_to_timeseries_separately_.*", _bold
                    )
                    for _bold in all_bolds
                ):
                    # this fork point should only result in 0 or 1 forks
                    unlabelled.remove("bold")
                del all_bolds
            all_forks = {
                key: set(
                    chain.from_iterable(
                        json_info["CpacVariant"][key]
                        for json_info in all_jsons
                        if "CpacVariant" in json_info
                        and key in json_info["CpacVariant"]
                    )
                )
                for key in unlabelled
            }
            # del all_jsons
            for key, forks in all_forks.items():
                if len(forks) < 2:  # noqa: PLR2004
                    # no int suffix needed if only one fork
                    unlabelled.remove(key)
            # del all_forks
            for pipe_idx in self.rpool[resource]:
                pipe_x = self.get_pipe_number(pipe_idx)
                json_info = self.rpool[resource][pipe_idx]["json"]
                out_dct = self.rpool[resource][pipe_idx]["out"]

                try:
                    if unlabelled:
                        num_variant += 1
                except TypeError:
                    pass

                try:
                    del json_info["subjson"]
                except KeyError:
                    pass

                if out_dct["subdir"] == "other" and not all:
                    continue

                unique_id = out_dct["unique_id"]
                resource_idx = resource

                if isinstance(num_variant, int):
                    resource_idx, out_dct = name_fork(
                        resource_idx, cfg, json_info, out_dct
                    )
                    if unlabelled:
                        if "desc-" in out_dct["filename"]:
                            for key in out_dct["filename"].split("_")[::-1]:
                                # final `desc` entity
                                if key.startswith("desc-"):
                                    out_dct["filename"] = out_dct["filename"].replace(
                                        key, f"{key}-{num_variant}"
                                    )
                                    resource_idx = resource_idx.replace(
                                        key, f"{key}-{num_variant}"
                                    )
                                    break
                        else:
                            suff = resource.split("_")[-1]
                            newdesc_suff = f"desc-{num_variant}_{suff}"
                            resource_idx = resource_idx.replace(suff, newdesc_suff)
                id_string = pe.Node(
                    Function(
                        input_names=[
                            "cfg",
                            "unique_id",
                            "resource",
                            "scan_id",
                            "template_desc",
                            "atlas_id",
                            "fwhm",
                            "subdir",
                            "extension",
                        ],
                        output_names=["out_filename"],
                        function=create_id_string,
                    ),
                    name=f"id_string_{resource_idx}_{pipe_x}",
                )
                id_string.inputs.cfg = self.cfg
                id_string.inputs.unique_id = unique_id
                id_string.inputs.resource = resource_idx
                id_string.inputs.subdir = out_dct["subdir"]

                # grab the iterable scan ID
                if out_dct["subdir"] == "func":
                    node, out = self.rpool["scan"]["['scan:func_ingress']"]["data"]
                    wf.connect(node, out, id_string, "scan_id")

                self.back_propogate_template_name(
                    wf, resource_idx, json_info, id_string
                )
                # grab the FWHM if smoothed
                for tag in resource.split("_"):
                    if "desc-" in tag and "-sm" in tag:
                        fwhm_idx = pipe_idx.replace(f"{resource}:", "fwhm:")
                        try:
                            node, out = self.rpool["fwhm"][fwhm_idx]["data"]
                            wf.connect(node, out, id_string, "fwhm")
                        except KeyError:
                            # smoothing was not done for this resource in the
                            # engine.py smoothing
                            pass
                        break
                atlas_suffixes = ["timeseries", "correlations", "statmap"]
                # grab the iterable atlas ID
                atlas_id = None
                if not resource.endswith("desc-confounds_timeseries"):
                    if resource.split("_")[-1] in atlas_suffixes:
                        atlas_idx = pipe_idx.replace(resource, "atlas_name")
                        # need the single quote and the colon inside the double
                        # quotes - it's the encoded pipe_idx
                        # atlas_idx = new_idx.replace(f"'{temp_rsc}:",
                        #                             "'atlas_name:")
                        if atlas_idx in self.rpool["atlas_name"]:
                            node, out = self.rpool["atlas_name"][atlas_idx]["data"]
                            wf.connect(node, out, id_string, "atlas_id")
                        elif "atlas-" in resource:
                            for tag in resource.split("_"):
                                if "atlas-" in tag:
                                    atlas_id = tag.replace("atlas-", "")
                            id_string.inputs.atlas_id = atlas_id
                        else:
                            warnings.warn(
                                str(
                                    LookupError(
                                        "\n[!] No atlas ID found for "
                                        f"{out_dct['filename']}.\n"
                                    )
                                )
                            )
                nii_name = pe.Node(Rename(), name=f"nii_{resource_idx}_{pipe_x}")
                nii_name.inputs.keep_ext = True

                if resource in Outputs.ciftis:
                    nii_name.inputs.keep_ext = False
                    id_string.inputs.extension = Outputs.ciftis[resource]
                else:
                    nii_name.inputs.keep_ext = True

                if resource in Outputs.giftis:
                    nii_name.inputs.keep_ext = False
                    id_string.inputs.extension = f"{Outputs.giftis[resource]}.gii"

                else:
                    nii_name.inputs.keep_ext = True

                wf.connect(id_string, "out_filename", nii_name, "format_string")

                node, out = self.rpool[resource][pipe_idx]["data"]
                try:
                    wf.connect(node, out, nii_name, "in_file")
                except OSError as os_error:
                    WFLOGGER.warning(os_error)
                    continue

                write_json_imports = ["import os", "import json"]
                write_json = pe.Node(
                    Function(
                        input_names=["json_data", "filename"],
                        output_names=["json_file"],
                        function=write_output_json,
                        imports=write_json_imports,
                    ),
                    name=f"json_{resource_idx}_{pipe_x}",
                )
                write_json.inputs.json_data = json_info

                wf.connect(id_string, "out_filename", write_json, "filename")
                ds = pe.Node(DataSink(), name=f"sinker_{resource_idx}_{pipe_x}")
                ds.inputs.parameterization = False
                ds.inputs.base_directory = out_dct["out_dir"]
                ds.inputs.encrypt_bucket_keys = cfg.pipeline_setup["Amazon-AWS"][
                    "s3_encryption"
                ]
                ds.inputs.container = out_dct["container"]

                if cfg.pipeline_setup["Amazon-AWS"]["aws_output_bucket_credentials"]:
                    ds.inputs.creds_path = cfg.pipeline_setup["Amazon-AWS"][
                        "aws_output_bucket_credentials"
                    ]
                expected_outputs += (
                    out_dct["subdir"],
                    create_id_string(
                        self.cfg,
                        unique_id,
                        resource_idx,
                        template_desc=id_string.inputs.template_desc,
                        atlas_id=atlas_id,
                        subdir=out_dct["subdir"],
                    ),
                )
                wf.connect(nii_name, "out_file", ds, f'{out_dct["subdir"]}.@data')
                wf.connect(write_json, "json_file", ds, f'{out_dct["subdir"]}.@json')
        outputs_logger.info(expected_outputs)

    def node_data(self, resource, **kwargs):
        """Create NodeData objects.

        Parameters
        ----------
        resource : str

        Returns
        -------
        NodeData
        """
        return NodeData(self, resource, **kwargs)

    def build_rpool(
            self, 
            data_paths, 
            default_CpacProvenance = ('ingress')):
        count = 1
        for index, row in data_paths[1].iterrows():
            # Check if 'meta__json' is not None and contains 'CpacProvenance'
            if row.get('meta__json') and row['meta__json'].get('CpacProvenance'):
                CpacProvenance = row['meta__json']['CpacProvenance']
            else:
                CpacProvenance = default_CpacProvenance
            resource = Resource(row, CpacProvenance)
            # making the rpool a list so that the duplicates are appended rather than overwritten
            self.rpool.setdefault(resource.suffix, [])
            self.rpool[resource.suffix].append(resource)
            # count += 1
            # if count >10:
            #     break
            
    
    def write_to_disk(self, path):
        for resources in self.rpool.values():
            for item in resources:
                print(item['resource'].write_to_disk(path))

    def get_resource(self, description):
        matching_resources = []
        for resources in self.rpool.get(description['suffix'], []):
            # Initialize a flag to True, assuming the resource matches until proven otherwise
            is_match = True
            for key, val in description.items():
                # Skip the 'suffix' key as it's used to select the pool, not to match resources
                if key == 'suffix':
                    continue
                # Check if the resource matches the description criteria
                # Use getattr for object attributes or resources.get for dictionary keys
                resource_val = getattr(resources, key, None)
                if resource_val.lower() != val.lower():
                    is_match = False
                    break  # Break out of the inner loop if any criteria does not match
            if is_match:
                # If the resource matches all criteria, append its name to the matching_resources list
                matching_resources.append(resources.name)
        for items in matching_resources: 
            print(items)
        return matching_resources
            

    def set_resource(self, name, value):
        self.rpool[name] = value

