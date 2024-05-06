# Copyright (C) 2021-2023  C-PAC Developers

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
import json
from itertools import chain
import logging
import os
import re
from typing import Any, Optional, Union
import warnings


from CPAC.pipeline import \
    nipype_pipeline_engine as pe  # pylint: disable=ungrouped-imports
from nipype import config, logging  # pylint: disable=wrong-import-order
from CPAC.pipeline.nodeblock import NodeBlockFunction  # pylint: disable=ungrouped-imports
from nipype.interfaces.utility import \
    Rename  # pylint: disable=wrong-import-order
from CPAC.image_utils.spatial_smoothing import spatial_smoothing
from CPAC.image_utils.statistical_transforms import z_score_standardize, \
    fisher_z_score_standardize
from CPAC.pipeline.check_outputs import ExpectedOutputs
from CPAC.pipeline.utils import MOVEMENT_FILTER_KEYS, name_fork, source_set
from CPAC.registration.registration import transform_derivative
from CPAC.utils.bids_utils import res_in_filename
from CPAC.utils.datasource import (
    create_anat_datasource,
    create_func_datasource,
    ingress_func_metadata,
    create_general_datasource,
    resolve_resolution
)
from CPAC.utils.interfaces.function import Function
from CPAC.utils.interfaces.datasink import DataSink
from CPAC.utils.monitoring import getLogger, LOGTAIL, \
                                  WARNING_FREESURFER_OFF_WITH_DATA
from CPAC.utils.outputs import Outputs
from CPAC.utils.typing import LIST_OR_STR, TUPLE
from CPAC.utils.utils import check_prov_for_regtool, \
    create_id_string, get_last_prov_entry, read_json, write_output_json

from CPAC.resources.templates.lookup_table import lookup_identifier

logger = getLogger('nipype.workflow')


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
            self.logdir = cfg.pipeline_setup['log_directory']['path']

            self.num_cpus = cfg.pipeline_setup['system_config'][
                'max_cores_per_participant']
            self.num_ants_cores = cfg.pipeline_setup['system_config'][
                'num_ants_threads']

            self.ants_interp = cfg.registration_workflows[
                'functional_registration']['func_registration_to_template'][
                'ANTs_pipelines']['interpolation']
            self.fsl_interp = cfg.registration_workflows[
                'functional_registration']['func_registration_to_template'][
                'FNIRT_pipelines']['interpolation']

            self.func_reg = cfg.registration_workflows[
                'functional_registration']['func_registration_to_template'][
                'run']

            self.run_smoothing = 'smoothed' in cfg.post_processing[
                'spatial_smoothing']['output']
            self.smoothing_bool = cfg.post_processing['spatial_smoothing']['run']
            self.run_zscoring = 'z-scored' in cfg.post_processing[
                'z-scoring']['output']
            self.zscoring_bool = cfg.post_processing['z-scoring']['run']
            self.fwhm = cfg.post_processing['spatial_smoothing']['fwhm']
            self.smooth_opts = cfg.post_processing['spatial_smoothing'][
                'smoothing_method']

        self.xfm = ['alff', 'desc-sm_alff', 'desc-zstd_alff', 
                    'desc-sm-zstd_alff',
                    'falff', 'desc-sm_falff', 'desc-zstd_falff',
                    'desc-sm-zstd_falff',
                    'reho', 'desc-sm_reho', 'desc-zstd_reho',
                    'desc-sm-zstd_reho']

    def __repr__(self) -> str:
        params = [f"{param}={getattr(self, param)}" for param in
                  ["rpool", "name", "cfg", "pipe_list"] if
                  getattr(self, param, None) is not None]
        return f'ResourcePool({", ".join(params)})'

    def __str__(self) -> str:
        if self.name:
            return f'ResourcePool({self.name}): {list(self.rpool)}'
        return f'ResourcePool: {list(self.rpool)}'

    def append_name(self, name):
        self.name.append(name)

    def back_propogate_template_name(self, wf, resource_idx: str, json_info: dict,
                                     id_string: 'pe.Node') -> None:
        """Find and apply the template name from a resource's provenance

        Parameters
        ----------
        resource_idx : str

        json_info : dict

        id_string : pe.Node

        Returns
        -------
        None
        """
        if ('template' in resource_idx and self.check_rpool('derivatives-dir')):
            if self.check_rpool('template'):
                node, out = self.get_data('template')
                wf.connect(node, out, id_string, 'template_desc')
        elif 'Template' in json_info:
            id_string.inputs.template_desc = json_info['Template']
        elif ('template' in resource_idx and
              len(json_info.get('CpacProvenance', [])) > 1):
            for resource in source_set(json_info['CpacProvenance']):
                source, value = resource.split(':', 1)
                if value.startswith('template_'
                                    ) and source != 'FSL-AFNI-bold-ref':
                    # 'FSL-AFNI-bold-ref' is currently allowed to be in
                    # a different space, so don't use it as the space for
                    # descendents
                    try:
                        anscestor_json = list(self.rpool.get(source).items()
                                              )[0][1].get('json', {})
                        if 'Description' in anscestor_json:
                            id_string.inputs.template_desc = anscestor_json[
                                'Description']
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
        return ResourcePool(rpool=copy.deepcopy(self.get_entire_rpool()),
                            name=self.name,
                            cfg=self.cfg,
                            pipe_list=copy.deepcopy(self.pipe_list))

    @staticmethod
    def get_raw_label(resource: str) -> str:
        """Removes ``desc-*`` label"""
        for tag in resource.split('_'):
            if 'desc-' in tag:
                resource = resource.replace(f'{tag}_', '')
                break
        return resource

    def get_strat_info(self, prov, label=None, logdir=None):
        strat_info = {}
        for entry in prov:
            if isinstance(entry, list):
                strat_info[entry[-1].split(':')[0]] = entry
            elif isinstance(entry, str):
                strat_info[entry.split(':')[0]] = entry.split(':')[1]
        if label:
            if not logdir:
                logdir = self.logdir
            print(f'\n\nPrinting out strategy info for {label} in {logdir}\n')
            write_output_json(strat_info, f'{label}_strat_info',
                              indent=4, basedir=logdir)

    def set_json_info(self, resource, pipe_idx, key, val):
        #TODO: actually should probably be able to inititialize resource/pipe_idx
        if pipe_idx not in self.rpool[resource]:
            raise Exception('\n[!] DEV: The pipeline/strat ID does not exist '
                            f'in the resource pool.\nResource: {resource}'
                            f'Pipe idx: {pipe_idx}\nKey: {key}\nVal: {val}\n')
        else:
            if 'json' not in self.rpool[resource][pipe_idx]:
                self.rpool[resource][pipe_idx]['json'] = {}
            self.rpool[resource][pipe_idx]['json'][key] = val

    def get_json_info(self, resource, pipe_idx, key):
        #TODO: key checks
        if not pipe_idx:
           for pipe_idx, val in self.rpool[resource].items():
                return val['json'][key]
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
            return prov[-1][-1].split(':')[0]
        elif isinstance(prov[-1], str):
            return prov[-1].split(':')[0]

    def regressor_dct(self, cfg) -> dict:
        """Returns the regressor dictionary for the current strategy if
        one exists. Raises KeyError otherwise."""
        # pylint: disable=attribute-defined-outside-init
        if hasattr(self, '_regressor_dct'):  # memoized
            # pylint: disable=access-member-before-definition
            return self._regressor_dct
        key_error = KeyError("[!] No regressors in resource pool. \n\n"
                             "Try turning on create_regressors or "
                             "ingress_regressors.")
        _nr = cfg['nuisance_corrections', '2-nuisance_regression']
        
        if not hasattr(self, 'timeseries'):
            if _nr['Regressors']:
                self.regressors = {reg["Name"]: reg for reg in _nr['Regressors']}
            else:
                self.regressors = []
        if self.check_rpool('parsed_regressors'):  # ingressed regressor
            # name regressor workflow without regressor_prov
            strat_name = _nr['ingress_regressors']['Regressors']['Name']
            if strat_name in self.regressors:
                self._regressor_dct = self.regressors[strat_name]
                return self._regressor_dct
            self.regressor_dct = _nr['ingress_regressors']['Regressors']
            return self.regressor_dct
        prov = self.get_cpac_provenance('desc-confounds_timeseries')
        strat_name_components = prov[-1].split('_')
        for _ in list(range(prov[-1].count('_'))):
            reg_name = '_'.join(strat_name_components[-_:])
            if reg_name in self.regressors:
                self._regressor_dct = self.regressors[reg_name]
                return self._regressor_dct
        raise key_error

    def set_data(self, resource, node, output, json_info, pipe_idx, node_name,
                 fork=False, inject=False):
        json_info = json_info.copy()
        cpac_prov = []
        if 'CpacProvenance' in json_info:
            cpac_prov = json_info['CpacProvenance']
        current_prov_list = list(cpac_prov)
        new_prov_list = list(cpac_prov)   # <---- making a copy, it was already a list
        if not inject:
            new_prov_list.append(f'{resource}:{node_name}')
        try:
            res, new_pipe_idx = self.generate_prov_string(new_prov_list)
        except IndexError:
            raise IndexError(f'\n\nThe set_data() call for {resource} has no '
                             'provenance information and should not be an '
                             'injection.')
        if not json_info:
            json_info = {'RawSources': [resource]}     # <---- this will be repopulated to the full file path at the end of the pipeline building, in gather_pipes()
        json_info['CpacProvenance'] = new_prov_list

        if resource not in self.rpool.keys():
            self.rpool[resource] = {}
        else:
            if not fork:     # <--- in the event of multiple strategies/options, this will run for every option; just keep in mind
                search = False
                if self.get_resource_from_prov(current_prov_list) == resource:
                    pipe_idx = self.generate_prov_string(current_prov_list)[1] # CHANGING PIPE_IDX, BE CAREFUL DOWNSTREAM IN THIS FUNCTION
                    if pipe_idx not in self.rpool[resource].keys():
                        search = True
                else:
                    search = True
                if search:
                    for idx in current_prov_list:
                        if self.get_resource_from_prov(idx) == resource:
                            if isinstance(idx, list):
                                pipe_idx = self.generate_prov_string(idx)[1] # CHANGING PIPE_IDX, BE CAREFUL DOWNSTREAM IN THIS FUNCTION
                            elif isinstance(idx, str):
                                pipe_idx = idx
                            break
                if pipe_idx in self.rpool[resource].keys():  # <--- in case the resource name is now new, and not the original
                    del self.rpool[resource][pipe_idx]  # <--- remove old keys so we don't end up with a new strat for every new node unit (unless we fork)
        if new_pipe_idx not in self.rpool[resource]:
            self.rpool[resource][new_pipe_idx] = {}
        if new_pipe_idx not in self.pipe_list:
            self.pipe_list.append(new_pipe_idx)

        self.rpool[resource][new_pipe_idx]['data'] = (node, output)
        self.rpool[resource][new_pipe_idx]['json'] = json_info

    def get(self, resource: LIST_OR_STR, pipe_idx: Optional[str] = None,
            report_fetched: Optional[bool] = False,
            optional: Optional[bool] = False) -> Union[
                TUPLE[Optional[dict], Optional[str]], Optional[dict]]:
        # NOTE!!!
        #   if this is the main rpool, this will return a dictionary of strats, and inside those, are dictionaries like {'data': (node, out), 'json': info}
        #   BUT, if this is a sub rpool (i.e. a strat_pool), this will return a one-level dictionary of {'data': (node, out), 'json': info} WITHOUT THE LEVEL OF STRAT KEYS ABOVE IT
        if not isinstance(resource, list):
            resource = [resource]
        # if a list of potential inputs are given, pick the first one
        # found
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
        raise LookupError(
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
            "https://fcp-indi.github.io/\n")

    def get_data(self, resource, pipe_idx=None, report_fetched=False,
                 quick_single=False):
        if report_fetched:
            if pipe_idx:
                connect, fetched = self.get(resource, pipe_idx=pipe_idx,
                                            report_fetched=report_fetched)
                return (connect['data'], fetched)
            connect, fetched =self.get(resource,
                                       report_fetched=report_fetched)
            return (connect['data'], fetched)
        elif pipe_idx:
            return self.get(resource, pipe_idx=pipe_idx)['data']
        elif quick_single or len(self.get(resource)) == 1:
            for key, val in self.get(resource).items():
                return val['data']
        return self.get(resource)['data']

    def copy_resource(self, resource, new_name):
        try:
            self.rpool[new_name] = self.rpool[resource]
        except KeyError:
            raise Exception(f"[!] {resource} not in the resource pool.")

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
        if 'json' in resource_strat_dct:
            strat_json = resource_strat_dct['json']
        else:
            raise Exception('\n[!] Developer info: the JSON '
                            f'information for {resource} and {strat} '
                            f'is incomplete.\n')
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
        return json_data['CpacProvenance']

    @staticmethod
    def generate_prov_string(prov):
        # this will generate a string from a SINGLE RESOURCE'S dictionary of
        # MULTIPLE PRECEDING RESOURCES (or single, if just one)
        #   NOTE: this DOES NOT merge multiple resources!!! (i.e. for merging-strat pipe_idx generation)
        if not isinstance(prov, list):
            raise Exception('\n[!] Developer info: the CpacProvenance '
                            f'entry for {prov} has to be a list.\n')
        last_entry = get_last_prov_entry(prov)
        resource = last_entry.split(':')[0]
        return (resource, str(prov))

    @staticmethod
    def generate_prov_list(prov_str):
        if not isinstance(prov_str, str):
            raise Exception('\n[!] Developer info: the CpacProvenance '
                            f'entry for {str(prov_str)} has to be a string.\n')
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
            resource = prov.split(':')[0]
            resource_strat_dct[resource] = prov
        else:
            for spot, entry in enumerate(prov):
                if isinstance(entry, list):
                    resource = entry[-1].split(':')[0]
                    resource_strat_dct[resource] = entry
                elif isinstance(entry, str):
                    resource = entry.split(':')[0]
                    resource_strat_dct[resource] = entry
        return resource_strat_dct

    def flatten_prov(self, prov):
        if isinstance(prov, str):
            return [prov]
        elif isinstance(prov, list):
            flat_prov = []
            for entry in prov:
                if isinstance(entry, list):
                    flat_prov += self.flatten_prov(entry)
                else:
                    flat_prov.append(entry)
            return flat_prov

    def get_strats(self, resources, debug=False):

        # TODO: NOTE: NOT COMPATIBLE WITH SUB-RPOOL/STRAT_POOLS
        # TODO: (and it doesn't have to be)

        import itertools

        linked_resources = []
        resource_list = []
        if debug:
            verbose_logger = getLogger('engine')
            verbose_logger.debug('\nresources: %s', resources)
        for resource in resources:
            # grab the linked-input tuples
            if isinstance(resource, tuple):
                linked = []
                for label in list(resource):
                    rp_dct, fetched_resource = self.get(label,
                                                        report_fetched=True,
                                                        optional=True)
                    if not rp_dct:
                        continue
                    linked.append(fetched_resource)
                resource_list += linked
                if len(linked) < 2:
                    continue
                linked_resources.append(linked)
            else:
                resource_list.append(resource)

        total_pool = []
        variant_pool = {}
        len_inputs = len(resource_list)
        if debug:
            verbose_logger = getLogger('engine')
            verbose_logger.debug('linked_resources: %s',
                                 linked_resources)
            verbose_logger.debug('resource_list: %s', resource_list)
        for resource in resource_list:
            rp_dct, fetched_resource = self.get(resource,
                                                report_fetched=True,             # <---- rp_dct has the strats/pipe_idxs as the keys on first level, then 'data' and 'json' on each strat level underneath
                                                optional=True)                   # oh, and we make the resource fetching in get_strats optional so we can have optional inputs, but they won't be optional in the node block unless we want them to be
            if not rp_dct:
                len_inputs -= 1
                continue
            sub_pool = []
            if debug:
                verbose_logger.debug('len(rp_dct): %s\n', len(rp_dct))
            for strat in rp_dct.keys():
                json_info = self.get_json(fetched_resource, strat)
                cpac_prov = json_info['CpacProvenance']
                sub_pool.append(cpac_prov)
                if fetched_resource not in variant_pool:
                    variant_pool[fetched_resource] = []
                if 'CpacVariant' in json_info:
                    for key, val in json_info['CpacVariant'].items():
                        if val not in variant_pool[fetched_resource]:
                            variant_pool[fetched_resource] += val
                            variant_pool[fetched_resource].append(
                                f'NO-{val[0]}')

            if debug:
                verbose_logger = getLogger('engine')
                verbose_logger.debug('%s sub_pool: %s\n', resource, sub_pool)
            total_pool.append(sub_pool)

        if not total_pool:
            raise LookupError('\n\n[!] C-PAC says: None of the listed '
                              'resources in the node block being connected '
                              'exist in the resource pool.\n\nResources:\n'
                              '%s\n\n' % resource_list)

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
                verbose_logger = getLogger('engine')
                verbose_logger.debug('len(strat_list_list): %s\n',
                                     len(strat_list_list))
            for strat_list in strat_list_list:

                json_dct = {}
                for strat in strat_list:
                    # strat is a prov list for a single resource/input
                    strat_resource, strat_idx = \
                        self.generate_prov_string(strat)
                    strat_json = self.get_json(strat_resource,
                                               strat=strat_idx)
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
                                
                                if 'CpacVariant' not in xjson:
                                    xjson['CpacVariant'] = {}
                                if 'CpacVariant' not in yjson:
                                    yjson['CpacVariant'] = {}
                                    
                                current_strat = []
                                for key, val in xjson['CpacVariant'].items():
                                    if isinstance(val, list):
                                        current_strat.append(val[0])
                                    else:
                                        current_strat.append(val)
                                current_spread = list(set(variant_pool[xlabel]))
                                for spread_label in current_spread:
                                    if 'NO-' in spread_label:
                                        continue
                                    if spread_label not in current_strat:
                                        current_strat.append(f'NO-{spread_label}')
                                
                                other_strat = []
                                for key, val in yjson['CpacVariant'].items():
                                    if isinstance(val, list):
                                        other_strat.append(val[0])
                                    else:
                                        other_strat.append(val)
                                other_spread = list(set(variant_pool[ylabel]))
                                for spread_label in other_spread:
                                    if 'NO-' in spread_label:
                                        continue
                                    if spread_label not in other_strat:
                                        other_strat.append(f'NO-{spread_label}')
                                
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
                new_strats[pipe_idx] = ResourcePool()     # <----- new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                
                # placing JSON info at one level higher only for copy convenience
                new_strats[pipe_idx].rpool['json'] = {}
                new_strats[pipe_idx].rpool['json']['subjson'] = {}
                new_strats[pipe_idx].rpool['json']['CpacProvenance'] = strat_list

                # now just invert resource:strat to strat:resource for each resource:strat
                for cpac_prov in strat_list:
                    resource, strat = self.generate_prov_string(cpac_prov)
                    resource_strat_dct = self.rpool[resource][strat]   # <----- remember, this is the dct of 'data' and 'json'.
                    new_strats[pipe_idx].rpool[resource] = resource_strat_dct   # <----- new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS! each one is a new slice of the resource pool combined together.
                    self.pipe_list.append(pipe_idx)
                    if 'CpacVariant' in resource_strat_dct['json']:
                        if 'CpacVariant' not in new_strats[pipe_idx].rpool['json']:
                            new_strats[pipe_idx].rpool['json']['CpacVariant'] = {}
                        for younger_resource, variant_list in resource_strat_dct['json']['CpacVariant'].items():
                            if younger_resource not in new_strats[pipe_idx].rpool['json']['CpacVariant']:
                                new_strats[pipe_idx].rpool['json']['CpacVariant'][younger_resource] = variant_list
                    # preserve each input's JSON info also
                    data_type = resource.split('_')[-1]
                    if data_type not in new_strats[pipe_idx].rpool['json']['subjson']:
                        new_strats[pipe_idx].rpool['json']['subjson'][data_type] = {}
                    new_strats[pipe_idx].rpool['json']['subjson'][data_type].update(copy.deepcopy(resource_strat_dct['json']))
        else:
            new_strats = {}
            for resource_strat_list in total_pool:       # total_pool will have only one list of strats, for the one input
                for cpac_prov in resource_strat_list:     # <------- cpac_prov here doesn't need to be modified, because it's not merging with other inputs
                    resource, pipe_idx = self.generate_prov_string(cpac_prov)
                    resource_strat_dct = self.rpool[resource][pipe_idx]   # <----- remember, this is the dct of 'data' and 'json'.
                    new_strats[pipe_idx] = ResourcePool(rpool={resource: resource_strat_dct})   # <----- again, new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                    # placing JSON info at one level higher only for copy convenience
                    new_strats[pipe_idx].rpool['json'] = resource_strat_dct['json']  # TODO: WARNING- THIS IS A LEVEL HIGHER THAN THE ORIGINAL 'JSON' FOR EASE OF ACCESS IN CONNECT_BLOCK WITH THE .GET(JSON)
                    new_strats[pipe_idx].rpool['json']['subjson'] = {}
                    new_strats[pipe_idx].rpool['json']['CpacProvenance'] = cpac_prov
                    # preserve each input's JSON info also
                    data_type = resource.split('_')[-1]                    
                    if data_type not in new_strats[pipe_idx].rpool['json']['subjson']:
                        new_strats[pipe_idx].rpool['json']['subjson'][data_type] = {}
                    new_strats[pipe_idx].rpool['json']['subjson'][data_type].update(copy.deepcopy(resource_strat_dct['json']))
        return new_strats

    def derivative_xfm(self, wf, label, connection, json_info, pipe_idx,
                       pipe_x):

        if label in self.xfm:

            json_info = dict(json_info)

            # get the bold-to-template transform from the current strat_pool
            # info
            xfm_idx = None
            xfm_label = 'from-bold_to-template_mode-image_xfm'
            for entry in json_info['CpacProvenance']:
                if isinstance(entry, list):
                    if entry[-1].split(':')[0] == xfm_label:
                        xfm_prov = entry
                        xfm_idx = self.generate_prov_string(xfm_prov)[1]
                        break

            # but if the resource doesn't have the bold-to-template transform
            # in its provenance/strategy, find the appropriate one for this
            # current pipe_idx/strat
            if not xfm_idx:
                xfm_info = []
                for pipe_idx, entry in self.get(xfm_label).items():
                    xfm_info.append((pipe_idx, entry['json']['CpacProvenance']))
            else:
                xfm_info = [(xfm_idx, xfm_prov)]

            for num, xfm_entry in enumerate(xfm_info):

                xfm_idx, xfm_prov = xfm_entry
                reg_tool = check_prov_for_regtool(xfm_prov)

                xfm = transform_derivative(f'{label}_xfm_{pipe_x}_{num}',
                                           label, reg_tool, self.num_cpus,
                                           self.num_ants_cores,
                                           ants_interp=self.ants_interp,
                                           fsl_interp=self.fsl_interp,
                                           opt=None)
                wf.connect(connection[0], connection[1],
                           xfm, 'inputspec.in_file')

                node, out = self.get_data("T1w-brain-template-deriv",
                                          quick_single=True)
                wf.connect(node, out, xfm, 'inputspec.reference')

                node, out = self.get_data('from-bold_to-template_mode-image_xfm',
                                          pipe_idx=xfm_idx)
                wf.connect(node, out, xfm, 'inputspec.transform')

                label = f'space-template_{label}'
                json_info['Template'] = self.get_json_info('T1w-brain-template-deriv',
                                                           None, 'Description')
                new_prov = json_info['CpacProvenance'] + xfm_prov
                json_info['CpacProvenance'] = new_prov
                new_pipe_idx = self.generate_prov_string(new_prov)
                self.set_data(label, xfm, 'outputspec.out_file', json_info,
                              new_pipe_idx, f'{label}_xfm_{num}', fork=True)

        return wf

    @property
    def filtered_movement(self) -> bool:
        """
        Check if the movement parameters have been filtered in this strat_pool

        Returns
        -------
        bool
        """
        try:
            return 'motion_estimate_filter' in str(self.get_cpac_provenance(
                'desc-movementParameters_motion'))
        except KeyError:
            # not a strat_pool or no movement parameters in strat_pool
            return False

    def filter_name(self, cfg) -> str:
        """
        In a strat_pool with filtered movement parameters, return the
        name of the filter for this strategy

        Returns
        -------
        str
        """
        motion_filters = cfg['functional_preproc',
                             'motion_estimates_and_correction',
                             'motion_estimate_filter', 'filters']
        if len(motion_filters) == 1 and cfg.switch_is_on([
            'functional_preproc', 'motion_estimates_and_correction',
            'motion_estimate_filter', 'run'], exclusive=True
        ):
            return motion_filters[0]['Name']
        try:
            key = 'motion'
            sidecar = self.get_json('desc-movementParameters_motion')
        except KeyError:
            sidecar = None
        if sidecar is not None and 'CpacVariant' in sidecar:
            if sidecar['CpacVariant'][key]:
                return sidecar['CpacVariant'][key][0][::-1].split('_',
                                                                  1)[0][::-1]
        return 'none'

    def post_process(self, wf, label, connection, json_info, pipe_idx, pipe_x,
                     outs):

        input_type = 'func_derivative'

        post_labels = [(label, connection[0], connection[1])]

        if re.match(r'(.*_)?[ed]c[bw]$', label) or re.match(r'(.*_)?lfcd[bw]$',
                                                            label):
            # suffix: [eigenvector or degree] centrality [binarized or weighted]
            # or lfcd [binarized or weighted]
            mask = 'template-specification-file'
        elif 'space-template' in label:
            if 'space-template_res-derivative_desc-bold_mask' in self.rpool.keys():
                mask = 'space-template_res-derivative_desc-bold_mask'
            else:
                mask = 'space-template_desc-bold_mask'
        else:
            mask = 'space-bold_desc-brain_mask'

        mask_idx = None
        for entry in json_info['CpacProvenance']:
            if isinstance(entry, list):
                if entry[-1].split(':')[0] == mask:
                    mask_prov = entry
                    mask_idx = self.generate_prov_string(mask_prov)[1]
                    break

        if self.smoothing_bool:
            if label in Outputs.to_smooth:
                for smooth_opt in self.smooth_opts:

                    sm = spatial_smoothing(f'{label}_smooth_{smooth_opt}_'
                                           f'{pipe_x}',
                                           self.fwhm, input_type, smooth_opt)
                    wf.connect(connection[0], connection[1],
                               sm, 'inputspec.in_file')
                    node, out = self.get_data(mask, pipe_idx=mask_idx,
                                              quick_single=mask_idx is None)
                    wf.connect(node, out, sm, 'inputspec.mask')

                    if 'desc-' not in label:
                        if 'space-' in label:
                            for tag in label.split('_'):
                                if 'space-' in tag:
                                    smlabel = label.replace(tag,
                                                            f'{tag}_desc-sm')
                                    break
                        else:
                            smlabel = f'desc-sm_{label}'
                    else:
                        for tag in label.split('_'):
                            if 'desc-' in tag:
                                newtag = f'{tag}-sm'
                                smlabel = label.replace(tag, newtag)
                                break

                    post_labels.append((smlabel, sm, 'outputspec.out_file'))

                    self.set_data(smlabel, sm, 'outputspec.out_file',
                                  json_info, pipe_idx,
                                  f'spatial_smoothing_{smooth_opt}',
                                  fork=True)
                    self.set_data('fwhm', sm, 'outputspec.fwhm', json_info,
                                  pipe_idx, f'spatial_smoothing_{smooth_opt}',
                                  fork=True)

        if self.zscoring_bool:            
            for label_con_tpl in post_labels:
                label = label_con_tpl[0]
                connection = (label_con_tpl[1], label_con_tpl[2])
                if label in Outputs.to_zstd:
                    zstd = z_score_standardize(f'{label}_zstd_{pipe_x}',
                                               input_type)

                    wf.connect(connection[0], connection[1],
                               zstd, 'inputspec.in_file')

                    node, out = self.get_data(mask, pipe_idx=mask_idx)
                    wf.connect(node, out, zstd, 'inputspec.mask')

                    if 'desc-' not in label:
                        if 'space-template' in label:
                            new_label = label.replace('space-template',
                                                      'space-template_desc-zstd')
                        else:
                            new_label = f'desc-zstd_{label}'
                    else:
                        for tag in label.split('_'):
                            if 'desc-' in tag:
                                newtag = f'{tag}-zstd'
                                new_label = label.replace(tag, newtag)
                                break

                    post_labels.append((new_label, zstd, 'outputspec.out_file'))

                    self.set_data(new_label, zstd, 'outputspec.out_file',
                                  json_info, pipe_idx, f'zscore_standardize',
                                  fork=True)

                elif label in Outputs.to_fisherz:

                    zstd = fisher_z_score_standardize(f'{label}_zstd_{pipe_x}',
                                                      label, input_type)

                    wf.connect(connection[0], connection[1],
                               zstd, 'inputspec.correlation_file')

                    # if the output is 'space-template_desc-MeanSCA_correlations', we want
                    # 'desc-MeanSCA_timeseries'
                    oned = label.replace('correlations', 'timeseries')

                    node, out = outs[oned]
                    wf.connect(node, out, zstd, 'inputspec.timeseries_oned')

                    post_labels.append((new_label, zstd, 'outputspec.out_file'))

                    self.set_data(new_label, zstd, 'outputspec.out_file',
                                  json_info, pipe_idx,
                                  'fisher_zscore_standardize',
                                  fork=True)

        return (wf, post_labels)

    def gather_pipes(self, wf, cfg, all=False, add_incl=None, add_excl=None):
        excl = []
        substring_excl = []
        outputs_logger = getLogger(f'{cfg["subject_id"]}_expectedOutputs')
        expected_outputs = ExpectedOutputs()

        if add_excl:
            excl += add_excl

        if 'nonsmoothed' not in cfg.post_processing['spatial_smoothing'][
                'output']:
            excl += Outputs.native_nonsmooth
            excl += Outputs.template_nonsmooth

        if 'raw' not in cfg.post_processing['z-scoring']['output']:
            excl += Outputs.native_raw
            excl += Outputs.template_raw

        if not cfg.pipeline_setup['output_directory']['write_debugging_outputs']:
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

            subdir = 'other'
            if resource in Outputs.anat:
                subdir = 'anat'
                #TODO: get acq- etc.
            elif resource in Outputs.func:
                subdir = 'func'
                #TODO: other stuff like acq- etc.

            for pipe_idx in self.rpool[resource]:
                unique_id = self.get_name()
                part_id = unique_id.split('_')[0]
                ses_id = unique_id.split('_')[1]

                if 'ses-' not in ses_id:
                    ses_id = f"ses-{ses_id}"

                out_dir = cfg.pipeline_setup['output_directory']['path']
                pipe_name = cfg.pipeline_setup['pipeline_name']
                container = os.path.join(f'pipeline_{pipe_name}', part_id,
                                         ses_id)
                filename = f'{unique_id}_{res_in_filename(self.cfg, resource)}'

                out_path = os.path.join(out_dir, container, subdir, filename)

                out_dct = {
                    'unique_id': unique_id,
                    'out_dir': out_dir,
                    'container': container,
                    'subdir': subdir,
                    'filename': filename,
                    'out_path': out_path
                }
                self.rpool[resource][pipe_idx]['out'] = out_dct

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
            all_jsons = [self.rpool[resource][pipe_idx]['json'] for pipe_idx in
                         self.rpool[resource]]
            unlabelled = set(key for json_info in all_jsons for key in
                             json_info.get('CpacVariant', {}).keys() if
                             key not in (*MOVEMENT_FILTER_KEYS, 'timeseries'))
            if 'bold' in unlabelled:
                all_bolds = list(
                    chain.from_iterable(json_info['CpacVariant']['bold'] for
                                        json_info in all_jsons if
                                        'CpacVariant' in json_info and
                                        'bold' in json_info['CpacVariant']))
                # not any(not) because all is overloaded as a parameter here
                if not any(not re.match(r'apply_(phasediff|blip)_to_'
                                        r'timeseries_separately_.*', _bold)
                           for _bold in all_bolds):
                    # this fork point should only result in 0 or 1 forks
                    unlabelled.remove('bold')
                del all_bolds
            all_forks = {key: set(
                chain.from_iterable(json_info['CpacVariant'][key] for
                                    json_info in all_jsons if
                                    'CpacVariant' in json_info and
                                    key in json_info['CpacVariant'])) for
                key in unlabelled}
            # del all_jsons
            for key, forks in all_forks.items():
                if len(forks) < 2:  # no int suffix needed if only one fork
                    unlabelled.remove(key)
            # del all_forks
            for pipe_idx in self.rpool[resource]:
                pipe_x = self.get_pipe_number(pipe_idx)
                json_info = self.rpool[resource][pipe_idx]['json']
                out_dct = self.rpool[resource][pipe_idx]['out']

                try:
                    if unlabelled:
                        num_variant += 1
                except TypeError:
                    pass

                try:
                    del json_info['subjson']
                except KeyError:
                    pass

                if out_dct['subdir'] == 'other' and not all:
                    continue

                unique_id = out_dct['unique_id']
                resource_idx = resource

                if isinstance(num_variant, int):
                    resource_idx, out_dct = name_fork(resource_idx, cfg,
                                                      json_info, out_dct)
                    if unlabelled:
                        if 'desc-' in out_dct['filename']:
                            for key in out_dct['filename'].split('_')[::-1]:
                                # final `desc` entity
                                if key.startswith('desc-'):
                                    out_dct['filename'] = out_dct['filename'
                                                                  ].replace(
                                        key, f'{key}-{num_variant}')
                                    resource_idx = resource_idx.replace(
                                        key, f'{key}-{num_variant}')
                                    break
                        else:
                            suff = resource.split('_')[-1]
                            newdesc_suff = f'desc-{num_variant}_{suff}'
                            resource_idx = resource_idx.replace(suff,
                                                                newdesc_suff)
                id_string = pe.Node(Function(input_names=['cfg', 'unique_id',
                                                          'resource',
                                                          'scan_id',
                                                          'template_desc',
                                                          'atlas_id',
                                                          'fwhm',
                                                          'subdir',
                                                          'extension'],
                                             output_names=['out_filename'],
                                             function=create_id_string),
                                    name=f'id_string_{resource_idx}_{pipe_x}')
                id_string.inputs.cfg = self.cfg
                id_string.inputs.unique_id = unique_id
                id_string.inputs.resource = resource_idx
                id_string.inputs.subdir = out_dct['subdir']

                # grab the iterable scan ID
                if out_dct['subdir'] == 'func':
                    node, out = self.rpool['scan']["['scan:func_ingress']"][
                            'data']
                    wf.connect(node, out, id_string, 'scan_id')
                
                self.back_propogate_template_name(wf, resource_idx, json_info,
                                                  id_string)
                # grab the FWHM if smoothed
                for tag in resource.split('_'):
                    if 'desc-' in tag and '-sm' in tag:
                        fwhm_idx = pipe_idx.replace(f'{resource}:', 'fwhm:')
                        try:
                            node, out = self.rpool['fwhm'][fwhm_idx]['data']
                            wf.connect(node, out, id_string, 'fwhm')
                        except KeyError:
                            # smoothing was not done for this resource in the
                            # engine.py smoothing
                            pass
                        break
                atlas_suffixes = ['timeseries', 'correlations', 'statmap']
                # grab the iterable atlas ID
                atlas_id = None
                if not resource.endswith('desc-confounds_timeseries'):
                    if resource.split('_')[-1] in atlas_suffixes:
                        atlas_idx = pipe_idx.replace(resource, 'atlas_name')
                        # need the single quote and the colon inside the double
                        # quotes - it's the encoded pipe_idx
                        #atlas_idx = new_idx.replace(f"'{temp_rsc}:",
                        #                            "'atlas_name:")
                        if atlas_idx in self.rpool['atlas_name']:
                            node, out = self.rpool['atlas_name'][atlas_idx][
                                'data']
                            wf.connect(node, out, id_string, 'atlas_id')
                        elif 'atlas-' in resource:
                            for tag in resource.split('_'):
                                if 'atlas-' in tag:
                                    atlas_id = tag.replace('atlas-', '')
                            id_string.inputs.atlas_id = atlas_id
                        else:
                            warnings.warn(str(
                                LookupError("\n[!] No atlas ID found for "
                                        f"{out_dct['filename']}.\n")))
                nii_name = pe.Node(Rename(), name=f'nii_{resource_idx}_'
                                                  f'{pipe_x}')
                nii_name.inputs.keep_ext = True
                
                if resource in Outputs.ciftis:
                   nii_name.inputs.keep_ext = False
                   id_string.inputs.extension = Outputs.ciftis[resource]
                else:
                   nii_name.inputs.keep_ext = True
                
               
                if resource in Outputs.giftis:

                   nii_name.inputs.keep_ext = False
                   id_string.inputs.extension = f'{Outputs.giftis[resource]}.gii'
                   
                else:
                   nii_name.inputs.keep_ext = True
                
                wf.connect(id_string, 'out_filename',
                           nii_name, 'format_string')
                
                node, out = self.rpool[resource][pipe_idx]['data']
                try:
                    wf.connect(node, out, nii_name, 'in_file')
                except OSError as os_error:
                    logger.warning(os_error)
                    continue

                write_json_imports = ['import os', 'import json']
                write_json = pe.Node(Function(input_names=['json_data',
                                                           'filename'],
                                              output_names=['json_file'],
                                              function=write_output_json,
                                              imports=write_json_imports),
                                     name=f'json_{resource_idx}_{pipe_x}')
                write_json.inputs.json_data = json_info

                wf.connect(id_string, 'out_filename', write_json, 'filename')
                ds = pe.Node(DataSink(), name=f'sinker_{resource_idx}_'
                                              f'{pipe_x}')
                ds.inputs.parameterization = False
                ds.inputs.base_directory = out_dct['out_dir']
                ds.inputs.encrypt_bucket_keys = cfg.pipeline_setup[
                    'Amazon-AWS']['s3_encryption']
                ds.inputs.container = out_dct['container']

                if cfg.pipeline_setup['Amazon-AWS'][
                    'aws_output_bucket_credentials']:
                    ds.inputs.creds_path = cfg.pipeline_setup['Amazon-AWS'][
                        'aws_output_bucket_credentials']
                expected_outputs += (out_dct['subdir'], create_id_string(
                    self.cfg, unique_id, resource_idx,
                    template_desc=id_string.inputs.template_desc,
                    atlas_id=atlas_id, subdir=out_dct['subdir']))
                wf.connect(nii_name, 'out_file',
                           ds, f'{out_dct["subdir"]}.@data')
                wf.connect(write_json, 'json_file',
                           ds, f'{out_dct["subdir"]}.@json')
        outputs_logger.info(expected_outputs)

    def node_data(self, resource, **kwargs):
        '''Factory function to create NodeData objects

        Parameters
        ----------
        resource : str

        Returns
        -------
        NodeData
        '''
        return NodeData(self, resource, **kwargs)


class NodeBlock:
    def __init__(self, node_block_functions, debug=False):
        if not isinstance(node_block_functions, list):
            node_block_functions = [node_block_functions]

        self.node_blocks = {}

        for node_block_function in node_block_functions:    # <---- sets up the NodeBlock object in case you gave it a list of node blocks instead of a single one - for option forking.
        
            self.input_interface = []
            if isinstance(node_block_function, tuple):
                self.input_interface = node_block_function[1]
                node_block_function = node_block_function[0]
                if not isinstance(self.input_interface, list):
                    self.input_interface = [self.input_interface]

            if not isinstance(node_block_function, NodeBlockFunction):
                # If the object is a plain function `__name__` will be more useful then `str()`
                obj_str = node_block_function.__name__ \
                    if hasattr(node_block_function, '__name__') else \
                    str(node_block_function)
                raise TypeError(f'Object is not a nodeblock: "{obj_str}"')

            name = node_block_function.name
            self.name = name
            self.node_blocks[name] = {}

            if self.input_interface:
                for interface in self.input_interface:
                    for orig_input in node_block_function.inputs:
                        if isinstance(orig_input, tuple):
                            list_tup = list(orig_input)
                            if interface[0] in list_tup:
                                list_tup.remove(interface[0])
                                list_tup.append(interface[1])
                                node_block_function.inputs.remove(orig_input)
                                node_block_function.inputs.append(tuple(list_tup))
                        else:
                            if orig_input == interface[0]:
                                node_block_function.inputs.remove(interface[0])
                                node_block_function.inputs.append(interface[1])

            for key, val in node_block_function.legacy_nodeblock_dict().items():
                self.node_blocks[name][key] = val

            self.node_blocks[name]['block_function'] = node_block_function

            #TODO: fix/replace below
            self.outputs = {}
            for out in node_block_function.outputs:
                self.outputs[out] = None

            self.options = ['base']
            if node_block_function.outputs is not None:
                self.options = node_block_function.outputs

            logger.info('Connecting %s...', name)
            if debug:
                config.update_config(
                    {'logging': {'workflow_level': 'DEBUG'}})
                logging.update_logging(config)
                logger.debug('"inputs": %s\n\t "outputs": %s%s',
                             node_block_function.inputs,
                             list(self.outputs.keys()),
                             f'\n\t"options": {self.options}'
                             if self.options != ['base'] else '')
                config.update_config(
                    {'logging': {'workflow_level': 'INFO'}})
                logging.update_logging(config)

    def get_name(self):
        return self.name

    def check_null(self, val):
        if isinstance(val, str):
            val = None if val.lower() == 'none' else val
        return val

    def check_output(self, outputs, label, name):
        if label not in outputs:
            raise NameError(f'\n[!] Output name "{label}" in the block '
                            'function does not match the outputs list '
                            f'{outputs} in Node Block "{name}"\n')

    def grab_tiered_dct(self, cfg, key_list):
        cfg_dct = cfg.dict()
        for key in key_list:
            try:
                cfg_dct = cfg_dct.get(key, {})
            except KeyError:
                raise Exception(f"[!] The config provided to the node block is not valid")  
        return cfg_dct

    def connect_block(self, wf, cfg, rpool):
        debug = cfg.pipeline_setup['Debugging']['verbose']
        all_opts = []
        for name, block_dct in self.node_blocks.items():
            opts = []
            config = self.check_null(block_dct['config'])
            option_key = self.check_null(block_dct['option_key'])
            option_val = self.check_null(block_dct['option_val'])
            if option_key and option_val:
                if not isinstance(option_key, list):
                    option_key = [option_key]
                if not isinstance(option_val, list):
                    option_val = [option_val]
                if config:
                    key_list = config + option_key
                else:
                    key_list = option_key
                if 'USER-DEFINED' in option_val:
                    # load custom config data into each 'opt'
                    opts = self.grab_tiered_dct(cfg, key_list)
                else:
                    for option in option_val:
                        try:
                            if option in self.grab_tiered_dct(cfg, key_list):   # <---- goes over the option_vals in the node block docstring, and checks if the user's pipeline config included it in the forking list
                                opts.append(option)
                        except AttributeError as err:
                            raise Exception(f"{err}\nNode Block: {name}")

                if opts is None:
                    opts = [opts]

            elif option_key and not option_val:
                # enables multiple config forking entries
                if not isinstance(option_key[0], list):
                    raise Exception(f'[!] The option_key field ({option_key}) '
                                    f'for {name} exists but there is no '
                                    'option_val.\n\nIf you are trying to '
                                    'populate multiple option keys, the '
                                    'option_val field must contain a list of '
                                    'a list.\n')
                for option_config in option_key:
                    # option_config is a list of pipe config levels down to the option
                    if config:
                        key_list = config + option_config
                    else:
                        key_list = option_config
                    option_val = option_config[-1]
                    if option_val in self.grab_tiered_dct(cfg, key_list[:-1]):
                        opts.append(option_val)                
            else:                                           #         AND, if there are multiple option-val's (in a list) in the docstring, it gets iterated below in 'for opt in option' etc. AND THAT'S WHEN YOU HAVE TO DELINEATE WITHIN THE NODE BLOCK CODE!!!
                opts = [None]
            all_opts += opts

        sidecar_additions = {
            'CpacConfigHash': hashlib.sha1(json.dumps(cfg.dict(), sort_keys=True).encode('utf-8')).hexdigest(),
            'CpacConfig': cfg.dict()
        }

        if cfg['pipeline_setup']['output_directory'].get('user_defined'):
            sidecar_additions['UserDefined'] = cfg['pipeline_setup']['output_directory']['user_defined']

        for name, block_dct in self.node_blocks.items():    # <--- iterates over either the single node block in the sequence, or a list of node blocks within the list of node blocks, i.e. for option forking.

            switch = self.check_null(block_dct['switch'])
            config = self.check_null(block_dct['config'])
            option_key = self.check_null(block_dct['option_key'])
            option_val = self.check_null(block_dct['option_val'])
            inputs = self.check_null(block_dct['inputs'])
            outputs = self.check_null(block_dct['outputs'])

            block_function = block_dct['block_function']

            opts = []
            if option_key and option_val:
                if not isinstance(option_key, list):
                    option_key = [option_key]
                if not isinstance(option_val, list):
                    option_val = [option_val]
                if config:
                    key_list = config + option_key
                else:
                    key_list = option_key
                if 'USER-DEFINED' in option_val:
                    # load custom config data into each 'opt'
                    opts = self.grab_tiered_dct(cfg, key_list)
                else:
                    for option in option_val:
                        if option in self.grab_tiered_dct(cfg, key_list):   # <---- goes over the option_vals in the node block docstring, and checks if the user's pipeline config included it in the forking list
                            opts.append(option)
            else:                                                           #         AND, if there are multiple option-val's (in a list) in the docstring, it gets iterated below in 'for opt in option' etc. AND THAT'S WHEN YOU HAVE TO DELINEATE WITHIN THE NODE BLOCK CODE!!!
                opts = [None]                                               #         THIS ALSO MEANS the multiple option-val's in docstring node blocks can be entered once in the entire node-block sequence, not in a list of multiples
            if not opts:
                # for node blocks where the options are split into different
                # block functions - opts will be empty for non-selected
                # options, and would waste the get_strats effort below
                continue

            if not switch:
                switch = [True]
            else:
                if config:
                    try:
                        key_list = config + switch
                    except TypeError:
                        raise Exception("\n\n[!] Developer info: Docstring error "
                                        f"for {name}, make sure the 'config' or "
                                        "'switch' fields are lists.\n\n")
                    switch = self.grab_tiered_dct(cfg, key_list)
                    
                else:
                    if isinstance(switch[0], list):
                        # we have multiple switches, which is designed to only work if
                        # config is set to "None"
                        switch_list = []
                        for key_list in switch:
                            val = self.grab_tiered_dct(cfg, key_list)
                            if isinstance(val, list):
                                # fork switches
                                if True in val:
                                    switch_list.append(True)
                                if False in val:
                                    switch_list.append(False)
                            else:
                                switch_list.append(val)
                        if False in switch_list:
                            switch = [False]
                        else:
                            switch = [True]
                    else:
                        # if config is set to "None"
                        key_list = switch
                        switch = self.grab_tiered_dct(cfg, key_list)
                if not isinstance(switch, list):
                    switch = [switch]
            if True in switch:
                for pipe_idx, strat_pool in rpool.get_strats(
                        inputs, debug).items():         # strat_pool is a ResourcePool like {'desc-preproc_T1w': { 'json': info, 'data': (node, out) }, 'desc-brain_mask': etc.}
                    fork = False in switch                                            #   keep in mind rpool.get_strats(inputs) = {pipe_idx1: {'desc-preproc_T1w': etc.}, pipe_idx2: {..} }
                    for opt in opts:                                            #   it's a dictionary of ResourcePools called strat_pools, except those sub-ResourcePools only have one level! no pipe_idx strat keys.
                        # remember, you can get 'data' or 'json' from strat_pool with member functions
                        # strat_pool has all of the JSON information of all the inputs!
                        # so when we set_data below for the TOP-LEVEL MAIN RPOOL (not the strat_pool), we can generate new merged JSON information for each output.
                        #    particularly, our custom 'CpacProvenance' field.
                        node_name = name
                        pipe_x = rpool.get_pipe_number(pipe_idx)

                        replaced_inputs = []
                        for interface in self.input_interface:
                            if isinstance(interface[1], list):
                                for input_name in interface[1]:
                                    if strat_pool.check_rpool(input_name):
                                        break
                            else:
                                input_name = interface[1]
                            strat_pool.copy_resource(input_name, interface[0])
                            replaced_inputs.append(interface[0])
                        try:
                            wf, outs = block_function(wf, cfg, strat_pool,
                                                      pipe_x, opt)
                        except IOError as e:  # duplicate node
                            logger.warning(e)
                            continue

                        if not outs:
                            if (block_function.__name__ == 'freesurfer_'
                                                           'postproc'):
                                logger.warning(
                                    WARNING_FREESURFER_OFF_WITH_DATA)
                                LOGTAIL['warnings'].append(
                                    WARNING_FREESURFER_OFF_WITH_DATA)
                            continue

                        if opt and len(option_val) > 1:
                            node_name = f'{node_name}_{opt}'
                        elif opt and 'USER-DEFINED' in option_val:
                            node_name = f'{node_name}_{opt["Name"]}'

                        if debug:
                            verbose_logger = getLogger('engine')
                            verbose_logger.debug('\n=======================')
                            verbose_logger.debug('Node name: %s', node_name)
                            prov_dct = \
                                rpool.get_resource_strats_from_prov(
                                    ast.literal_eval(pipe_idx))
                            for key, val in prov_dct.items():
                                verbose_logger.debug('-------------------')
                                verbose_logger.debug('Input - %s:', key)
                                sub_prov_dct = \
                                    rpool.get_resource_strats_from_prov(val)
                                for sub_key, sub_val in sub_prov_dct.items():
                                    sub_sub_dct = \
                                        rpool.get_resource_strats_from_prov(
                                            sub_val)
                                    verbose_logger.debug('  sub-input - %s:',
                                                         sub_key)
                                    verbose_logger.debug('    prov = %s',
                                                         sub_val)
                                    verbose_logger.debug(
                                        '    sub_sub_inputs = %s',
                                        sub_sub_dct.keys())

                        for label, connection in outs.items():
                            self.check_output(outputs, label, name)
                            new_json_info = copy.deepcopy(strat_pool.get('json'))

                            # transfer over data-specific json info
                            #   for example, if the input data json is _bold and the output is also _bold
                            data_type = label.split('_')[-1]
                            if data_type in new_json_info['subjson']:
                                if 'SkullStripped' in new_json_info['subjson'][data_type]:
                                    new_json_info['SkullStripped'] = new_json_info['subjson'][data_type]['SkullStripped']

                            # determine sources for the outputs, i.e. all input data into the node block                   
                            new_json_info['Sources'] = [x for x in strat_pool.get_entire_rpool() if x != 'json' and x not in replaced_inputs]
                            
                            if isinstance(outputs, dict):
                                new_json_info.update(outputs[label])
                                if 'Description' not in outputs[label]:
                                    # don't propagate old Description
                                    try:
                                        del new_json_info['Description']
                                    except KeyError:
                                        pass
                                if 'Template' in outputs[label]:
                                    template_key = outputs[label]['Template']
                                    if template_key in new_json_info['Sources']:
                                        # only if the pipeline config template key is entered as the 'Template' field
                                        # otherwise, skip this and take in the literal 'Template' string
                                        try:
                                            new_json_info['Template'] = new_json_info['subjson'][template_key]['Description']
                                        except KeyError:
                                            pass
                                    try:
                                        new_json_info['Resolution'] = new_json_info['subjson'][template_key]['Resolution']
                                    except KeyError:
                                        pass
                            else:
                                # don't propagate old Description
                                try:
                                    del new_json_info['Description']
                                except KeyError:
                                    pass

                            if 'Description' in new_json_info:
                                new_json_info['Description'] = ' '.join(new_json_info['Description'].split())

                            for sidecar_key, sidecar_value in sidecar_additions.items():
                                if sidecar_key not in new_json_info:
                                    new_json_info[sidecar_key] = sidecar_value

                            try:
                                del new_json_info['subjson']
                            except KeyError:
                                pass

                            if fork or len(opts) > 1 or len(all_opts) > 1:
                                if 'CpacVariant' not in new_json_info:
                                    new_json_info['CpacVariant'] = {}
                                raw_label = rpool.get_raw_label(label)
                                if raw_label not in new_json_info['CpacVariant']:
                                    new_json_info['CpacVariant'][raw_label] = []
                                new_json_info['CpacVariant'][raw_label].append(node_name)
 
                            rpool.set_data(label,
                                           connection[0],
                                           connection[1],
                                           new_json_info,
                                           pipe_idx, node_name, fork)

                            wf, post_labels = rpool.post_process(
                                wf, label, connection, new_json_info, pipe_idx,
                                pipe_x, outs)

                            if rpool.func_reg:
                                for postlabel in post_labels:
                                    connection = (postlabel[1], postlabel[2])
                                    wf = rpool.derivative_xfm(wf, postlabel[0],
                                                              connection,
                                                              new_json_info,
                                                              pipe_idx,
                                                              pipe_x)
        return wf


def wrap_block(node_blocks, interface, wf, cfg, strat_pool, pipe_num, opt):
    """Wrap a list of node block functions to make them easier to use within
    other node blocks.

    Example usage:

        # This calls the 'bold_mask_afni' and 'bold_masking' node blocks to
        # skull-strip an EPI field map, without having to invoke the NodeBlock
        # connection system.

        # The interface dictionary tells wrap_block to set the EPI field map
        # in the parent node block's throw-away strat_pool as 'bold', so that
        # the 'bold_mask_afni' and 'bold_masking' node blocks will see that as
        # the 'bold' input.

        # It also tells wrap_block to set the 'desc-brain_bold' output of
        # the 'bold_masking' node block to 'opposite_pe_epi_brain' (what it
        # actually is) in the parent node block's strat_pool, which gets
        # returned.

        # Note 'bold' and 'desc-brain_bold' (all on the left side) are the
        # labels that 'bold_mask_afni' and 'bold_masking' understand/expect
        # through their interfaces and docstrings.

        # The right-hand side (the values of the 'interface' dictionary) are
        # what 'make sense' within the current parent node block - in this
        # case, the distortion correction node block dealing with field maps.

        interface = {'bold': (match_epi_fmaps_node, 'opposite_pe_epi'),
                     'desc-brain_bold': 'opposite_pe_epi_brain'}
        wf, strat_pool = wrap_block([bold_mask_afni, bold_masking],
                                    interface, wf, cfg, strat_pool,
                                    pipe_num, opt)

        ...further downstream in the parent node block:

        node, out = strat_pool.get_data('opposite_pe_epi_brain')

        # The above line will connect the output of the 'bold_masking' node
        # block (which is the skull-stripped version of 'opposite_pe_epi') to
        # the next node.

    """
    for block in node_blocks:
        #new_pool = copy.deepcopy(strat_pool)
        for in_resource, val in interface.items():
            if isinstance(val, tuple):
                strat_pool.set_data(in_resource, val[0], val[1], {}, "", "",
                                    fork=True)#
        if 'sub_num' not in strat_pool.get_pool_info():
            strat_pool.set_pool_info({'sub_num': 0})
        sub_num = strat_pool.get_pool_info()['sub_num']
        
        wf, outputs = block(wf, cfg, strat_pool, f'{pipe_num}-{sub_num}', opt)#
        for out, val in outputs.items():
            if out in interface and isinstance(interface[out], str):
                strat_pool.set_data(interface[out], outputs[out][0], outputs[out][1],
                                    {}, "", "")
            else:
                strat_pool.set_data(out, outputs[out][0], outputs[out][1],
                                    {}, "", "")
        sub_num += 1
        strat_pool.set_pool_info({'sub_num': sub_num})

    return (wf, strat_pool)

def ingress_raw_anat_data(wf, rpool, cfg, data_paths, unique_id, part_id,
                          ses_id):

    if 'anat' not in data_paths:
        print('No anatomical data present.')
        return rpool

    if 'creds_path' not in data_paths:
        data_paths['creds_path'] = None

    anat_flow = create_anat_datasource(f'anat_T1w_gather_{part_id}_{ses_id}')

    anat = {}
    if type(data_paths['anat']) is str:
        anat['T1']=data_paths['anat']
    elif 'T1w' in data_paths['anat']:
        anat['T1']=data_paths['anat']['T1w']

    if 'T1' in anat:
        anat_flow.inputs.inputnode.set(
            subject=part_id,
            anat=anat['T1'],
            creds_path=data_paths['creds_path'],
            dl_dir=cfg.pipeline_setup['working_directory']['path'],
            img_type='anat'
        )
        rpool.set_data('T1w', anat_flow, 'outputspec.anat', {},
                    "", "anat_ingress")
    
    if 'T2w' in data_paths['anat']: 
        anat_flow_T2 = create_anat_datasource(f'anat_T2w_gather_{part_id}_{ses_id}')
        anat_flow_T2.inputs.inputnode.set(
            subject=part_id,
            anat=data_paths['anat']['T2w'],
            creds_path=data_paths['creds_path'],
            dl_dir=cfg.pipeline_setup['working_directory']['path'],
            img_type='anat'
        )
        rpool.set_data('T2w', anat_flow_T2, 'outputspec.anat', {},
                    "", "anat_ingress")

    if cfg.surface_analysis['freesurfer']['ingress_reconall']:
        rpool = ingress_freesurfer(wf, rpool, cfg, data_paths, unique_id, part_id,
                          ses_id)
                
    return rpool

def ingress_freesurfer(wf, rpool, cfg, data_paths, unique_id, part_id,
                          ses_id):
    
    try: 
        fs_path = os.path.join(cfg.pipeline_setup['freesurfer_dir'], part_id)
    except KeyError:
        print('No FreeSurfer data present.')
        return rpool
    
    #fs_path = os.path.join(cfg.pipeline_setup['freesurfer_dir'], part_id)
    if not os.path.exists(fs_path):
        if 'sub' in part_id:
            fs_path = os.path.join(cfg.pipeline_setup['freesurfer_dir'], part_id.replace('sub-', ''))
        else:
            fs_path = os.path.join(cfg.pipeline_setup['freesurfer_dir'], ('sub-' + part_id))
        
        # patch for flo-specific data
        if not os.path.exists(fs_path):
            subj_ses = part_id + '-' + ses_id
            fs_path = os.path.join(cfg.pipeline_setup['freesurfer_dir'], subj_ses)
            if not os.path.exists(fs_path):
                print(f'No FreeSurfer data found for subject {part_id}')
                return rpool
    
    # Check for double nested subj names
    if os.path.exists(os.path.join(fs_path, os.path.basename(fs_path))): 
        fs_path = os.path.join(fs_path, part_id)

    fs_ingress = create_general_datasource('gather_freesurfer_dir') 
    fs_ingress.inputs.inputnode.set(
        unique_id=unique_id,
        data=fs_path,
        creds_path=data_paths['creds_path'],
        dl_dir=cfg.pipeline_setup['working_directory']['path'])
    rpool.set_data("freesurfer-subject-dir", fs_ingress, 'outputspec.data',
                    {}, "", "freesurfer_config_ingress")

    recon_outs = {
        'pipeline-fs_raw-average': 'mri/rawavg.mgz',
        'pipeline-fs_subcortical-seg': 'mri/aseg.mgz',
        'pipeline-fs_brainmask': 'mri/brainmask.mgz',
        'pipeline-fs_wmparc': 'mri/wmparc.mgz',
        'pipeline-fs_T1': 'mri/T1.mgz',
        'pipeline-fs_hemi-L_desc-surface_curv': 'surf/lh.curv',
        'pipeline-fs_hemi-R_desc-surface_curv': 'surf/rh.curv',
        'pipeline-fs_hemi-L_desc-surfaceMesh_pial': 'surf/lh.pial',
        'pipeline-fs_hemi-R_desc-surfaceMesh_pial': 'surf/rh.pial',
        'pipeline-fs_hemi-L_desc-surfaceMesh_smoothwm': 'surf/lh.smoothwm',
        'pipeline-fs_hemi-R_desc-surfaceMesh_smoothwm': 'surf/rh.smoothwm',
        'pipeline-fs_hemi-L_desc-surfaceMesh_sphere': 'surf/lh.sphere',
        'pipeline-fs_hemi-R_desc-surfaceMesh_sphere': 'surf/rh.sphere',
        'pipeline-fs_hemi-L_desc-surfaceMap_sulc': 'surf/lh.sulc',
        'pipeline-fs_hemi-R_desc-surfaceMap_sulc': 'surf/rh.sulc',
        'pipeline-fs_hemi-L_desc-surfaceMap_thickness': 'surf/lh.thickness',
        'pipeline-fs_hemi-R_desc-surfaceMap_thickness': 'surf/rh.thickness',
        'pipeline-fs_hemi-L_desc-surfaceMap_volume': 'surf/lh.volume',
        'pipeline-fs_hemi-R_desc-surfaceMap_volume': 'surf/rh.volume',
        'pipeline-fs_hemi-L_desc-surfaceMesh_white': 'surf/lh.white',
        'pipeline-fs_hemi-R_desc-surfaceMesh_white': 'surf/rh.white',
        'pipeline-fs_xfm': 'mri/transforms/talairach.lta'
    }
    
    for key, outfile in recon_outs.items():
        fullpath = os.path.join(fs_path, outfile)
        if os.path.exists(fullpath):
            fs_ingress = create_general_datasource(f'gather_fs_{key}_dir')
            fs_ingress.inputs.inputnode.set(
                unique_id=unique_id,
                data=fullpath,
                creds_path=data_paths['creds_path'],
                dl_dir=cfg.pipeline_setup['working_directory']['path'])
            rpool.set_data(key, fs_ingress, 'outputspec.data',
                            {}, "", f"fs_{key}_ingress")
        else:
            warnings.warn(str(
                    LookupError("\n[!] Path does not exist for "
                                    f"{fullpath}.\n")))
                
    return rpool

def ingress_raw_func_data(wf, rpool, cfg, data_paths, unique_id, part_id,
                          ses_id):

    func_paths_dct = data_paths['func']

    func_wf = create_func_datasource(func_paths_dct, rpool,
                                     f'func_ingress_{part_id}_{ses_id}')
    func_wf.inputs.inputnode.set(
        subject=part_id,
        creds_path=data_paths['creds_path'],
        dl_dir=cfg.pipeline_setup['working_directory']['path']
    )
    func_wf.get_node('inputnode').iterables = \
        ("scan", list(func_paths_dct.keys()))   

    rpool.set_data('subject', func_wf, 'outputspec.subject', {}, "",
                   "func_ingress")
    rpool.set_data('bold', func_wf, 'outputspec.rest', {}, "", "func_ingress")
    rpool.set_data('scan', func_wf, 'outputspec.scan', {}, "", "func_ingress")
    rpool.set_data('scan-params', func_wf, 'outputspec.scan_params', {}, "",
                   "scan_params_ingress")
    
    # TODO: CHECK FOR PARAMETERS

    wf, rpool, diff, blip, fmap_rp_list = \
        ingress_func_metadata(wf, cfg, rpool, data_paths, part_id,
                              data_paths['creds_path'], ses_id)

    # Memoize list of local functional scans
    # TODO: handle S3 files
    # Skip S3 files for now

    local_func_scans = [
        func_paths_dct[scan]['scan'] for scan in func_paths_dct.keys() if not
        func_paths_dct[scan]['scan'].startswith('s3://')]
    if local_func_scans:
        # pylint: disable=protected-access
        wf._local_func_scans = local_func_scans
        if cfg.pipeline_setup['Debugging']['verbose']:
            verbose_logger = getLogger('engine')
            verbose_logger.debug('local_func_scans: %s', local_func_scans)
    del local_func_scans

    return (wf, rpool, diff, blip, fmap_rp_list)


def ingress_output_dir(wf, cfg, rpool, unique_id, data_paths, part_id, ses_id, creds_path=None):

    dir_path = data_paths['derivatives_dir']

    print(f"\nPulling outputs from {dir_path}.\n")

    anat = os.path.join(dir_path, 'anat')
    func = os.path.join(dir_path, 'func')

    exts = ['.nii', '.gz', '.mat', '.1D', '.txt', '.csv', '.rms', '.tsv']

    outdir_anat = []
    outdir_func = []
    func_paths = {}
    func_dict = {}

    for subdir in [anat, func]:
        if os.path.isdir(subdir):
            for filename in os.listdir(subdir):
                for ext in exts:
                    if ext in filename:
                        if subdir == anat:
                            outdir_anat.append(os.path.join(subdir,
                                                    filename))
                        else:
                            outdir_func.append(os.path.join(subdir,
                                                    filename))

     # Add derivatives directory to rpool
    ingress = create_general_datasource(f'gather_derivatives_dir')
    ingress.inputs.inputnode.set(
            unique_id=unique_id,
            data=dir_path,
            creds_path=creds_path,
            dl_dir=cfg.pipeline_setup['working_directory']['path']
        )
    rpool.set_data("derivatives-dir", ingress, 'outputspec.data',
                {}, "", "outdir_config_ingress")

    for subdir in [outdir_anat, outdir_func]:
        for filepath in subdir:
            filename = str(filepath)
            for ext in exts:
                filename = filename.split("/")[-1].replace(ext, '')

            data_label = filename.split(unique_id)[1].lstrip('_')

            if len(filename) == len(data_label):
                raise Exception('\n\n[!] Possibly wrong participant or '
                                'session in this directory?\n\n'
                                f'Filepath: {filepath}\n\n')

            bidstag = ''
            for tag in data_label.split('_'):
                for prefix in ['task-', 'run-', 'acq-', 'rec']:
                    if tag.startswith(prefix):
                        bidstag += f'{tag}_'
                        data_label = data_label.replace(f'{tag}_', '')
            data_label, json = strip_template(data_label, dir_path, filename)

            rpool, json_info, pipe_idx, node_name, data_label = \
                json_outdir_ingress(rpool, filepath, \
                exts, data_label, json)

            if ('template' in data_label and not json_info['Template'] == \
                    cfg.pipeline_setup['outdir_ingress']['Template']):
                continue
            # Rename confounds to avoid confusion in nuisance regression
            if data_label.endswith('desc-confounds_timeseries'):
                data_label = 'pipeline-ingress_desc-confounds_timeseries'

            if len(bidstag) > 1:
                # Remove tail symbol
                bidstag = bidstag[:-1]
                if bidstag.startswith('task-'):
                    bidstag = bidstag.replace('task-', '')

            # Rename bold mask for CPAC naming convention
            # and to avoid collision with anat brain mask
            if data_label.endswith('desc-brain_mask') and filepath in outdir_func: 
                data_label = data_label.replace('brain_mask', 'bold_mask')

            try:
                pipe_x = rpool.get_pipe_number(pipe_idx)
            except ValueError:
                pipe_x = len(rpool.pipe_list)
            if filepath in outdir_anat:
                ingress = create_general_datasource(f'gather_anat_outdir_{str(data_label)}_{pipe_x}')
                ingress.inputs.inputnode.set(
                    unique_id=unique_id,
                    data=filepath,
                    creds_path=creds_path,
                    dl_dir=cfg.pipeline_setup['working_directory']['path']
                )
                rpool.set_data(data_label, ingress, 'outputspec.data', json_info,
                    pipe_idx, node_name, f"outdir_{data_label}_ingress", inject=True)
            else:
                if data_label.endswith('desc-preproc_bold'): 
                    func_key = data_label
                    func_dict[bidstag] = {}
                    func_dict[bidstag]['scan'] = str(filepath)
                    func_dict[bidstag]['scan_parameters'] = json_info
                    func_dict[bidstag]['pipe_idx'] = pipe_idx
                if data_label.endswith('desc-brain_mask'): 
                    data_label = data_label.replace('brain_mask', 'bold_mask')
                try:
                    func_paths[data_label].append(filepath)
                except:
                    func_paths[data_label] = []
                    func_paths[data_label].append(filepath)

    if func_dict:
        wf, rpool = func_outdir_ingress(wf, cfg, func_dict, rpool, unique_id, \
            creds_path, part_id, func_key, func_paths)

    if cfg.surface_analysis['freesurfer']['ingress_reconall']:
        rpool = ingress_freesurfer(wf, rpool, cfg, data_paths, unique_id, part_id,
                          ses_id)
    return wf, rpool

def json_outdir_ingress(rpool, filepath, exts, data_label, json):
    
    desc_val = None
    for tag in data_label.split('_'):
        if 'desc-' in tag:
            desc_val = tag
            break
    jsonpath = str(filepath)
    for ext in exts:
        jsonpath = jsonpath.replace(ext, '')
    jsonpath = f"{jsonpath}.json"

    if not os.path.exists(jsonpath):
        print(f'\n\n[!] No JSON found for file {filepath}.\nCreating '
            f'{jsonpath}..\n\n')
        json_info = {
            'Description': 'This data was generated elsewhere and '
                        'supplied by the user into this C-PAC run\'s '
                        'output directory. This JSON file was '
                        'automatically generated by C-PAC because a '
                        'JSON file was not supplied with the data.'
        }
        json_info = {**json_info, **json}
        write_output_json(json_info, jsonpath)
    else:
        json_info = read_json(jsonpath)
        json_info = {**json_info, **json}
    if 'CpacProvenance' in json_info:
        if desc_val:
            # it's a C-PAC output, let's check for pipe_idx/strat integer
            # suffixes in the desc- entries.
            only_desc = str(desc_val)
        
            if only_desc[-1].isdigit():
                for idx in range(0, 3):
                    # let's stop at 3, please don't run >999 strategies okay?
                    if only_desc[-1].isdigit():
                        only_desc = only_desc[:-1]
        
                if only_desc[-1] == '-':
                    only_desc = only_desc.rstrip('-')
                else:
                    raise Exception('\n[!] Something went wrong with either '
                                    'reading in the output directory or when '
                                    'it was written out previously.\n\nGive '
                                    'this to your friendly local C-PAC '
                                    f'developer:\n\n{str(data_label)}\n')

            # remove the integer at the end of the desc-* variant, we will 
            # get the unique pipe_idx from the CpacProvenance below
            data_label = data_label.replace(desc_val, only_desc)

        # preserve cpac provenance/pipe_idx
        pipe_idx = rpool.generate_prov_string(json_info['CpacProvenance'])
        node_name = ""
        
    else:
        json_info['CpacProvenance'] = [f'{data_label}:Non-C-PAC Origin: {filepath}']
        if not 'Description' in json_info:
            json_info['Description'] = 'This data was generated elsewhere and ' \
                                    'supplied by the user into this C-PAC run\'s '\
                                    'output directory. This JSON file was '\
                                    'automatically generated by C-PAC because a '\
                                    'JSON file was not supplied with the data.'
        pipe_idx = rpool.generate_prov_string(json_info['CpacProvenance'])
        node_name = f"{data_label}_ingress"

    return rpool, json_info, pipe_idx, node_name, data_label

def func_outdir_ingress(wf, cfg, func_dict, rpool, unique_id, creds_path, part_id, key, \
                            func_paths):
    pipe_x = len(rpool.pipe_list)
    exts = ['.nii', '.gz', '.mat', '.1D', '.txt', '.csv', '.rms', '.tsv']
    ingress = create_func_datasource(func_dict, rpool, f'gather_func_outdir_{key}_{pipe_x}')
    ingress.inputs.inputnode.set(
        subject=unique_id,
        creds_path=creds_path,
        dl_dir=cfg.pipeline_setup['working_directory']['path']
    )
    rpool.set_data('subject', ingress, 'outputspec.subject', {}, "",
        "func_ingress")
    ingress.get_node('inputnode').iterables = \
        ("scan", list(func_dict.keys())) 
    rpool.set_data(key, ingress, 'outputspec.rest', {}, "",
            "func_ingress")
    
    rpool.set_data('scan', ingress, 'outputspec.scan', {}, "", 'func_ingress')
    rpool.set_data('scan-params', ingress, 'outputspec.scan_params', {}, "",
        "scan_params_ingress")
    wf, rpool, diff, blip, fmap_rp_list = ingress_func_metadata(wf, cfg, \
            rpool, func_dict, part_id, creds_path, key)
    
    # Have to do it this weird way to save the parsed BIDS tag & filepath
    mask_paths_key = 'desc-bold_mask' if 'desc-bold_mask' in func_paths else \
                                    'space-template_desc-bold_mask'
    ts_paths_key = 'pipeline-ingress_desc-confounds_timeseries'

    # Connect func data with approproate scan name
    iterables = pe.Node(Function(input_names=['scan',
                                              'mask_paths',
                                              'ts_paths'],
                                output_names=['out_scan', 
                                              'mask',
                                              'confounds'],
                                function=set_iterables),
                                name=f'set_iterables_{pipe_x}')
    iterables.inputs.mask_paths = func_paths[mask_paths_key]
    iterables.inputs.ts_paths = func_paths[ts_paths_key]
    wf.connect(ingress, 'outputspec.scan', iterables, 'scan')

    for key in func_paths:
        if key == mask_paths_key or key == ts_paths_key:
            ingress_func = create_general_datasource(f'ingress_func_data_{key}')
            ingress_func.inputs.inputnode.set(
                unique_id=unique_id,
                creds_path=creds_path,
                dl_dir=cfg.pipeline_setup['working_directory']['path'])
            wf.connect(iterables, 'out_scan', ingress_func, 'inputnode.scan')
            if key == mask_paths_key:
                wf.connect(iterables, 'mask', ingress_func, 'inputnode.data')
                rpool.set_data(key, ingress_func, 'inputnode.data', {}, "", f"outdir_{key}_ingress")
            elif key == ts_paths_key:
                wf.connect(iterables, 'confounds', ingress_func, 'inputnode.data')
                rpool.set_data(key, ingress_func, 'inputnode.data', {}, "", f"outdir_{key}_ingress")

    return wf, rpool

def set_iterables(scan, mask_paths=None, ts_paths=None):
    
    # match scan with filepath to get filepath
    mask_path = [path for path in mask_paths if scan in path]
    ts_path = [path for path in ts_paths if scan in path]

    return (scan, mask_path[0], ts_path[0]) 

def strip_template(data_label, dir_path, filename):
    
    json = {}
    # rename to template 
    for prefix in ['space-', 'from-', 'to-']: 
        for bidstag in data_label.split('_'):
            if bidstag.startswith(prefix):
                template_key, template_val = bidstag.split('-')
                template_name, _template_desc = lookup_identifier(template_val)
                if template_name:
                    json['Template'] = template_val
                    data_label = data_label.replace(template_val, 'template')
            elif bidstag.startswith('res-'):
                res_key, res_val = bidstag.split('-')
                json['Resolution'] = res_val
                data_label = data_label.replace(bidstag, '')
    if data_label.find('__'): data_label = data_label.replace('__', '_') 
    return data_label, json


def ingress_pipeconfig_paths(cfg, rpool, unique_id, creds_path=None):
    # ingress config file paths
    # TODO: may want to change the resource keys for each to include one level up in the YAML as well

    import pkg_resources as p
    import pandas as pd
    import ast

    template_csv = p.resource_filename('CPAC', 'resources/cpac_templates.csv')
    template_df = pd.read_csv(template_csv, keep_default_na=False)
    
    for row in template_df.itertuples():
    
        key = row.Key
        val = row.Pipeline_Config_Entry
        val = cfg.get_nested(cfg, [x.lstrip() for x in val.split(',')])
        resolution = row.Intended_Resolution_Config_Entry
        desc = row.Description

        if not val:
            continue

        if resolution:
            res_keys = [x.lstrip() for x in resolution.split(',')]
            tag = res_keys[-1]
        json_info = {} 

        if '$FSLDIR' in val:
            val = val.replace('$FSLDIR', cfg.pipeline_setup[
                'system_config']['FSLDIR'])
        if '$priors_path' in val:
            priors_path = cfg.segmentation['tissue_segmentation']['FSL-FAST']['use_priors']['priors_path'] or ''
            if '$FSLDIR' in priors_path:
                priors_path = priors_path.replace('$FSLDIR', cfg.pipeline_setup['system_config']['FSLDIR'])
            val = val.replace('$priors_path', priors_path)
        if '${resolution_for_anat}' in val:
            val = val.replace('${resolution_for_anat}', cfg.registration_workflows['anatomical_registration']['resolution_for_anat'])               
        if '${func_resolution}' in val:
            val = val.replace('${func_resolution}', cfg.registration_workflows[
                'functional_registration']['func_registration_to_template'][
                'output_resolution'][tag])

        if desc:
            template_name, _template_desc = lookup_identifier(val)
            if template_name:
                desc = f"{template_name} - {desc}"
            json_info['Description'] = f"{desc} - {val}"
        if resolution:
            resolution = cfg.get_nested(cfg, res_keys)
            json_info['Resolution'] = resolution

            resampled_template = pe.Node(Function(input_names=['resolution',
                                                               'template',
                                                               'template_name',
                                                               'tag'],
                                                  output_names=['resampled_template'],
                                                  function=resolve_resolution,
                                                  as_module=True),
                                         name='resampled_' + key)

            resampled_template.inputs.resolution = resolution
            resampled_template.inputs.template = val
            resampled_template.inputs.template_name = key
            resampled_template.inputs.tag = tag
            
            # the set_data below is set up a little differently, because we are
            # injecting and also over-writing already-existing entries
            #   other alternative would have been to ingress into the
            #   resampled_template node from the already existing entries, but we
            #   didn't do that here
            rpool.set_data(key,
                           resampled_template,
                           'resampled_template',
                           json_info, "",
                           "template_resample") #, inject=True)   # pipe_idx (after the blank json {}) should be the previous strat that you want deleted! because you're not connecting this the regular way, you have to do it manually

        else:
            if val:
                config_ingress = create_general_datasource(f'gather_{key}')
                config_ingress.inputs.inputnode.set(
                    unique_id=unique_id,
                    data=val,
                    creds_path=creds_path,
                    dl_dir=cfg.pipeline_setup['working_directory']['path']
                )
                rpool.set_data(key, config_ingress, 'outputspec.data',
                               json_info, "", f"{key}_config_ingress")
    # templates, resampling from config
    '''
    template_keys = [
        ("anat", ["network_centrality", "template_specification_file"]),
        ("anat", ["nuisance_corrections", "2-nuisance_regression",
                  "lateral_ventricles_mask"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "FSL-FAST", "use_priors",
          "CSF_path"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "FSL-FAST", "use_priors",
          "GM_path"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "FSL-FAST", "use_priors",
          "WM_path"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "Template_Based", "CSF"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "Template_Based", "GRAY"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "Template_Based", "WHITE"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T1w_ACPC_template"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T1w_brain_ACPC_template"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T2w_ACPC_template"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T2w_brain_ACPC_template"])]

    def get_nested_attr(c, template_key):
        attr = getattr(c, template_key[0])
        keys = template_key[1:]

        def _get_nested(attr, keys):
            if len(keys) > 1:
                return (_get_nested(attr[keys[0]], keys[1:]))
            elif len(keys):
                return (attr[keys[0]])
            else:
                return (attr)

        return (_get_nested(attr, keys))

    def set_nested_attr(c, template_key, value):
        attr = getattr(c, template_key[0])
        keys = template_key[1:]

        def _set_nested(attr, keys):
            if len(keys) > 1:
                return (_set_nested(attr[keys[0]], keys[1:]))
            elif len(keys):
                attr[keys[0]] = value
            else:
                return (attr)

        return (_set_nested(attr, keys))

    for key_type, key in template_keys:
        attr = cfg.get_nested(cfg, key)
        if isinstance(attr, str) or attr == None:
            node = create_check_for_s3_node(
                key[-1],
                attr, key_type,
                data_paths['creds_path'],
                cfg.pipeline_setup['working_directory']['path'],
                map_node=False
            )
            cfg.set_nested(cfg, key, node)

    template_keys_in_list = [
        ("anat",
         ["segmentation", "tissue_segmentation", "ANTs_Prior_Based",
          "template_brain_list"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "ANTs_Prior_Based",
          "template_segmentation_list"]),
    ]

    for key_type, key in template_keys_in_list:
        node = create_check_for_s3_node(
            key[-1],
            cfg.get_nested(cfg, key), key_type,
            data_paths['creds_path'],
            cfg.pipeline_setup['working_directory']['path'],
            map_node=True
        )
        cfg.set_nested(cfg, key, node)
    '''

    return rpool


def initiate_rpool(wf, cfg, data_paths=None, part_id=None):
    '''

    data_paths format:
      {'anat': {
            'T1w': '{T1w path}',
            'T2w': '{T2w path}'
        },
       'creds_path': {None OR path to credentials CSV},
       'func': {
           '{scan ID}':
               {
                   'scan': '{path to BOLD}',
                   'scan_parameters': {scan parameter dictionary}
               }
       },
       'site_id': 'site-ID',
       'subject_id': 'sub-01',
       'unique_id': 'ses-1',
       'derivatives_dir': '{derivatives_dir path}'}
    '''

    # TODO: refactor further, integrate with the ingress_data functionality
    # TODO: used for BIDS-Derivatives (below), and possible refactoring of
    # TODO: the raw data config to use 'T1w' label instead of 'anat' etc.

    if data_paths:
        part_id = data_paths['subject_id']
        ses_id = data_paths['unique_id']
        if 'creds_path' not in data_paths:
            creds_path = None
        else:
            creds_path = data_paths['creds_path']
        unique_id = f'{part_id}_{ses_id}'
    
    elif part_id:
        unique_id = part_id
        creds_path = None

    rpool = ResourcePool(name=unique_id, cfg=cfg)

    if data_paths:
        # ingress outdir
        try: 
            if data_paths['derivatives_dir'] and cfg.pipeline_setup['outdir_ingress']['run']:
                wf, rpool = \
                     ingress_output_dir(wf, cfg, rpool, unique_id, data_paths, part_id, \
                    ses_id, creds_path=None)
        except:
            rpool = ingress_raw_anat_data(wf, rpool, cfg, data_paths, unique_id,
                                        part_id, ses_id)
            if 'func' in data_paths:
                wf, rpool, diff, blip, fmap_rp_list = \
                    ingress_raw_func_data(wf, rpool, cfg, data_paths, unique_id,
                                        part_id, ses_id)

    # grab any file paths from the pipeline config YAML
    rpool = ingress_pipeconfig_paths(cfg, rpool, unique_id, creds_path)

    # output files with 4 different scans

    return (wf, rpool)


def run_node_blocks(blocks, data_paths, cfg=None):
    import os
    from CPAC.pipeline import nipype_pipeline_engine as pe
    from CPAC.pipeline.engine import NodeBlock

    if not cfg:
        cfg = {
            'pipeline_setup': {
                'working_directory': {
                    'path': os.getcwd()
                },
                'log_directory': {
                    'path': os.getcwd()
                }
            }
        }

    # TODO: WE HAVE TO PARSE OVER UNIQUE ID'S!!!
    _, rpool = initiate_rpool(cfg, data_paths)

    wf = pe.Workflow(name='node_blocks')
    wf.base_dir = cfg.pipeline_setup['working_directory']['path']
    wf.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': cfg.pipeline_setup['log_directory']['path']
    }

    run_blocks = []
    if rpool.check_rpool('desc-preproc_T1w'):
        print("Preprocessed T1w found, skipping anatomical preprocessing.")
    else:
        run_blocks += blocks[0]
    if rpool.check_rpool('desc-preproc_bold'):
        print("Preprocessed BOLD found, skipping functional preprocessing.")
    else:
        run_blocks += blocks[1]

    for block in run_blocks:
        wf = NodeBlock(block, debug=cfg['pipeline_setup', 'Debugging',
                                        'verbose']).connect_block(
                                            wf, cfg, rpool)
    rpool.gather_pipes(wf, cfg)

    wf.run()


class NodeData:
    r"""Class to hold outputs of
    CPAC.pipeline.engine.ResourcePool().get_data(), so one can do

    ``node_data = strat_pool.node_data(resource)`` and have
    ``node_data.node`` and ``node_data.out`` instead of doing
    ``node, out = strat_pool.get_data(resource)`` and needing two
    variables (``node`` and ``out``) to store that information.

    Also includes ``variant`` attribute providing the resource's self-
    keyed value within its ``CpacVariant`` dictionary.

    Examples
    --------
    >>> rp = ResourcePool()
    >>> rp.node_data(None)
    NotImplemented (NotImplemented)

    >>> rp.set_data('test',
    ...             pe.Node(Function(input_names=[]), 'test'),
    ...             'b', [], 0, 'test')
    >>> rp.node_data('test')
    test (b)
    >>> rp.node_data('test').out
    'b'

    >>> try:
    ...     rp.node_data('b')
    ... except LookupError as lookup_error:
    ...     print(str(lookup_error).strip().split('\n')[0].strip())
    [!] C-PAC says: None of the listed resources are in the resource pool:
    """
    # pylint: disable=too-few-public-methods
    def __init__(self, strat_pool=None, resource=None, **kwargs):
        self.node = NotImplemented
        self.out = NotImplemented
        if strat_pool is not None and resource is not None:
            self.node, self.out = strat_pool.get_data(resource, **kwargs)

    def __repr__(self):
        return f'{getattr(self.node, "name", str(self.node))} ({self.out})'
