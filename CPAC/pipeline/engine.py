import os
import ast
import six
import json
import warnings
import logging
import copy
from unittest import TestCase

from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.utility import Rename
from CPAC.utils.interfaces.function import Function
from CPAC.utils.interfaces.datasink import DataSink

from CPAC.registration.registration import transform_derivative
from CPAC.nuisance import NuisanceRegressor

from CPAC.utils.utils import read_json, create_id_string, write_output_json, \
    get_last_prov_entry, ordereddict_to_dict, check_prov_for_regtool
from CPAC.utils.datasource import (
    create_anat_datasource,
    create_func_datasource,
    ingress_func_metadata,
    create_general_datasource,
    create_check_for_s3_node,
    resolve_resolution
)
from CPAC.image_utils.spatial_smoothing import spatial_smoothing
from CPAC.image_utils.statistical_transforms import z_score_standardize, \
    fisher_z_score_standardize

logger = logging.getLogger('workflow')


class ResourcePool(object):
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
            self.run_zscoring = 'z-scored' in cfg.post_processing[
                'z-scoring']['output']
            self.fwhm = cfg.post_processing['spatial_smoothing']['fwhm']
            self.smooth_opts = cfg.post_processing['spatial_smoothing'][
                'smoothing_method']

        self.xfm = ['alff', 'falff', 'reho']

        self.smooth = ['alff', 'falff', 'reho', 'vmhc',
                       'space-template_alff',
                       'space-template_falff',
                       'space-template_reho',
                       'degree-centrality',
                       'eigen-centrality',
                       'lfcd']
        self.zscore = self.smooth + ['alff', 'falff', 'reho',
                                     'desc-sm_alff',
                                     'desc-sm_falff',
                                     'desc-sm_reho',
                                     'space-template_alff',
                                     'space-template_falff',
                                     'space-template_reho',
                                     'space-template_degree-centrality',
                                     'space-template_eigen-centrality',
                                     'space-template_lfcd'
                                     'space-template_desc-sm_alff',
                                     'space-template_desc-sm_falff',
                                     'space-template_desc-sm_reho',
                                     'space-template_desc-sm_degree-centrality',
                                     'space-template_desc-sm_eigen-centrality',
                                     'space-template_desc-sm_lfcd']
        self.fisher_zscore = ['desc-MeanSCA_correlations']

    def append_name(self, name):
        self.name.append(name)

    def get_name(self):
        return self.name

    def check_rpool(self, resource):
        if resource not in self.rpool:
            return False
        return True

    def get_pipe_number(self, pipe_idx):
        return self.pipe_list.index(pipe_idx)

    def get_entire_rpool(self):
        return self.rpool

    def get_resources(self):
        return self.rpool.keys()

    def copy_rpool(self):
        return ResourcePool(rpool=copy.deepcopy(self.get_entire_rpool()),
                            name=self.name,
                            cfg=self.cfg,
                            pipe_list=copy.deepcopy(self.pipe_list))

    def get_raw_label(self, resource):
        # remove desc-* label
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
        return self.rpool[resource][pipe_idx][key]

    def get_resource_from_prov(self, prov):
        # each resource (i.e. "desc-cleaned_bold" AKA nuisance-regressed BOLD
        # data) has its own provenance list. the name of the resource, and
        # the node that produced it, is always the last item in the provenance
        # list, with the two separated by a colon :
        if isinstance(prov[-1], list):
            return prov[-1][-1].split(':')[0]
        elif isinstance(prov[-1], str):
            return prov[-1].split(':')[0]

    def set_data(self, resource, node, output, json_info, pipe_idx, node_name,
                 fork=False, inject=False):
        '''
        pipe_idx, node_name = new_id
        if f';{resource}:' not in pipe_idx:
            pipe_idx = f'{pipe_idx};{resource}:'   # <--- doing this up here, now, because the del self.rpool[resource][pipe_idx] below should only get deleted for the same input/output tag!
        if resource not in self.rpool.keys():
            self.rpool[resource] = {}
        else:
            if not fork:     # <--- in the event of multiple strategies/options, this will run for every option; just keep in mind
                if pipe_idx in self.rpool[resource].keys():  # <--- in case the resource name is now new, and not the original
                    del self.rpool[resource][pipe_idx]  # <--- remove old keys so we don't end up with a new strat for every new node unit (unless we fork)
        if pipe_idx[-1] == ';' or pipe_idx[-1] == ':':  # <--- if the ':', this kicks off when the pipe_idx is only something like 'T1w:', at the beginning
            new_name = node_name                        #      but how do we manage new threads, mid-pipeline?
        else:
            new_name = f',{node_name}'
        new_pipe_idx = f'{pipe_idx}{new_name}'
        '''
        json_info = json_info.copy()

        cpac_prov = []
        if 'CpacProvenance' in json_info:
            cpac_prov = json_info['CpacProvenance']
        new_prov_list = list(cpac_prov)   # <---- making a copy, it was already a list
        if not inject:
            new_prov_list.append(f'{resource}:{node_name}')
        res, new_pipe_idx = self.generate_prov_string(new_prov_list)

        if not json_info:
            json_info = {'RawSources': [resource]}     # <---- this will be repopulated to the full file path at the end of the pipeline building, in gather_pipes()
        json_info['CpacProvenance'] = new_prov_list

        if resource not in self.rpool.keys():
            self.rpool[resource] = {}
        else:
            if not fork:     # <--- in the event of multiple strategies/options, this will run for every option; just keep in mind
                if pipe_idx in self.rpool[resource].keys():  # <--- in case the resource name is now new, and not the original
                    del self.rpool[resource][pipe_idx]  # <--- remove old keys so we don't end up with a new strat for every new node unit (unless we fork)

        if new_pipe_idx not in self.rpool[resource]:
            self.rpool[resource][new_pipe_idx] = {}
        if new_pipe_idx not in self.pipe_list:
            self.pipe_list.append(new_pipe_idx)
        self.rpool[resource][new_pipe_idx]['data'] = (node, output)
        self.rpool[resource][new_pipe_idx]['json'] = json_info

    def get(self, resource, pipe_idx=None, report_fetched=False,
            optional=False):
        # NOTE!!!
        #   if this is the main rpool, this will return a dictionary of strats, and inside those, are dictionaries like {'data': (node, out), 'json': info}
        #   BUT, if this is a sub rpool (i.e. a strat_pool), this will return a one-level dictionary of {'data': (node, out), 'json': info} WITHOUT THE LEVEL OF STRAT KEYS ABOVE IT
        if isinstance(resource, list):
            # if a list of potential inputs are given, pick the first one
            # found
            for label in resource:
                if label in self.rpool.keys():
                    if report_fetched:
                        return (self.rpool[label], label)
                    return self.rpool[label]
            else:
                if optional:
                    if report_fetched:
                        return (None, None)
                    return None
                raise Exception("\n[!] C-PAC says: None of the listed "
                                "resources are in the resource pool:\n"
                                f"{resource}\n")
        else:
            if resource not in self.rpool.keys():
                if optional:
                    if report_fetched:
                        return (None, None)
                    return None
                raise LookupError("\n\n[!] C-PAC says: The listed resource is "
                                  f"not in the resource pool:\n{resource}\n\n"
                                  "Developer Note: This may be due to a mis"
                                  "match between the node block's docstring "
                                  "'input' field and a strat_pool.get_data() "
                                  "call within the block function.\n")
            if report_fetched:
                if pipe_idx:
                    return (self.rpool[resource][pipe_idx], resource)
                return (self.rpool[resource], resource)
            if pipe_idx:
                return self.rpool[resource][pipe_idx]
            return self.rpool[resource]

    def get_data(self, resource, pipe_idx=None, report_fetched=False,
                 quick_single=False):
        if quick_single:
            for key, val in self.get(resource).items():
                return val['data']
        if report_fetched:
            if pipe_idx:
                connect, fetched = self.get(resource, pipe_idx=pipe_idx,
                                            report_fetched=report_fetched)
                return (connect['data'], fetched)
            connect, fetched =self.get(resource,
                                       report_fetched=report_fetched)
            return (connect['data'], fetched)
        if pipe_idx:
            return self.get(resource, pipe_idx=pipe_idx)['data']
        return self.get(resource)['data']

    def copy_resource(self, resource, new_name):
        self.rpool[new_name] = self.rpool[resource]

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
                            f'is  incomplete.\n')
        return strat_json

    def get_cpac_provenance(self, resource, strat=None):
        # NOTE: resource_strat_dct has to be entered properly by the developer
        # it has to either be rpool[resource][strat] or strat_pool[resource]
        json_data = self.get_json(resource, strat)
        return json_data['CpacProvenance']

    def generate_prov_string(self, prov):
        # this will generate a string from a SINGLE RESOURCE'S dictionary of
        # MULTIPLE PRECEDING RESOURCES (or single, if just one)
        #   NOTE: this DOES NOT merge multiple resources!!! (i.e. for merging-strat pipe_idx generation)
        if not isinstance(prov, list):
            raise Exception('\n[!] Developer info: the CpacProvenance '
                            f'entry for {prov} has to be a list.\n')
        last_entry = get_last_prov_entry(prov)
        resource = last_entry.split(':')[0]
        return (resource, str(prov))

    def get_resource_strats_from_prov(self, prov):
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

    def get_strats(self, resources):

        # TODO: NOTE: NOT COMPATIBLE WITH SUB-RPOOL/STRAT_POOLS
        # TODO: (and it doesn't have to be)

        import itertools
        
        linked_resources = []
        resource_list = []
        for resource in resources:
            # grab the linked-input tuples
            if isinstance(resource, tuple):
                linked = []
                for label in list(resource):
                    #if isinstance(label, list):
                    rp_dct, fetched_resource = self.get(label,
                                                        report_fetched=True,
                                                        optional=True)
                    if not rp_dct:
                        continue
                    linked.append(fetched_resource)
                    #else:
                    #    linked.append(label)
                linked_resources.append(linked)
                resource_list += linked
            else:
                resource_list.append(resource)

        total_pool = []
        len_inputs = len(resource_list)
        for resource in resource_list:
            rp_dct, fetched_resource = self.get(resource,
                                                report_fetched=True,             # <---- rp_dct has the strats/pipe_idxs as the keys on first level, then 'data' and 'json' on each strat level underneath
                                                optional=True)                   # oh, and we make the resource fetching in get_strats optional so we can have optional inputs, but they won't be optional in the node block unless we want them to be
            if not rp_dct:
                len_inputs -= 1
                continue
            sub_pool = []
            for strat in rp_dct.keys():
                json_info = self.get_json(fetched_resource, strat)
                cpac_prov = json_info['CpacProvenance']
                sub_pool.append(cpac_prov)
            total_pool.append(sub_pool)

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
            for strat_tuple in strats:
                # the strats are still in provenance-list form
                #   not string-based pipe_idx's yet
                strat_list = list(strat_tuple)     # <------- strat_list is now a list of strats all combined together, one of the permutations. keep in mind each strat in the combo comes from a different data source/input
                                                   #          like this:   strat_list = [desc-preproc_T1w:pipe_idx_anat_1, desc-brain_mask:pipe_idx_mask_1]

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
                        for xlabel in linked:
                            if drop:
                                break
                            xjson = json_dct[xlabel]
                            for ylabel in linked:
                                if xlabel == ylabel:
                                    continue
                                yjson = json_dct[ylabel]

                                if 'CpacVariant' not in xjson and 'CpacVariant' in yjson:
                                    raw_xlabel = self.get_raw_label(xlabel)
                                    if raw_xlabel in yjson['CpacVariant']:
                                        drop = True
                                        #print(f'xlabel: {xlabel}')
                                        #print(f'ylabel: {ylabel}')
                                        #print('yjson had cpacvariant and xjson didnt, but something in xjson exists in the yjson variant list')
                                        break
                                    for entry in self.get_resource_strats_from_prov(xjson['CpacProvenance']):
                                        raw_label = self.get_raw_label(entry)
                                        if raw_label in yjson['CpacVariant']:
                                            drop = True
                                            break
                                        else:
                                            for sub_entry in self.get_resource_strats_from_prov(entry):
                                                raw_sub_entry = self.get_raw_label(sub_entry)
                                                if raw_sub_entry in yjson['CpacVariant']:
                                                    drop = True
                                                    break

                                elif 'CpacVariant' not in yjson and 'CpacVariant' in xjson:
                                    raw_ylabel = self.get_raw_label(ylabel)
                                    if raw_ylabel in xjson['CpacVariant']:
                                        drop = True
                                        #print(f'xlabel: {xlabel}')
                                        #print(f'ylabel: {ylabel}')
                                        #print('xjson had cpacvariant and yjson didnt, but something in yjson exists in the xjson variant list')
                                        break
                                    for entry in self.get_resource_strats_from_prov(yjson['CpacProvenance']):
                                        raw_label = self.get_raw_label(entry)
                                        if raw_label in xjson['CpacVariant']:
                                            drop = True
                                            break
                                        else:
                                            for sub_entry in self.get_resource_strats_from_prov(entry):
                                                raw_sub_entry = self.get_raw_label(sub_entry)
                                                if raw_sub_entry in xjson['CpacVariant']:
                                                    drop = True
                                                    break

                                elif 'CpacVariant' not in xjson and 'CpacVariant' not in yjson:
                                    continue

                                else:
                                    for younger_resource in xjson['CpacVariant']:
                                        if younger_resource in yjson['CpacVariant']:
                                            if xjson['CpacVariant'][younger_resource] != yjson['CpacVariant'][younger_resource]:
                                                drop = True
                                                break
                                    else:
                                        # if no drops - but we still need to
                                        # look for absences of CpacVariant
                                        # younger resources
                                        for younger_resource, variant_list in xjson['CpacVariant'].items():
                                            yprov = copy.deepcopy(yjson['CpacProvenance'])
                                            flat_yprov = self.flatten_prov(yprov)
                                            raw_labels = []
                                            flat_raw_yprov = []
                                            for entry in flat_yprov:
                                                raw_resource = self.get_raw_label(entry.split(':')[0])
                                                entry = entry.replace(entry.split(':')[0],
                                                                      raw_resource)
                                                raw_labels.append(raw_resource)
                                                flat_raw_yprov.append(entry)
                                            if younger_resource in raw_labels:
                                                for variant in variant_list:
                                                    if f'{younger_resource}:{variant}' not in flat_raw_yprov:
                                                        drop = True
                                                        break

                                    for younger_resource in yjson['CpacVariant']:
                                        if younger_resource in xjson['CpacVariant']:
                                            if xjson['CpacVariant'][younger_resource] != yjson['CpacVariant'][younger_resource]:
                                                drop = True
                                                break
                                    else:
                                        # if no drops - but we still need to
                                        # look for absences of CpacVariant
                                        # younger resources
                                        for younger_resource, variant_list in yjson['CpacVariant'].items():
                                            xprov = copy.deepcopy(xjson['CpacProvenance'])
                                            flat_xprov = self.flatten_prov(xprov)
                                            raw_labels = []
                                            flat_raw_xprov = []
                                            for entry in flat_xprov:
                                                raw_resource = self.get_raw_label(entry.split(':')[0])
                                                entry = entry.replace(entry.split(':')[0],
                                                                      raw_resource)
                                                raw_labels.append(raw_resource)
                                                flat_raw_xprov.append(entry)
                                            if younger_resource in raw_labels:
                                                for variant in variant_list:
                                                    if f'{younger_resource}:{variant}' not in flat_raw_xprov:
                                                        drop = True
                                                        break
                if drop:
                    continue

                # make the merged strat label from the multiple inputs
                # strat_list is actually the merged CpacProvenance lists
                pipe_idx = str(strat_list)
                new_strats[pipe_idx] = ResourcePool()     # <----- new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                new_strats[pipe_idx].rpool['json'] = {}
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
        else:
            new_strats = {}
            for resource_strat_list in total_pool:       # total_pool will have only one list of strats, for the one input
                for cpac_prov in resource_strat_list:     # <------- cpac_prov here doesn't need to be modified, because it's not merging with other inputs
                    resource, pipe_idx = self.generate_prov_string(cpac_prov)
                    resource_strat_dct = self.rpool[resource][pipe_idx]   # <----- remember, this is the dct of 'data' and 'json'.
                    new_strats[pipe_idx] = ResourcePool(rpool={resource: resource_strat_dct})   # <----- again, new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                    new_strats[pipe_idx].rpool['json'] = resource_strat_dct['json']  # TODO: WARNING- THIS IS A LEVEL HIGHER THAN THE ORIGINAL 'JSON' FOR EASE OF ACCESS IN CONNECT_BLOCK WITH THE .GET(JSON)
                    new_strats[pipe_idx].rpool['json']['CpacProvenance'] = cpac_prov

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

                node, out = self.get_data("T1w_brain_template_deriv",
                                          quick_single=True)
                wf.connect(node, out, xfm, 'inputspec.reference')

                node, out = self.get_data('from-bold_to-template_mode-image_xfm',
                                          pipe_idx=xfm_idx)
                wf.connect(node, out, xfm, 'inputspec.transform')

                label = f'space-template_{label}'

                new_prov = json_info['CpacProvenance'] + xfm_prov
                json_info['CpacProvenance'] = new_prov
                new_pipe_idx = self.generate_prov_string(new_prov)

                self.set_data(label, xfm, 'outputspec.out_file', json_info,
                              new_pipe_idx, f'{label}_xfm_{num}', fork=True)

        return wf

    def post_process(self, wf, label, connection, json_info, pipe_idx, pipe_x,
                     outs):

        input_type = 'func_derivative'
        if 'centrality' in label or 'lfcd' in label:
            input_type = 'func_derivative_multi'

        if 'centrality' in label or 'lfcd' in label:
            mask = 'template_specification_file'
        elif 'space-template' in label:
            mask = 'space-template_res-derivative_desc-bold_mask'
        else:
            mask = 'space-bold_desc-brain_mask'

        mask_idx = None
        for entry in json_info['CpacProvenance']:
            if isinstance(entry, list):
                if entry[-1].split(':')[0] == mask:
                    mask_prov = entry
                    mask_idx = self.generate_prov_string(mask_prov)[1]
                    break

        if self.run_smoothing:
            if label in self.smooth:
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

                    self.set_data(smlabel, sm, 'outputspec.out_file',
                                  json_info, pipe_idx,
                                  f'spatial_smoothing_{smooth_opt}',
                                  fork=True)
                    self.set_data('fwhm', sm, 'outputspec.fwhm', json_info,
                                  pipe_idx, f'spatial_smoothing_{smooth_opt}',
                                  fork=True)

        if self.run_zscoring:

            if 'desc-' not in label:
                if 'space-template' in label:
                    label = label.replace('space-template',
                                          'space-template_desc-zstd')
                else:
                    label = f'desc-zstd_{label}'
            else:
                for tag in label.split('_'):
                    if 'desc-' in tag:
                        newtag = f'{tag}-zstd'
                        new_label = label.replace(tag, newtag)
                        break

            if label in self.zscore:

                zstd = z_score_standardize(f'{label}_zstd_{pipe_x}',
                                           input_type)

                wf.connect(connection[0], connection[1],
                           zstd, 'inputspec.in_file')

                node, out = self.get_data(mask, pipe_idx=mask_idx)
                wf.connect(node, out, zstd, 'inputspec.mask')

                self.set_data(new_label, zstd, 'outputspec.out_file',
                              json_info, pipe_idx, f'zscore_standardize',
                              fork=True)

            elif label in self.fisher_zscore:

                zstd = fisher_z_score_standardize(f'{label}_zstd_{pipe_x}',
                                                  label, input_type)

                wf.connect(connection[0], connection[1],
                           zstd, 'inputspec.correlation_file')

                # if the output is 'desc-MeanSCA_correlations', we want
                # 'desc-MeanSCA_timeseries'
                oned = label.replace('correlations', 'timeseries')

                node, out = outs[oned]
                wf.connect(node, out, zstd, 'inputspec.timeseries_oned')

                self.set_data(new_label, zstd, 'outputspec.out_file',
                              json_info, pipe_idx,
                              'fisher_zscore_standardize',
                              fork=True)

        return wf

    def gather_pipes(self, wf, cfg, all=False, add_incl=None, add_excl=None):
        # TODO: cpac_outputs.csv etc
        # TODO: might be better to do an inclusion instead
        non_sink = ['scan', 'TR', 'tpattern', 'start_tr', 'stop_tr',
                    'pe_direction', 'subject', 'atlas_name', 'scan_params',
                    'deltaTE', 'diff_phase_dwell', 'dwell_asym_ratio',
                    'diffphase_scan_params', 'diffmag_scan_params',
                    'motion-basefile']
        excl = ['T1w', 'bold', 'motion_basefile', 'bold_coreg_input',
                'diffphase', 'diffmag', 'epi']
        config_paths = ['T1w_ACPC_template', 'T1w_brain_ACPC_template',
                        'unet_model', 'T1w_brain_template', 'T1w_template',
                        'T1w_brain_template_mask',
                        'T1w_brain_template_symmetric',
                        'T1w_template_symmetric',
                        'dilated_symmetric_brain_mask',
                        'dilated_symmetric_brain_mask_for_template',
                        'T1w_brain_template_symmetric_for_resample',
                        'T1w_template_symmetric_for_resample', 'ref_mask',
                        'template_for_resample',
                        'T1w_brain_template_for_func',
                        'T1w_template_for_func',
                        'template_epi', 'template_epi_mask',
                        'lateral_ventricles_mask', 'eye_mask_path',
                        'EPI_template', 'EPI_template_mask',
                        'EPI_template_deriv', 'EPI_template_for_resample',
                        'EPI_template_funcreg', 'T1w_brain_template_deriv',
                        'T1w_template_deriv', 'T1w_brain_template_funcreg',
                        'T1w_template_funcreg',
                        'T1w_template_symmetric_funcreg',
                        'T1w_brain_template_symmetric_funcreg',
                        'T1w_brain_template_for_resample',
                        'T1w_template_for_resample']
        excl += non_sink
        excl += config_paths

        if add_excl:
            excl += add_excl

        anat = ['T1w', 'probseg', 'T1w-template']
        func = ['bold', 'timeseries', 'alff', 'falff', 'reho', 'vmhc',
                'correlations', 'statmap', 'regressors', 'degree-centrality',
                'eigen-centrality', 'lfcd']

        motions = ['motion', 'movement', 'coordinate', 'displacement',
                   'dvars', 'power-params']

        qc_anat = ['T1w-axial-qc',
                   'T1w-sagittal-qc',
                   'dseg-axial-qc',
                   'dseg-sagittal-qc']
        anat += qc_anat

        qc_func = ['bold-axial-qc',
                   'bold-sagittal-qc',
                   'bold-carpet-qc',
                   'framewise-displacement-jenkinson-plot-qc',
                   'movement-parameters-trans-qc',
                   'movement-parameters-rot-qc',
                   'bold-snr-axial-qc',
                   'bold-snr-sagittal-qc',
                   'bold-snr-hist-qc',
                   'bold-snr-qc']
        func += qc_func

        if all:
            excl = non_sink

        for resource in self.rpool.keys():
            # TODO: cpac_outputs.csv etc
            if resource in excl:
                continue
            subdir = 'other'
            if resource.split('_')[-1] in anat:
                subdir = 'anat'
                #TODO: get acq- etc.
            elif resource.split('_')[-1] in func:
                subdir = 'func'
                #TODO: other stuff like acq- etc.
            elif resource.split('_')[-1] == 'mask':
                if 'space-T1w' in resource:
                    subdir = 'anat'
                if 'label-CSF' in resource or 'label-GM' in resource or \
                        'label-WM' in resource:
                    subdir = 'anat'
                if 'space-bold' in resource:
                    subdir = 'func'
            elif resource.split('_')[-1] == 'xfm':
                if 'from-T1w' in resource:
                    subdir = 'anat'
                if 'from-bold' in resource:
                    subdir = 'func'
            else:
                for tag in motions:
                    if tag in resource:
                        subdir = 'func'

            for pipe_idx in self.rpool[resource]:
                unique_id = self.get_name()

                out_dir = cfg.pipeline_setup['output_directory']['path']
                pipe_name = cfg.pipeline_setup['pipeline_name']
                container = os.path.join(f'cpac_{pipe_name}', unique_id)
                filename = f'{unique_id}_{resource}'

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
            # TODO: cpac_outputs.csv etc
            if resource in excl:
                continue
                
            if not all:
                if 'symtemplate' in resource:
                    continue

            num_variant = 0
            if len(self.rpool[resource]) == 1:
                num_variant = ""
            for pipe_idx in self.rpool[resource]:

                pipe_x = self.get_pipe_number(pipe_idx)

                try:
                    num_variant += 1
                except TypeError:
                    pass

                json_info = self.rpool[resource][pipe_idx]['json']
                out_dct = self.rpool[resource][pipe_idx]['out']

                if out_dct['subdir'] == 'other' and not all:
                    continue

                unique_id = out_dct['unique_id']

                if num_variant:
                    for key in out_dct['filename'].split('_'):
                        if 'desc-' in key:
                            out_dct['filename'] = out_dct['filename'
                            ].replace(key, f'{key}-{num_variant}')
                            resource_idx = resource.replace(key,
                                                            f'{key}-{num_variant}')
                            break
                        else:
                            suff = resource.split('_')[-1]
                            newdesc_suff = f'desc-{num_variant}_{suff}'
                            resource_idx = resource.replace(suff,
                                                            newdesc_suff)
                else:
                    resource_idx = resource

                id_string = pe.Node(Function(input_names=['unique_id',
                                                          'resource',
                                                          'scan_id',
                                                          'atlas_id',
                                                          'fwhm'],
                                             output_names=['out_filename'],
                                             function=create_id_string),
                                    name=f'id_string_{resource_idx}_{pipe_x}')
                id_string.inputs.unique_id = unique_id
                id_string.inputs.resource = resource_idx

                # grab the iterable scan ID
                if out_dct['subdir'] == 'func':
                    node, out = self.rpool['scan']["['scan:func_ingress']"][
                        'data']
                    wf.connect(node, out, id_string, 'scan_id')
                    
                # grab the FWHM if smoothed
                for tag in resource.split('_'):
                    if 'desc-' in tag and '-sm' in tag:
                        fwhm_idx = pipe_idx.replace(f'{resource}:', 'fwhm:')
                        node, out = self.rpool['fwhm'][fwhm_idx]['data']
                        wf.connect(node, out, id_string, 'fwhm')
                        break

                atlas_suffixes = ['timeseries', 'correlations', 'statmap']
                # grab the iterable atlas ID
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
                wf.connect(id_string, 'out_filename',
                           nii_name, 'format_string')

                node, out = self.rpool[resource][pipe_idx]['data']
                wf.connect(node, out, nii_name, 'in_file')

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

                wf.connect(nii_name, 'out_file',
                           ds, f'{out_dct["subdir"]}.@data')
                wf.connect(write_json, 'json_file',
                           ds, f'{out_dct["subdir"]}.@json')


class NodeBlock(object):
    def __init__(self, node_block_functions):

        if not isinstance(node_block_functions, list):
            node_block_functions = [node_block_functions]

        self.node_blocks = {}

        for node_block_function in node_block_functions:    # <---- sets up the NodeBlock object in case you gave it a list of node blocks instead of a single one - for option forking.
            init_dct = self.grab_docstring_dct(node_block_function.__doc__)
            name = init_dct['name']
            self.name = name
            self.node_blocks[name] = {}

            for key, val in init_dct.items():
                self.node_blocks[name][key] = val

            self.node_blocks[name]['block_function'] = node_block_function

            #TODO: fix/replace below
            self.outputs = {}
            for out in init_dct['outputs']:
                self.outputs[out] = None

            self.options = ['base']
            if 'options' in init_dct:
                self.options = init_dct['options']

    def get_name(self):
        return self.name

    def grab_docstring_dct(self, fn_docstring):
        init_dct_schema = ['name', 'config', 'switch', 'option_key',
                           'option_val', 'inputs', 'outputs']
        if 'Node Block:' in fn_docstring:
            fn_docstring = fn_docstring.split('Node Block:')[1]
        fn_docstring = fn_docstring.replace('\n', '').replace(' ', '')
        #dct = json.loads(fn_docstring.replace('\n', '').replace(' ', ''))
        dct = ast.literal_eval(fn_docstring)
        #try:
        #    dct = json.loads(fn_docstring.replace('\n', '').replace(' ', ''))
        #except Exception as e:
        #    raise Exception('\n\n[!] Node block docstring error.\n\n'
        #                    f'Docstring:\n{fn_docstring}\n\n')
        for key in init_dct_schema:
            if key not in dct.keys():
                raise Exception('\n[!] Developer info: At least one of the '
                                'required docstring keys in your node block '
                                'is missing.\n\nNode block docstring keys:\n'
                                f'{init_dct_schema}\n\nYou provided:\n'
                                f'{dct.keys()}\n\nDocstring:\n{fn_docstring}'
                                '\n\n')
        return dct

    def check_null(self, val):
        if isinstance(val, str):
            val = None if val.lower() == 'none' else val
        return val

    def check_output(self, outputs, label, name):
        if label not in outputs:
            raise Exception('\n[!] Output name in the block function does '
                            'not match the outputs list in Node Block '
                            f'{name}\n')

    def grab_tiered_dct(self, cfg, key_list):
        cfg_dct = cfg
        for key in key_list:
            cfg_dct = cfg_dct.__getitem__(key)
        return cfg_dct

    def connect_block(self, wf, cfg, rpool):
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
                        if option in self.grab_tiered_dct(cfg, key_list):   # <---- goes over the option_vals in the node block docstring, and checks if the user's pipeline config included it in the forking list
                            opts.append(option)
            else:                                                           #         AND, if there are multiple option-val's (in a list) in the docstring, it gets iterated below in 'for opt in option' etc. AND THAT'S WHEN YOU HAVE TO DELINEATE WITHIN THE NODE BLOCK CODE!!!
                opts = [None]
            all_opts += opts

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
                else:
                    key_list = switch
                switch = self.grab_tiered_dct(cfg, key_list)
                if not isinstance(switch, list):
                    switch = [switch]

            #print(f'switch and opts for {name}: {switch} --- {opts}')
            if True in switch:
                print(f"Connecting {name}...\n")
                for pipe_idx, strat_pool in rpool.get_strats(inputs).items():         # strat_pool is a ResourcePool like {'desc-preproc_T1w': { 'json': info, 'data': (node, out) }, 'desc-brain_mask': etc.}
                    fork = False in switch                                            #   keep in mind rpool.get_strats(inputs) = {pipe_idx1: {'desc-preproc_T1w': etc.}, pipe_idx2: {..} }
                    for opt in opts:                                            #   it's a dictionary of ResourcePools called strat_pools, except those sub-ResourcePools only have one level! no pipe_idx strat keys.
                        # remember, you can get 'data' or 'json' from strat_pool with member functions
                        # strat_pool has all of the JSON information of all the inputs!
                        # so when we set_data below for the TOP-LEVEL MAIN RPOOL (not the strat_pool), we can generate new merged JSON information for each output.
                        #    particularly, our custom 'CpacProvenance' field.
                        node_name = name
                        pipe_x = rpool.get_pipe_number(pipe_idx)
                        wf, outs = block_function(wf, cfg, strat_pool,
                                                  pipe_x, opt)

                        if not outs:
                            continue

                        if opt and len(option_val) > 1:
                            node_name = f'{node_name}_{opt}'
                        elif opt and 'USER-DEFINED' in option_val:
                            node_name = f'{node_name}_{opt["Name"]}'

                        if cfg.pipeline_setup['Debugging']['verbose']:
                            print('\n=======================')
                            print(f'Node name: {node_name}')
                            prov_dct = \
                                rpool.get_resource_strats_from_prov(ast.literal_eval(pipe_idx))
                            for key, val in prov_dct.items():
                                print('-------------------')
                                print(f'Input - {key}:')
                                sub_prov_dct = \
                                    rpool.get_resource_strats_from_prov(val)
                                for sub_key, sub_val in sub_prov_dct.items():
                                    sub_sub_dct = \
                                    rpool.get_resource_strats_from_prov(sub_val)
                                    print(f'  sub-input - {sub_key}:')
                                    print(f'    prov = {sub_val}')
                                    print(f'    sub_sub_inputs = {sub_sub_dct.keys()}')

                        for label, connection in outs.items():
                            self.check_output(outputs, label, name)
                            new_json_info = copy.deepcopy(strat_pool.get('json'))
                            new_json_info['Sources'] = [x for x in strat_pool.get_entire_rpool() if x != 'json']

                            if strat_pool.check_rpool(label):
                                # so we won't get extra forks if we are
                                # merging strats (multiple inputs) plus the
                                # output name is one of the input names
                                old_pipe_prov = list(strat_pool.get_cpac_provenance(label))
                                new_json_info['CpacProvenance'] = old_pipe_prov
                                pipe_idx = strat_pool.generate_prov_string(old_pipe_prov)[1]

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

                            if rpool.func_reg:
                                wf = rpool.derivative_xfm(wf, label,
                                                          connection,
                                                          new_json_info,
                                                          pipe_idx,
                                                          pipe_x)

                            wf = rpool.post_process(wf, label, connection,
                                                    new_json_info, pipe_idx,
                                                    pipe_x, outs)

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
        new_pool = copy.deepcopy(strat_pool)
        for input, val in interface.items():
            if isinstance(val, tuple):
                new_pool.set_data(input, val[0], val[1], {}, "", "")
        wf, outputs = block(wf, cfg, new_pool, pipe_num, opt)
        for out, val in outputs.items():
            if out in interface and isinstance(val, str):
                strat_pool.set_data(out, outputs[out][0], outputs[out][1],
                                    {}, "", "")

    return (wf, strat_pool)


def ingress_raw_anat_data(wf, rpool, cfg, data_paths, unique_id, part_id,
                          ses_id):
    if 'creds_path' not in data_paths:
        data_paths['creds_path'] = None

    anat_flow = create_anat_datasource(f'anat_gather_{part_id}_{ses_id}')
    anat_flow.inputs.inputnode.set(
        subject=part_id,
        anat=data_paths['anat'],
        creds_path=data_paths['creds_path'],
        dl_dir=cfg.pipeline_setup['working_directory']['path'],
        img_type='anat'
    )
    rpool.set_data('T1w', anat_flow, 'outputspec.anat', {},
                   "", "anat_ingress")

    return rpool


def ingress_raw_func_data(wf, rpool, cfg, data_paths, unique_id, part_id,
                          ses_id):

    func_paths_dct = data_paths['func']

    func_wf = create_func_datasource(func_paths_dct,
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
    rpool.set_data('scan_params', func_wf, 'outputspec.scan_params', {}, "",
                   "scan_params_ingress")

    wf, rpool, diff, blip, fmap_rp_list = \
        ingress_func_metadata(wf, cfg, rpool, data_paths, part_id,
                              data_paths['creds_path'], ses_id)

    return (wf, rpool, diff, blip, fmap_rp_list)


def ingress_output_dir(cfg, rpool, unique_id, creds_path=None):

    out_dir = cfg.pipeline_setup['output_directory']['path']

    if cfg.pipeline_setup['output_directory']['pull_source_once']:
        if os.path.isdir(cfg.pipeline_setup['output_directory']['path']):
            if not os.listdir(cfg.pipeline_setup['output_directory']['path']):
                if cfg.pipeline_setup['output_directory']['source_outputs_dir']:
                    out_dir = cfg.pipeline_setup['output_directory'][
                        'source_outputs_dir']
                else:
                    out_dir = cfg.pipeline_setup['output_directory']['path']
            else:
                out_dir = cfg.pipeline_setup['output_directory']['path']
        else:
            if cfg.pipeline_setup['output_directory']['source_outputs_dir']:
                out_dir = cfg.pipeline_setup['output_directory'][
                    'source_outputs_dir']
    else:
        if cfg.pipeline_setup['output_directory']['source_outputs_dir']:
            out_dir = cfg.pipeline_setup['output_directory'][
                'source_outputs_dir']
        else:
            out_dir = cfg.pipeline_setup['output_directory']['path']

    if os.path.isdir(out_dir):
        if not os.listdir(out_dir):
            print(f"\nOutput directory {out_dir} does not exist yet, "
                  f"initializing.")
            return rpool
    else:
        print(f"\nOutput directory {out_dir} does not exist yet, "
              f"initializing.")
        return rpool

    print(f"\nPulling outputs from {out_dir}.\n")

    cpac_dir = os.path.join(out_dir,
                            f'cpac_{cfg.pipeline_setup["pipeline_name"]}',
                            unique_id)
    cpac_dir_anat = os.path.join(cpac_dir, 'anat')
    cpac_dir_func = os.path.join(cpac_dir, 'func')

    exts = ['.nii', '.gz', '.mat', '.1D', '.txt', '.csv', '.rms']

    all_output_dir = []
    if os.path.isdir(cpac_dir_anat):
        for filename in os.listdir(cpac_dir_anat):
            for ext in exts:
                if ext in filename:
                    all_output_dir.append(os.path.join(cpac_dir_anat,
                                                       filename))
    if os.path.isdir(cpac_dir_func):
        for filename in os.listdir(cpac_dir_func):
            for ext in exts:
                if ext in filename:
                    all_output_dir.append(os.path.join(cpac_dir_func,
                                                       filename))

    for filepath in all_output_dir:
        filename = str(filepath)
        for ext in exts:
            filename = filename.split("/")[-1].replace(ext, '')
        data_label = filename.split(unique_id)[1].lstrip('_')

        if len(filename) == len(data_label):
            raise Exception('\n\n[!] Possibly wrong participant or '
                            'session in this directory?\n\n'
                            f'Filepath: {filepath}\n\n')

        if 'task-' in data_label:
            for tag in data_label.split('_'):
                if 'task-' in tag:
                    break
            runtag = None
            if 'run-' in data_label:
                for runtag in data_label.split('_'):
                    if 'run-' in runtag:
                        break
            data_label = data_label.replace(f'{tag}_', '')
            if runtag:
                data_label = data_label.replace(f'{runtag}_', '')

        unique_data_label = str(data_label)

        #if 'sub-' in data_label or 'ses-' in data_label:
        #    raise Exception('\n\n[!] Possibly wrong participant or '
        #                    'session in this directory?\n\nDirectory: '
        #                    f'{cpac_dir_anat}\nFilepath: {filepath}\n\n')
        suffix = data_label.split('_')[-1]
        for tag in data_label.split('_'):
            if 'desc-' in tag:
                desc_val = tag
                break
        jsonpath = str(filepath)
        for ext in exts:
            jsonpath = jsonpath.replace(ext, '')
        jsonpath = f"{jsonpath}.json"

        if not os.path.exists(jsonpath):
            raise Exception('\n\n[!] No JSON found for file '
                            f'{filepath}.\n\n')

        json_info = read_json(jsonpath)
        if 'CpacProvenance' in json_info:
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
                                    f'developer:\n\n{unique_data_label}\n')

                # remove the integer at the end of the desc-* variant, we will get
                # the unique pipe_idx from the CpacProvenance below
                data_label = data_label.replace(desc_val, only_desc)

            # preserve cpac provenance/pipe_idx
            pipe_idx = rpool.generate_prov_string(json_info['CpacProvenance'])
            node_name = ""
        else:
            pipe_idx = ""
            node_name = f"{data_label}_ingress"

        resource = data_label

        ingress = create_general_datasource(f'gather_{unique_data_label}')
        ingress.inputs.inputnode.set(
            unique_id=unique_id,
            data=filepath,
            creds_path=creds_path,
            dl_dir=cfg.pipeline_setup['working_directory']['path']
        )
        rpool.set_data(resource, ingress, 'outputspec.data', json_info,
                       pipe_idx, node_name, inject=True)

    return rpool


def ingress_pipeconfig_paths(cfg, rpool, unique_id, creds_path=None):
    # ingress config file paths
    # TODO: pull this from some external list instead
    # TODO: nah, even better: just loop through the config for .nii's
    # TODO: may want to change the resource keys for each to include one level up in the YAML as well

    templates_for_resampling = [
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.registration_workflows['anatomical_registration']['T1w_brain_template'], 'T1w_brain_template', 'resolution_for_anat'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.registration_workflows['anatomical_registration']['T1w_template'], 'T1w_template', 'resolution_for_anat'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg'], 'T1w_brain_template_funcreg', 'func_preproc_outputs'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg'], 'T1w_brain_template_deriv', 'func_derivative_outputs'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg'], 'T1w_template_funcreg', 'func_preproc_outputs'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg'], 'T1w_template_deriv', 'func_derivative_outputs'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_brain_template_symmetric'], 'T1w_brain_template_symmetric', 'resolution_for_anat'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg'], 'T1w_brain_template_symmetric_deriv', 'func_derivative_outputs'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_template_symmetric'], 'T1w_template_symmetric', 'resolution_for_anat'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg'], 'T1w_template_symmetric_deriv', 'func_derivative_outputs'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['dilated_symmetric_brain_mask'], 'template_dilated_symmetric_brain_mask', 'resolution_for_anat'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['ref_mask'], 'template_ref_mask', 'resolution_for_anat'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg'], 'T1w_brain_template_funcreg', 'func_preproc_outputs'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg'], 'T1w_template_funcreg', 'func_preproc_outputs'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['EPI_template']['EPI_template_funcreg'], 'EPI_template_funcreg', 'func_preproc_outputs'),  # no difference of skull and only brain
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['EPI_template']['EPI_template_funcreg'], 'EPI_template_deriv', 'func_derivative_outputs'),  # no difference of skull and only brain
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg'], 'T1w_brain_template_deriv', 'func_derivative_outputs'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg'], 'T1w_template_deriv', 'func_derivative_outputs')
    ]

    if cfg.PyPEER['run']:
        templates_for_resampling.append((cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.PyPEER['eye_mask_path'], 'template_eye_mask', 'func_preproc_outputs'))
        #Outputs.any.append("template_eye_mask")

    # update resampled template to resource pool
    for resolution, template, template_name, tag in templates_for_resampling:

        if '$FSLDIR' in template:
            template = template.replace('$FSLDIR', cfg.pipeline_setup[
                'system_config']['FSLDIR'])
        #if '${resolution_for_anat}' in template:
        #    template = template.replace('${resolution_for_anat}',
        #                                cfg.registration_workflows[
        #                                    'anatomical_registration'][
        #                                    'resolution_for_anat'])
        if '${func_resolution}' in template:
            template = template.replace('func_resolution', tag)

        resampled_template = pe.Node(Function(input_names=['resolution',
                                                           'template',
                                                           'template_name',
                                                           'tag'],
                                              output_names=['resampled_template'],
                                              function=resolve_resolution,
                                              as_module=True),
                                     name='resampled_' + template_name)

        resampled_template.inputs.resolution = resolution
        resampled_template.inputs.template = template
        resampled_template.inputs.template_name = template_name
        resampled_template.inputs.tag = tag

        # the set_data below is set up a little differently, because we are
        # injecting and also over-writing already-existing entries
        #   other alternative would have been to ingress into the
        #   resampled_template node from the already existing entries, but we
        #   didn't do that here
        rpool.set_data(template_name,
                       resampled_template,
                       'resampled_template',
                       #{'CpacProvenance': [f'{template_name}:{template_name}_config_ingress']},
                       #f"['{template_name}:{template_name}_config_ingress']",
                       {}, "",
                       "template_resample") #, inject=True)   # pipe_idx (after the blank json {}) should be the previous strat that you want deleted! because you're not connecting this the regular way, you have to do it manually

    config_resource_paths = [
        ('CSF_path', cfg.segmentation['tissue_segmentation']['FSL-FAST']['use_priors']['CSF_path']),
        ('WM_path', cfg.segmentation['tissue_segmentation']['FSL-FAST']['use_priors']['WM_path']),
        ('GM_path', cfg.segmentation['tissue_segmentation']['FSL-FAST']['use_priors']['GM_path']),
        ('T1w_ACPC_template', cfg.anatomical_preproc['acpc_alignment']['T1w_ACPC_template']),
        ('T1w_brain_ACPC_template', cfg.anatomical_preproc['acpc_alignment']['T1w_brain_ACPC_template']),
        ('unet_model', cfg.anatomical_preproc['brain_extraction']['UNet']['unet_model']),
        ('T1w_brain_template', cfg.registration_workflows['anatomical_registration']['T1w_brain_template']),
        ('T1w_template', cfg.registration_workflows['anatomical_registration']['T1w_template']),
        ('T1w_brain_template_mask', cfg.registration_workflows['anatomical_registration']['T1w_brain_template_mask']),
        ('T1w_brain_template_symmetric', cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_brain_template_symmetric']),
        ('T1w_template_symmetric', cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_template_symmetric']),
        ('dilated_symmetric_brain_mask', cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['dilated_symmetric_brain_mask']),
        ('T1w_brain_template_symmetric_for_resample', cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_brain_template_symmetric_for_resample']),
        ('T1w_template_symmetric_for_resample', cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_template_symmetric_for_resample']),
        ('dilated_symmetric_brain_mask_for_resample', cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['dilated_symmetric_brain_mask_for_resample']),
        ('ref_mask', cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['ref_mask']),
        ('T1w_template_for_resample', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_for_resample']),
        ('EPI_template_for_resample', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['EPI_template']['EPI_template_for_resample']),
        ('T1w_brain_template_funcreg', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg']),
        ('T1w_brain_template_deriv', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg']),
        ('T1w_template_funcreg', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg']),
        ('T1w_template_deriv', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg']),
        ('T1w_brain_template_symmetric_deriv', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template_funcreg']),
        ('T1w_template_symmetric_deriv', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template_funcreg']),
        ('EPI_template_funcreg', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['EPI_template']['EPI_template_funcreg']),
        ('EPI_template_deriv', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['EPI_template']['EPI_template_funcreg']),
        ('EPI_template', cfg.registration_workflows['functional_registration']['EPI_registration']['EPI_template']),
        ('EPI_template_mask', cfg.registration_workflows['functional_registration']['EPI_registration']['EPI_template_mask']),
        ('lateral_ventricles_mask', cfg.nuisance_corrections['2-nuisance_regression']['lateral_ventricles_mask']),
        ('template_specification_file', cfg.network_centrality['template_specification_file'])
    ]

    if cfg.PyPEER['run']:
        config_resource_paths.append(
            ('eye_mask_path', cfg.PyPEER['eye_mask_path']))

    for resource in config_resource_paths:
        key = resource[0]
        val = resource[1]

        if rpool.check_rpool(key):
            continue

        if not val:
            continue

        if '$FSLDIR' in val:
            val = val.replace('$FSLDIR', cfg.pipeline_setup['system_config']['FSLDIR'])
        if '$priors_path' in val:
            priors_path = cfg.segmentation['tissue_segmentation']['FSL-FAST']['use_priors']['priors_path']
            if '$FSLDIR' in priors_path:
                priors_path = priors_path.replace('$FSLDIR', cfg.pipeline_setup['system_config']['FSLDIR'])
            val = val.replace('$priors_path', priors_path)
        if '${resolution_for_anat}' in val:
            val = val.replace('${resolution_for_anat}', cfg.registration_workflows['anatomical_registration']['resolution_for_anat'])
        if '${func_resolution}' in val:
            # functional registration
            if 'funcreg' in key:
                out_res = 'func_preproc_outputs'
            # functional derivatives
            else:
                out_res = 'func_derivative_outputs'
            val = val.replace('${func_resolution}', cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution'][out_res])

        if val:
            config_ingress = create_general_datasource(f'gather_{key}')
            config_ingress.inputs.inputnode.set(
                unique_id=unique_id,
                data=val,
                creds_path=creds_path,
                dl_dir=cfg.pipeline_setup['working_directory']['path']
            )
            rpool.set_data(key, config_ingress, 'outputspec.data', {}, "",
                           f"{key}_config_ingress")

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
        ("anat", ["anatomical_preproc", "acpc_alignment", "T1w_brain_ACPC_template"])]

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
      {'anat': '{T1w path}',
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
       'unique_id': 'ses-1'}
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
        rpool = ingress_raw_anat_data(wf, rpool, cfg, data_paths, unique_id,
                                      part_id, ses_id)

        wf, rpool, diff, blip, fmap_rp_list = \
            ingress_raw_func_data(wf, rpool, cfg, data_paths, unique_id,
                                  part_id, ses_id)

    # grab already-processed data from the output directory
    rpool = ingress_output_dir(cfg, rpool, unique_id, creds_path)

    # grab any file paths from the pipeline config YAML
    rpool = ingress_pipeconfig_paths(cfg, rpool, unique_id, creds_path)

    return (wf, rpool)


def run_node_blocks(blocks, data_paths, cfg=None):
    import os
    from CPAC.pipeline import nipype_pipeline_engine as pe
    from CPAC.utils.strategy import NodeBlock

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
    rpool = initiate_rpool(cfg, data_paths)

    wf = pe.Workflow(name=f'node_blocks')
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
        wf = NodeBlock(block).connect_block(wf, cfg, rpool)
    rpool.gather_pipes(wf, cfg)

    wf.run()
