import os
import six
import json
import warnings
import logging

import nipype.pipeline.engine as pe
import nipype.interfaces.utility
from nipype.interfaces.utility import Rename
from CPAC.utils.interfaces.function import Function
from CPAC.utils.interfaces.datasink import DataSink

from CPAC.utils.utils import read_json, create_id_string, write_output_json, \
    get_last_prov_entry
from CPAC.utils.datasource import (
    create_anat_datasource,
    create_func_datasource,
    ingress_func_metadata,
    create_general_datasource,
    create_check_for_s3_node,
    resolve_resolution
)

logger = logging.getLogger('workflow')


class ResourcePool(object):
    def __init__(self, rpool=None, name=None):
        if not rpool:
            self.rpool = {}
        else:
            self.rpool = rpool
        self.pipe_list = []
        self.name = name

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

    def set_json_info(self, resource, pipe_idx, key, val):
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

    def set_data(self, resource, node, output, json_info, pipe_idx, node_name,
                 fork=False, inject=False):
        cpac_prov = []
        if 'CpacProvenance' in json_info:
            cpac_prov = json_info['CpacProvenance']
        new_prov_list = list(cpac_prov)   # <---- making a copy
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
                print(f'old pipe idx: {pipe_idx}')
                print(self.rpool[resource].keys())
                if pipe_idx in self.rpool[resource].keys():  # <--- in case the resource name is now new, and not the original
                    del self.rpool[resource][pipe_idx]  # <--- remove old keys so we don't end up with a new strat for every new node unit (unless we fork)

        if new_pipe_idx not in self.rpool[resource]:
            self.rpool[resource][new_pipe_idx] = {}
        if new_pipe_idx not in self.pipe_list:
            self.pipe_list.append(new_pipe_idx)
        self.rpool[resource][new_pipe_idx]['data'] = (node, output)
        self.rpool[resource][new_pipe_idx]['json'] = json_info

    def get(self, resource, report_fetched=False, optional=False):
        # NOTE!!!
        #   if this is the main rpool, this will return a dictionary of strats, and inside those, are dictionaries like {'data': (node, out), 'json': info}
        #   BUT, if this is a sub rpool (i.e. a strat_pool), this will return a one-level dictionary of {'data': (node, out), 'json': info} WITHOUT THE LEVEL OF STRAT KEYS ABOVE IT
        if isinstance(resource, list):
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
                raise Exception("\n\n[!] C-PAC says: The listed resource is "
                                f"not in the resource pool:\n{resource}\n\n"
                                "Developer Note: This may be due to a mis"
                                "match between the node block's docstring "
                                "'input' field and a strat_pool.get_data() "
                                "call within the block function.\n")
            if report_fetched:
                return (self.rpool[resource], resource)
            return self.rpool[resource]

    def get_data(self, resource, report_fetched=False):
        return self.get(resource, report_fetched)['data']

    def copy_resource(self, resource, new_name):
        self.rpool[new_name] = self.rpool[resource]

    def get_json(self, resource, strat=None):
        # NOTE: resource_strat_dct has to be entered properly by the developer
        # it has to either be rpool[resource][strat] or strat_pool[resource]
        if strat:
            resource_strat_dct = self.rpool[resource][strat]
        else:
            # for strat_pools mainly, where there is no 'strat' key level
            resource_strat_dct = self.rpool[resource]

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

    def update_cpac_provenance(self, new_resource, prov, inputs, new_node_id):
        if new_resource not in prov:
            prov[new_resource] = [inputs]
        prov[new_resource].append(new_node_id)
        self.set_json_info(new_resource, )

    def get_strats(self, resources):

        # TODO: NOTE: NOT COMPATIBLE WITH SUB-RPOOL/STRAT_POOLS
        # TODO: (and it doesn't have to be)

        import itertools

        total_pool = []
        len_inputs = len(resources)
        for resource in resources:
            rp_dct, fetched_resource = self.get(resource, True, True)   # <---- rp_dct has the strats/pipe_idxs as the keys on first level, then 'data' and 'json' on each strat level underneath
                                                                        # oh, and we make the resource fetching in get_strats optional so we can have optional inputs, but they won't be optional in the node block unless we want them to be
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

            # we now currently have "strats", the combined permutations of all the strategies, as a list of tuples, each tuple being a permutation.
            # OF ALL THE DIFFERENT INPUTS. and they are tagged by their fetched inputs with {name}:{strat}.
            # so, each tuple has ONE STRAT FOR EACH INPUT, so if there are three inputs, each tuple will have 3 items.

            new_strats = {}
            for strat_tuple in strats:
                strat_list = list(strat_tuple)     # <------- strat_list is now a list of combined strats, one of the permutations. keep in mind each strat in the combo comes from a different data source/input
                                                   #          like this:   strat_list = [desc-preproc_T1w:pipe_idx_anat_1, desc-brain_mask:pipe_idx_mask_1]
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

    def gather_pipes(self, wf, cfg):
        # TODO: cpac_outputs.csv etc
        # TODO: might be better to do an inclusion instead
        excl = ['T1w', 'bold', 'scan', 'scan_params', 'TR', 'tpattern',
                'start_tr', 'stop_tr', 'pe_direction', 'subject',
                'motion_basefile']
        config_paths = ['T1w_ACPC_template', 'T1w_brain_ACPC_template',
                        'unet_model', 'T1w_brain_template', 'T1w_template',
                        'T1w_brain_template_mask',
                        'T1w_brain_template_symmetric',
                        'T1w_template_symmetric',
                        'dilated_symmetric_brain_mask',
                        'T1w_brain_template_symmetric_for_resample',
                        'T1w_template_symmetric_for_resample', 'ref_mask',
                        'template_for_resample',
                        'T1w_brain_template_for_func',
                        'T1w_template_for_func',
                        'template_epi', 'template_epi_mask',
                        'lateral_ventricles_mask', 'eye_mask_path']
        excl += config_paths
        anat = ['T1w', 'probseg']
        func = ['bold']
        motions = ['motion', 'movement', 'coordinate', 'displacement']

        for resource in self.rpool.keys():
            # TODO: cpac_outputs.csv etc
            if resource in excl:
                continue

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
                subdir = 'other'
                for tag in motions:
                    if tag in resource:
                        subdir = 'func'

            for pipe_idx in self.rpool[resource]:
                pipe_num = self.get_pipe_number(pipe_idx)
                unique_id = self.get_name()

                out_dir = cfg.pipeline_setup['output_directory']['path']
                container = os.path.join('cpac', unique_id)
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

            num_variant = 0
            if len(self.rpool[resource]) == 1:
                num_variant = ""
            for pipe_idx in self.rpool[resource]:
                try:
                    num_variant += 1
                except TypeError:
                    pass
                json_info = self.rpool[resource][pipe_idx]['json']
                out_dct = self.rpool[resource][pipe_idx]['out']

                unique_id = out_dct['unique_id']

                for key in out_dct['filename'].split('_'):
                    if 'desc-' in key:
                        out_dct['filename'] = out_dct['filename'
                        ].replace(key, f'{key}{num_variant}')
                        resource_idx = resource.replace(key,
                                                        f'{key}{num_variant}')
                        break
                    else:
                        resource_idx = f'{out_dct["filename"]}{num_variant}'

                id_string = pe.Node(Function(input_names=['unique_id',
                                                          'resource',
                                                          'scan_id'],
                                             output_names=['out_filename'],
                                             function=create_id_string),
                                    name=f'id_string_{resource_idx}')
                id_string.inputs.unique_id = unique_id
                id_string.inputs.resource = resource_idx

                if out_dct['subdir'] == 'func':
                    node, out = self.rpool['scan']["['scan:func_ingress']"][
                        'data']
                    wf.connect(node, out, id_string, 'scan_id')

                nii_name = pe.Node(Rename(), name=f'{resource_idx}')
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
                                     name=f'json_{resource_idx}')
                write_json.inputs.json_data = json_info

                wf.connect(id_string, 'out_filename', write_json, 'filename')

                ds = pe.Node(DataSink(), name=f'sinker_{resource_idx}')
                ds.inputs.parameterization = False
                ds.inputs.base_directory = out_dct['out_dir']
                ds.inputs.encrypt_bucket_keys = cfg.pipeline_setup['Amazon-AWS']['s3_encryption']
                ds.inputs.container = out_dct['container']

                if cfg.pipeline_setup['Amazon-AWS']['aws_output_bucket_credentials']:
                    ds.inputs.creds_path = cfg.pipeline_setup['Amazon-AWS']['aws_output_bucket_credentials']

                wf.connect(nii_name, 'out_file', ds, f'{out_dct["subdir"]}.@data')
                wf.connect(write_json, 'json_file', ds, f'{out_dct["subdir"]}.@json')

