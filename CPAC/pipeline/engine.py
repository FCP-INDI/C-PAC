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
        import json
        init_dct_schema = ['name', 'config', 'switch', 'option_key',
                           'option_val', 'inputs', 'outputs']
        try:
            dct = json.loads(fn_docstring.replace('\n', '').replace(' ', ''))
        except Exception as e:
            raise Exception('\n\n[!] Node block docstring error.\n\n'
                            f'Docstring:\n{fn_docstring}\n\n')
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
        for name, block_dct in self.node_blocks.items():    # <--- iterates over either the single node block in the sequence, or a list of node blocks within the list of node blocks, i.e. for option forking.

            switch = self.check_null(block_dct['switch'])
            config = self.check_null(block_dct['config'])
            option_key = self.check_null(block_dct['option_key'])
            option_val = self.check_null(block_dct['option_val'])
            inputs = self.check_null(block_dct['inputs'])
            outputs = self.check_null(block_dct['outputs'])

            block_function = block_dct['block_function']

            # NOTE: this has been moved up to the node-block-pipeline level
            #exists = [rpool.check_rpool(out) for out in outputs]
            #if False not in exists:
            #    print(f'{outputs} already exist, bypassing {name}')
            #    return wf

            opts = []
            if option_key and option_val:
                if not isinstance(option_key, list):
                    option_key = [option_key]
                if not isinstance(option_val, list):
                    option_val = [option_val]
                key_list = config + option_key
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
                try:
                    key_list = config + switch
                except TypeError:
                    raise Exception("\n\n[!] Developer info: Docstring error "
                                    f"for {name}, make sure the 'config' or "
                                    "'switch' fields are lists.\n\n")
                switch = self.grab_tiered_dct(cfg, key_list)
                if not isinstance(switch, list):
                    switch = [switch]

            print(f'switch and opts for {name}: {switch} --- {opts}')

            if True in switch:
                for pipe_idx, strat_pool in rpool.get_strats(inputs).items():         # strat_pool is a ResourcePool like {'desc-preproc_T1w': { 'json': info, 'data': (node, out) }, 'desc-brain_mask': etc.}
                    fork = False in switch                                            #   keep in mind rpool.get_strats(inputs) = {pipe_idx1: {'desc-preproc_T1w': etc.}, pipe_idx2: {..} }
                    for opt in opts:                                            #   it's a dictionary of ResourcePools called strat_pools, except those sub-ResourcePools only have one level! no pipe_idx strat keys.
                        # remember, you can get 'data' or 'json' from strat_pool with member functions
                        # strat_pool has all of the JSON information of all the inputs!
                        # so when we set_data below for the TOP-LEVEL MAIN RPOOL (not the strat_pool), we can generate new merged JSON information for each output.
                        #    particularly, our custom 'CpacProvenance' field.

                        pipe_x = rpool.get_pipe_number(pipe_idx)
                        wf, outs = block_function(wf, cfg, strat_pool,
                                                  pipe_x, opt)

                        if opt and len(option_val) > 1:
                            name = f'{name}_{opt}'

                        for label, connection in outs.items():
                            self.check_output(outputs, label, name)
                            json_info = dict(strat_pool.get('json'))
                            json_info['Sources'] = [x for x in strat_pool.get_entire_rpool() if x != 'json']
                            if strat_pool.check_rpool(label):
                                # so we won't get extra forks if we are
                                # merging strats (multiple inputs) plus the
                                # output name is one of the input names
                                old_pipe_prov = strat_pool.get_cpac_provenance(label)
                                pipe_idx = strat_pool.generate_prov_string(old_pipe_prov)[1]

                            rpool.set_data(label,
                                           connection[0],
                                           connection[1],
                                           json_info,
                                           pipe_idx, name, fork)

        return wf


def initiate_rpool(wf, cfg, data_paths):

    # TODO: refactor further, integrate with the ingress_data functionality
    # TODO: used for BIDS-Derivatives (below), and possible refactoring of
    # TODO: the raw data config to use 'T1w' label instead of 'anat' etc.

    part_id = data_paths['subject_id']
    ses_id = data_paths['unique_id']

    unique_id = f'{part_id}_{ses_id}'

    rpool = ResourcePool(name=unique_id)

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

    func_paths_dict = data_paths['func']

    func_wf = create_func_datasource(func_paths_dict,
                                     f'func_ingress_{part_id}_{ses_id}')
    func_wf.inputs.inputnode.set(
        subject=part_id,
        creds_path=data_paths['creds_path'],
        dl_dir=cfg.pipeline_setup['working_directory']['path']
    )
    func_wf.get_node('inputnode').iterables = \
        ("scan", list(func_paths_dict.keys()))

    rpool.set_data('subject', func_wf, 'outputspec.subject', {}, "",
                   "func_ingress")
    rpool.set_data('bold', func_wf, 'outputspec.rest', {}, "", "func_ingress")
    rpool.set_data('scan', func_wf, 'outputspec.scan', {}, "", "func_ingress")
    rpool.set_data('scan_params', func_wf, 'outputspec.scan_params', {}, "",
                   "scan_params_ingress")

    wf, rpool, diff, blip, fmap_rp_list = \
        ingress_func_metadata(wf, cfg, rpool, data_paths, part_id,
                              data_paths['creds_path'])

    # grab already-processed data from the output directory
    out_dir = cfg.pipeline_setup['output_directory']['path']
    cpac_dir = os.path.join(out_dir, 'cpac', unique_id)
    cpac_dir_anat = os.path.join(cpac_dir, 'anat')
    cpac_dir_func = os.path.join(cpac_dir, 'func')

    all_output_dir = []
    if os.path.isdir(cpac_dir_anat):
        for filename in os.listdir(cpac_dir_anat):
            if '.nii' in filename:
                all_output_dir.append(os.path.join(cpac_dir_anat, filename))
    if os.path.isdir(cpac_dir_func):
        for filename in os.listdir(cpac_dir_func):
            if '.nii' in filename:
                all_output_dir.append(os.path.join(cpac_dir_func, filename))

    for filepath in all_output_dir:
        filename = filepath.split("/")[-1].replace('.nii', '').replace('.gz', '')
        data_label = filename.split(unique_id)[1].lstrip('_')
        if 'sub-' in data_label or 'ses-' in data_label:
            raise Exception('\n\n[!] Possibly wrong participant or '
                            'session in this directory?\n\nDirectory: '
                            f'{cpac_dir_anat}\nFilepath: {filepath}\n\n')
        suffix = data_label.split('_')[-1]
        for tag in data_label.split('_'):
            if 'desc-' in tag:
                desc_val = tag
        jsonpath = filepath.replace('.gz', '').replace('.nii', '.json')
        if not os.path.exists(jsonpath):
            raise Exception('\n\n[!] No JSON found for file '
                            f'{filepath}.\n\n')

        json_info = read_json(jsonpath)
        if 'CpacProvenance' in json_info:
            # it's a C-PAC output, let's check for pipe_idx/strat integer
            # suffixes in the desc- entries.
            for idx in range(0, 3):
                # let's stop at 3, please don't run >999 strategies okay?
                if desc_val[-1].isdigit():
                    desc_val = desc_val[:-1]

            # preserve cpac provenance/pipe_idx
            pipe_idx = rpool.generate_prov_string(json_info['CpacProvenance'])
            node_name = ""
        else:
            pipe_idx = ""
            node_name = f"{data_label}_ingress"

        resource = data_label

        ingress = create_general_datasource(f'gather_{data_label}')
        ingress.inputs.inputnode.set(
            unique_id=unique_id,
            data=filepath,
            creds_path=data_paths['creds_path'],
            dl_dir=cfg.pipeline_setup['working_directory']['path']
        )
        rpool.set_data(resource, ingress, 'outputspec.data', json_info,
                       pipe_idx, node_name, inject=True)

    # ingress config file paths
    # TODO: pull this from some external list instead
    # TODO: nah, even better: just loop through the config for .nii's
    # TODO: may want to change the resource keys for each to include one level up in the YAML as well
    config_resource_paths = [
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
        ('template_for_resample', cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['template_for_resample']),
        ('T1w_brain_template_for_func', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1w_brain_template']),
        ('T1w_template_for_func', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1w_template']),
        ('template_epi', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['template_epi']),
        ('template_epi_mask', cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['template_epi_mask']),
        ('lateral_ventricles_mask', cfg.nuisance_corrections['2-nuisance_regression']['lateral_ventricles_mask'])
    ]

    if cfg.PyPeer['run']:
        config_resource_paths.append(('eye_mask_path', cfg.PyPeer['eye_mask_path']))

    for resource in config_resource_paths:
        key = resource[0]
        val = resource[1]

        if not val:
            continue

        if '$FSLDIR' in val:
            val = val.replace('$FSLDIR', cfg.pipeline_setup['system_config']['FSLDIR'])
        if '${resolution_for_anat}' in val:
            val = val.replace('${resolution_for_anat}', cfg.registration_workflows['anatomical_registration']['resolution_for_anat'])

        if val:
            config_ingress = create_general_datasource(f'gather_{key}')
            config_ingress.inputs.inputnode.set(
                unique_id=unique_id,
                data=val,
                creds_path=data_paths['creds_path'],
                dl_dir=cfg.pipeline_setup['working_directory']['path']
            )
            rpool.set_data(key, config_ingress, 'outputspec.data', {}, "",
                           f"{key}_config_ingress")

    # templates, resampling from config
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
        ("other",
         ["voxel_mirrored_homotopic_connectivity", "symmetric_registration",
          "FNIRT_pipelines", "config_file"]),
    ]

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

    templates_for_resampling = [
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.registration_workflows['anatomical_registration']['T1w_brain_template'], 'T1w_brain_template', 'resolution_for_anat'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.registration_workflows['anatomical_registration']['T1w_template'], 'T1w_template', 'resolution_for_anat'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_brain_template_symmetric'], 'T1w_brain_template_symmetric', 'resolution_for_anat'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['T1w_template_symmetric'], 'T1w_template_symmetric', 'resolution_for_anat'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.voxel_mirrored_homotopic_connectivity['symmetric_registration']['dilated_symmetric_brain_mask'], 'template_dilated_symmetric_brain_mask', 'resolution_for_anat'),
        (cfg.registration_workflows['anatomical_registration']['resolution_for_anat'], cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['ref_mask'], 'template_ref_mask', 'resolution_for_anat'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template'], 'template_brain_for_func_preproc', 'resolution_for_func_preproc'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template'], 'template_skull_for_func_preproc', 'resolution_for_func_preproc'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['EPI_template']['template_epi'], 'template_epi', 'resolution_for_func_preproc'),  # no difference of skull and only brain
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['EPI_template']['template_epi'], 'template_epi_derivative', 'resolution_for_func_derivative'),  # no difference of skull and only brain
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_brain_template'], 'template_brain_for_func_derivative', 'resolution_for_func_preproc'),
        (cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_derivative_outputs'], cfg.registration_workflows['functional_registration']['func_registration_to_template']['target_template']['T1_template']['T1w_template'], 'template_skull_for_func_derivative', 'resolution_for_func_preproc'),
    ]

    if True in cfg.PyPEER['run']:
        templates_for_resampling.append((cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'], cfg.PyPEER['eye_mask_path'], 'template_eye_mask', 'resolution_for_func_preproc'))
        #Outputs.any.append("template_eye_mask")

    # update resampled template to resource pool
    for resolution, template, template_name, tag in templates_for_resampling:

        if '$FSLDIR' in template:
            template = template.replace('$FSLDIR', cfg.pipeline_setup[
                'system_config']['FSLDIR'])

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

        rpool.set_data(template_name,
                       resampled_template,
                       'resampled_template', {},
                       f"['{template_name}:{template_name}_config_ingress']",
                       "template_resample")   # pipe_idx (after the blank json {}) should be the previous strat that you want deleted! because you're not connecting this the regular way, you have to do it manually

    return (wf, rpool)