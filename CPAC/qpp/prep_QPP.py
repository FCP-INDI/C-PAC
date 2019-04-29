import os
import glob
import shutil
import fnmatch
import itertools
import pandas as pd
import pkg_resources as p
import numpy as np
import nibabel as nib

from CPAC.pipeline.cpac_group_runner import (
    load_text_file,
    grab_pipeline_dir_subs,
    load_config_yml,
    gather_nifti_globs,
    create_output_df_dict,
    gather_outputs
)
from CPAC.pipeline.cpac_ga_model_generator import create_merged_copefile, create_merge_mask
from CPAC.utils import Outputs

def split_subdfs(output_df_dict, sess_inclusion=None, scan_inclusion=None,
                 grp_by_strat=None):

    for unique_resource in output_df_dict.keys():

        resource_id = unique_resource[0]
        strat_info = unique_resource[1]

        output_df = output_df_dict[unique_resource]

        if "Series" in output_df:
            output_df.rename(columns={"Series": "Scan"},
                             inplace=True)

        join_columns = ["participant_id"]
        col_names = output_df.columns.tolist()

        # We're then going to reduce the Output directory to contain only those scans and or the sessions which are expressed by the user.
        # If the user answers all to the option, then we're obviously not going to do any repeated measures.
        # add a qpp dict so that you don't make stupid af errors again!
        qpp_dict = {}

        grp_by_sessions = False
        grp_by_scans = False
        grp_by_both = False
        repeated_measures = False

        if sess_inclusion:
            # multiple sessions, so you're going group by scans
            if len(sess_inclusion) > 0:
                grp_by_scans = True

        # Multiple scans so you're going to group by sessions
        if scan_inclusion:
            if len(scan_inclusion) > 0:
                grp_by_sessions = True

        if grp_by_scans or grp_by_sessions:
            repeated_measures = True


        if grp_by_scans and grp_by_sessions:
            grp_by_both = True

        print(repeated_measures)
        # PSA #both multiple sessions and scans, youre going to do nothing
        # if neither, the output directory will not have any level of grouping, it will be ses1_scan1 --> nuisance strat, etc
        if repeated_measures:

            if grp_by_scans:
                new_output_df = output_df[output_df["Sessions"].isin(sess_inclusion)]
                join_columns.append("Sessions")
            if grp_by_sessions:
                # drop all the scans that are not in the scan list
                new_output_df = output_df[output_df["Scan"].isin(scan_inclusion)]
                join_columns.append("Scan")
            newer_output_df, dropped_parts = balance_df(new_output_df,sess_inclusion,scan_inclusion)


        pre_qpp_dict = {}

        if grp_by_strat and repeated_measures == True:

            if grp_by_strat == ["Scan"]:
                for scan_df_tuple in newer_output_df.groupby("Scan"):
                    scans = scan_df_tuple[0]
                    scan_df = scan_df_tuple[1]
                    if 'Sessions' in scan_df.columns:
                        for ses_df_tuple in scan_df.groupby('Sessions'):
                            session = 'ses-{0}'.format(ses_df_tuple[0])
                            ses_df = ses_df_tuple[1]
                            if not ses_df_tuple[0] in qpp_dict:
                                pre_qpp_dict[ses_df_tuple[0]] = []
                            pre_qpp_dict[ses_df_tuple[0]].append(ses_df)
                        for key in pre_qpp_dict:
                            list_df = pre_qpp_dict[key]
                            concat_df = list_df[0]
                            for df in list_df[1:]:
                                concat_df = pd.concat(concat_df, df)
                            qpp_dict[key] = concat_df

            elif grp_by_strat == ["Session"]:
                session = 'ses-1'
                qpp_dict['ses-1'] = scan_df
            elif grp_by_strat == None:
                print(newer_output_df)
                qpp_dict['output_df'] = newer_output_df

        if repeated_measures == True and grp_by_strat == ['None']:
            print(output_df)
            qpp_dict['output_df'] = newer_output_df
        
        if repeated_measures == False and grp_by_strat == ['None']:
            print(output_df)
            qpp_dict['output_df'] = output_df

        if len(qpp_dict) == 0:
            err = '\n\n[!]C-PAC says:Could not find match betterrn the ' \
                  'particpants in your pipeline output directory that were ' \
                  'included in your analysis, and the particpants in the phenotype ' \
                  'phenotype file provided.\n\n'
            raise Exception(err)
    print(qpp_dict)
    return qpp_dict,resource_id,strat_info


def prep_inputs(group_config_file):
    keys_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.csv')
    try:
        keys = pd.read_csv(keys_csv)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.csv " \
              "resource file:\n{0}\n\nError details {1}\n".format(keys_csv, e)
        raise Exception(err)

    group_config_obj = load_config_yml(group_config_file)
    pipeline_folder = group_config_obj.pipeline_dir

    try:
        grp_by_strat = group_config_obj.qpp_grpby_strat
    except AttributeError:
        grp_by_strat = False

    if not group_config_obj.participant_list:
        inclusion_list = grab_pipeline_dir_subs(pipeline_folder)
    elif '.' in group_config_obj.participant_list:

        if not os.path.isfile(group_config_obj.participant_list):

            raise Exception('\n[!] C-PAC says: Your participant '
                            'inclusion list is not a valid file!\n\n'
                            'File path: {0}'
                            '\n'.format(group_config_obj.participant_list))
        else:
            inclusion_list = load_text_file(group_config_obj.participant_list, "group-level analysis participant list")

    else:
        inclusion_list = grab_pipeline_dir_subs(pipeline_folder)

    resource_list = ['functional_nuisance_residuals']
    output_df_dict = gather_outputs(pipeline_folder, resource_list, inclusion_list)

    if not output_df_dict:
        err = '\n\n[!] For QPP, C-PAC requires the \'functional_nuisance_residuals\' outputs ' \
              'in the individual-level analysis pipeline output directory! But none were ' \
              'found.\n\nPipeline directory:\n{0}\n\n'.format(pipeline_folder)
        raise Exception(err)

    qpp_dict,resource_id,strat_info = split_subdfs(output_df_dict, group_config_obj.qpp_sess_inclusion,
                            group_config_obj.qpp_scan_inclusion, grp_by_strat)

    return qpp_dict, inclusion_list, resource_id, strat_info


def use_inputs(group_config_file):


    qpp_dict, inclusion_list, resource_id, strat_info = prep_inputs(group_config_file)
    group_config_obj = load_config_yml(group_config_file)

    for key in qpp_dict:

        newer_output_df = qpp_dict[key]

        pipeline_ID = group_config_obj.pipeline_dir.rstrip('/').split('/')[-1]

        model_dir = os.path.join(group_config_obj.output_dir,
                                 'cpac_group_analysis',
                                 'CPAC_QPP')
        out_dir = os.path.join(model_dir,
                               resource_id,  # nuisance strat to initialize
                               strat_info, 'model_files')

        participant_id = newer_output_df['participant_id'].tolist()
        for element in participant_id:
            if "Scan" in newer_output_df.columns:
                scan_list = newer_output_df["Scan"].tolist()
                nrn = len(set(scan_list))
            elif "Session" in newer_output_df.columns:
                session_list = newer_output_df["Session"].tolist()
                nrn = len(set(session_list))
            else:
                nrn = 1

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        subject_list = newer_output_df["Filepath"].tolist()
        merge_outfile = os.path.join(out_dir, '_merged.nii.gz')
        merge_file = create_merged_copefile(newer_output_df["Filepath"].tolist(), merge_outfile)
        merge_mask_outfile = os.path.join(out_dir, '_merged_mask.nii.gz')
        merge_mask = create_merge_mask(merge_file, merge_mask_outfile)

    return merge_file, merge_mask, subject_list, inclusion_list, out_dir, nrn


def balance_df(new_output_df, session_list, scan_list):
    import pandas as pd
    from collections import Counter

    part_ID_count = Counter(new_output_df["participant_id"])


    if scan_list and session_list:
        sessions_x_scans = len(session_list) * len(scan_list)
    elif session_list:
        sessions_x_scans = len(session_list)
    elif scan_list:
        sessions_x_scans = len(scan_list)
    dropped_parts = []
    for part_ID in part_ID_count.keys():
        if part_ID_count[part_ID] != sessions_x_scans:
            newer_output_df = new_output_df[new_output_df.participant_id != part_ID]
            dropped_parts.append(part_ID)
        else:
            newer_output_df=new_output_df

    return newer_output_df, dropped_parts
