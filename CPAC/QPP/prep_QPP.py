def create_output_dict_list(nifti_globs,pipeline_output_folder,resource_list,search_dirs):
    import os
    import glob
    import fnmatch
    import itertools
    import pandas as pd
    import pkg_resources as p

    if len(resource_list) == 0:
        err= "\n\n[!] No derivatives selected!\n\n"
        raise Exception(err)
    exts=['nii','nii.gz']
    exts = ['.'+ ext.lstrip('.')for ext in exts]
    output_dict_list={}
    for root, _,files in os.walk(pipeline_output_folder):
        for filename in files:
            filepath=os.path.join(root,filename)

            if not any(fnmatch.fnmatch(filepath,pattern)for pattern in nifti_globs):
                continue
            if not any(filepath.endswith(ext)for ext in exts):
                continue
            relative_filepath=filepath.split(pipeline_output_folder)[1]
            filepath_pieces=filter(None, relative_filepath.split("/"))
            resource_id = filepath_pieces[1]

            if resource_id not in search_dirs:
                continue

            scan_id_string = filepath_pieces[2]
            strat_info = "_".join(filepath_pieces[3:])[:-len(ext)]
            unique_resource_id=(resource_id,strat_info)
            if unique_resource_id not in output_dict_list.keys():

                output_dict_list[unique_resource_id] = []

            unique_id = filepath_pieces[0]

            scan_id = scan_id_string.replace("_scan_", "")
            scan_id = scan_id.replace("_rest", "")

            new_row_dict = {}
            new_row_dict["participant_session_id"] = unique_id
            new_row_dict["participant_id"], new_row_dict["Sessions"] = \
                unique_id.split('_')

            new_row_dict["Scan"] = scan_id
            new_row_dict["Filepath"] = filepath

            print('{0} - {1} - {2}'.format(unique_id, scan_id,
                                           resource_id))
            output_dict_list[unique_resource_id].append(new_row_dict)
        #analysis, grouped either by sessions or scans.
    return output_dict_list


def gather_outputs(pipeline_folder,resource_list,inclusion_list):

    from CPAC.pipeline.cpac_group_runner import gather_nifti_globs
    from CPAC.pipeline.cpac_group_runner import create_output_df_dict
    import pandas as pd
    import pkg_resources as p

    keys_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.csv')
    try:
        keys=pd.read_csv(keys_csv)

    except Exception as e:
        err= "\n[!] Could not access or read the cpac_outputs.csv " \
                "resource file:\n{0}\n\nError details {1}\n".format(keys_csv,e)
        raise Exception(err)
    resource_list = ['functional_nuisance_residuals']
    derivatives =list(keys[keys['Space'] == 'functional'][keys['Functional timeseries'] == 'yes']['Resource'])

    nifti_globs, search_dirs = gather_nifti_globs(pipeline_folder, resource_list, derivatives)

    output_dict_list = create_output_dict_list(nifti_globs,pipeline_folder,resource_list,search_dirs)
    # now we have a good dictionary which contains all the filepaths of the files we need to merge later on.
    # Steps after this: 1. This is only a dictionary so let's convert it to a data frame.
    # 2. In the data frame, we're going to only include whatever is in the participant list
    # 3. From the list of included participants, we're going to further prune the output dataframe to only contain the
    # scans included, and/or the sessions included
    # 4. Our final output df will contain, file paths for the .nii files of all the participants that are included in the
    output_df_dict=create_output_df_dict(output_dict_list,inclusion_list)

    return output_df_dict


def prep_inputs(group_config_file):

    import os
    import shutil
    import pandas as pd
    import pkg_resources as p
    from CPAC.pipeline.cpac_group_runner import load_text_file
    from CPAC.pipeline.cpac_group_runner import grab_pipeline_dir_subs
    from CPAC.pipeline.cpac_group_runner import load_config_yml
    from CPAC.pipeline.cpac_ga_model_generator import create_merged_copefile,create_merge_mask

    keys_csv = p.resource_filename('CPAC','resources/cpac_outputs.csv')
    try:
        keys = pd.read_csv(keys_csv)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.csv " \
              "resource file:\n{0}\n\nError details {1}\n".format(keys_csv,e)
        raise Exception(err)

    group_config_obj=load_config_yml(group_config_file)
    pipeline_folder = group_config_obj.pipeline_dir
    #inclusion list function
    if not group_config_obj.participant_list:
        inclusion_list = grab_pipeline_dir_subs(pipeline_folder)
    elif '.' in group_config_obj.participant_list:

        if not os.path.isfile(group_config_obj.participant_list):

            raise Exception('\n[!] C-PAC says: Your participant '
                            'inclusion list is not a valid file!\n\n'
                            'File path: {0}'
                            '\n'.format(group_config_obj.participant_list))
        else:
            inclusion_list = load_text_file(group_config_obj.participant_list,"group-level analysis participant list")

    else:
        inclusion_list = grab_pipeline_dir_subs(pipeline_folder)

    resource_list = ['functional_nuisance_residuals']
    output_df_dict=gather_outputs(pipeline_folder,resource_list,inclusion_list)


    if not output_df_dict:
        err = '\n\n[!] For QPP, C-PAC requires the \'functional_nuisance_residuals\' outputs '\
              'in the individual-level analysis pipeline output directory! But none were '\
              'found.\n\nPipeline directory:\n{0}\n\n'.format(pipeline_folder)
        raise Exception(err)

    for unique_resource in output_df_dict.keys():
        resource_id = unique_resource[0]
        strat_info=unique_resource[1]
        output_df=output_df_dict[unique_resource]
        #We're going to reduce the size of the output df based on nuisance strat and the
        #participant list that actually is included.
        if not group_config_obj.participant_list:
            inclusion_list = grab_pipeline_dir_subs(pipeline_folder)
            output_df = output_df[output_df["participant_id"].isin(inclusion_list)]
        elif '.' in group_config_obj.participant_list:
            if os.path.isfile(group_config_obj.participant_list):
                inclusion_list = load_text_file(group_config_obj.participant_list,
                                        "group-level analysis "
                                        "participant list")
                output_df = output_df[output_df["participant_id"].isin(inclusion_list)]
            else:
                raise Exception('\nCannot read group-level analysis participant ' \
                        'list.\n')

        if len(output_df) == 0:
            err = "\n\n[!]The output data frame has not been compiled correctly. Please " \
                  "recheck the participants in the participant inclusion list and the " \
                  "pipeline output directory to troubleshoot.\n\n"
            raise Exception(err)

        join_columns = ["participant_id"]
        col_names = output_df.columns.tolist()
        # We're then going to reduce the Output directory to contain only those scans and or the sessions which are expressed by the user.
        # If the user answers all to the option, then we're obviously not going to do any repeated measures.
        #add a qpp dict so that you don't make stupid af errors again!
        qpp_dict={}
        grpby_strat = "None"
        for value in group_config_obj.qpp_grpby_strat:
            if value == 'Session':
                grpby_strat="session"
            if value == 'Scan':
                grpby_strat="scan"
            if value == 'None':
                grpby_strat="None"


        grp_by_sessions = False
        grp_by_scans = False
        grp_by_both = False
        repeated_measures=False
        if group_config_obj.qpp_sess_inclusion:
            if len(group_config_obj.qpp_sess_inclusion) > 0:
                grp_by_scans = True
        if group_config_obj.qpp_scan_inclusion:
            if len(group_config_obj.qpp_scan_inclusion) > 0:
                grp_by_sessions = True
        if grp_by_scans or grp_by_sessions:
            repeated_measures=True
        if grp_by_scans and grp_by_sessions:
            grp_by_both=True
        session_list=group_config_obj.qpp_sess_inclusion
        scan_list=group_config_obj.qpp_scan_inclusion

        #PSA #both multiple sessions and scans, youre going to do nothing
             #if neither, the output directory will not have any level of grouping, it will be ses1_scan1 --> nuisance strat, etc
        if repeated_measures:
            if grp_by_scans:
                new_output_df, dropped_parts = balance_df(output_df, session_list, scan_list)
                new_output_df = new_output_df[new_output_df["Sessions"].isin(session_list)]
                join_columns.append("Sessions")
            if grp_by_sessions:
                new_output_df, dropped_parts = balance_df(output_df, session_list, scan_list)
                #drop all the scans that are not in the scan list
                new_output_df = new_output_df[new_output_df["Scan"].isin(scan_list)]
                join_columns.append("Scan")
            if grp_by_both:
                new_output_df, dropped_parts = balance_df(output_df, session_list, scan_list)
                new_output_df = new_output_df[new_output_df["Scan"].isin(scan_list)]
                join_columns.append("Scan")
                new_output_df = new_output_df[new_output_df["Sessions"].isin(session_list)]
                join_columns.append("Sessions")

            pre_qpp_dict={}
            for scan_df_tuple in new_output_df.groupby("Scan"):
                scans = scan_df_tuple[0]
                scan_df=scan_df_tuple[1]
                    #print("scan_df:{0}".format(scan_df))
                    #if you have multiple sessions
                if grpby_strat == "session":
                    if 'Sessions' in scan_df.columns:
                        print("grouping by session")
                        for ses_df_tuple in scan_df.groupby('Sessions'):
                            session = 'ses-{0}'.format(ses_df_tuple[0])
                            ses_df = ses_df_tuple[1]
                            if not ses_df_tuple[0] in qpp_dict:
                                pre_qpp_dict[ses_df_tuple[0]] = []
                            pre_qpp_dict[ses_df_tuple[0]].append(ses_df)
                        for key in pre_qpp_dict:
                            list_df=pre_qpp_dict[key]
                            concat_df = list_df[0]
                            for df in list_df[1:]:
                                concat_df=pd.concat(concat_df,df)
                            qpp_dict[key]=concat_df
                if grpby_strat == "scan":
                    print("grouping by scan")
                    session='ses-1'
                    qpp_dict['ses-1'] = scan_df
                if grpby_strat == "None":
                    print("no grouping stratergy")
                    qpp_dict['output_df'] = new_output_df

        if repeated_measures == False:
            new_output_df, dropped_parts = balance_df(output_df, session_list,scan_list)
            qpp_dict['output_df'] = new_output_df

                #qpp_dict['output_df'] = new_output_df

        if len(qpp_dict) == 0:
            err = '\n\n[!]C-PAC says:Could not find the ' \
                      'particpants in your pipeline output directory that were ' \
                      'included in your analysis, and the particpants in the phenotype ' \
                      'phenotype file provided.\n\n'
            raise Exception(err)

    print(qpp_dict)
    return qpp_dict,inclusion_list,resource_id,strat_info

def use_inputs(group_config_file):
    import os
    import shutil
    from CPAC.pipeline.cpac_group_runner import load_config_yml
    from CPAC.pipeline.cpac_ga_model_generator import create_merged_copefile, create_merge_mask
    import nibabel as nib
    import numpy as np


    qpp_dict,inclusion_list,resource_id,strat_info = prep_inputs(group_config_file)
    group_config_obj=load_config_yml(group_config_file)
    use_other_function = False
    for key in qpp_dict:

        newer_output_df = qpp_dict[key]

        pipeline_ID = group_config_obj.pipeline_dir.rstrip('/').split('/')[-1]

        model_dir = os.path.join(group_config_obj.output_dir,
                               'cpac_group_analysis',
                               'CPAC_QPP')
        out_dir = os.path.join(model_dir,
                               resource_id,  # nuisance strat to initialize
                               strat_info,'model_files')

        participant_id = newer_output_df['participant_id'].tolist()
        for element in participant_id:
            if "Scan" in newer_output_df.columns:
                scan_list=newer_output_df["Scan"].tolist()
                nrn=len(set(scan_list))
            elif "Session" in newer_output_df.columns:
                session_list=newer_output_df["Session"].tolist()
                nrn=len(set(session_list))
            else:
                nrn=1
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        subject_list=newer_output_df["Filepath"].tolist()
        for subject in subject_list:
            sub_img=nib.load(subject)
            if sub_img.shape==3:
               use_other_function=False
               merge_outfile = os.path.join(out_dir,'_merged.nii.gz')
               merge_file = create_merged_copefile(newer_output_df["Filepath"].tolist(), merge_outfile)
               merge_mask_outfile=os.path.join(out_dir,'_merged_mask.nii.gz')
               merge_mask = create_merge_mask(merge_file, merge_mask_outfile)
                #return merge_file, merge_mask
            else:
               use_other_function=True
               merge_outfile = os.path.join(out_dir, '_merged.nii.gz')
               merge_file = create_merged_copefile(newer_output_df["Filepath"].tolist(), merge_outfile)
               merge_mask_outfile = os.path.join(out_dir, '_merged_mask.nii.gz')
               merge_mask = create_merge_mask(merge_file, merge_mask_outfile)



    return use_other_function,merge_mask,subject_list,inclusion_list,out_dir,nrn

def balance_df(new_output_df,session_list,scan_list):
    import pandas as pd
    from collections import Counter

    part_ID_count = Counter(new_output_df["participant_id"])
    print(part_ID_count)
    if scan_list and session_list:
        sessions_x_scans= len(session_list)*len(scan_list)
    if session_list:
            sessions_x_scans = len(session_list)
    if scan_list:
            sessions_x_scans = len(scan_list)
    dropped_parts = []
    for part_ID in part_ID_count.keys():
        if part_ID_count[part_ID] > 1:
            if part_ID_count[part_ID] % 2 != 0:
                new_output_df = new_output_df[new_output_df.participant_id != part_ID]
                print(new_output_df)
                dropped_parts.append(part_ID)


    return new_output_df, dropped_parts














