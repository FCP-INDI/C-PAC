"""Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>."""
import re
import os
import sys
import glob

from CPAC.utils.datasource import create_grp_analysis_dataflow
from CPAC.utils import Configuration


def write_new_sub_file(current_mod_path, subject_list, new_participant_list):
    # write the new participant list
    new_sub_file = os.path.join(current_mod_path,os.path.basename(subject_list))

    try:
        with open(new_sub_file, "w") as f:
            for part_ID in new_participant_list:
                print(part_ID, file=f)
    except Exception as e:
        err = "\n\n[!] CPAC says: Could not write new participant list for " \
              "current model and derivative in group-level analysis. Ensure "\
              "you have write access to the directory:\n%s\n\nError " \
              "details: %s\n\n" % (current_mod_path, e)
        raise Exception(err)

    return new_sub_file


def create_dir(dir_path, description):

    if not os.path.isdir(dir_path):
        try:
            os.makedirs(dir_path)
        except Exception as e:
            err = "\n\n[!] Could not create the %s directory.\n\n" \
                  "Attempted directory creation: %s\n\n" \
                  "Error details: %s\n\n" % (description, dir_path, e)
            raise Exception(err)


def create_merged_copefile(list_of_output_files, merged_outfile):

    import subprocess

    merge_string = ["fslmerge", "-t", merged_outfile]
    merge_string = merge_string + list_of_output_files

    try:
        retcode = subprocess.check_output(merge_string)
    except Exception as e:
        err = "\n\n[!] Something went wrong with FSL's fslmerge during the " \
              "creation of the 4D merged file for group analysis.\n\n" \
              "Attempted to create file: %s\n\nLength of list of files to " \
              "merge: %d\n\nError details: %s\n\n" \
              % (merged_outfile, len(list_of_output_files), e)
        raise Exception(err)

    return merged_outfile


def create_merge_mask(merged_file, mask_outfile):

    import subprocess

    mask_string = ["fslmaths", merged_file, "-abs", "-Tmin", "-bin",
                   mask_outfile]

    try:
        retcode = subprocess.check_output(mask_string)
    except Exception as e:
        err = "\n\n[!] Something went wrong with FSL's fslmaths during the " \
              "creation of the merged copefile group mask.\n\nAttempted to " \
              "create file: %s\n\nMerged file: %s\n\nError details: %s\n\n" \
              % (mask_outfile, merged_file, e)
        raise Exception(err)

    return mask_outfile


def check_merged_file(list_of_output_files, merged_outfile):

    import subprocess

    # make sure the order is correct
    #   we are ensuring each volume of the merge file correlates perfectly
    #   with the output file it should correspond to
    i = 0
    for output_file in list_of_output_files:
        test_string = ["3ddot", "-demean",
                       "{0}[{1}]".format(merged_outfile, str(i)), output_file]

        try:
            retcode = subprocess.check_output(test_string)
        except Exception as e:
            err = "\n\n[!] Something went wrong while trying to run AFNI's " \
                  "3ddot for the purpose of testing the merge file output." \
                  "\n\nError details: %s\n\n" % e
            raise Exception(err)

        retcode = retcode.rstrip("\n").rstrip("\t")

        if retcode != "1":
            err = "\n\n[!] The volumes of the merged file do not correspond "\
                  "to the correct order of output files as described in the "\
                  "phenotype matrix. If you are seeing this error, " \
                  "something possibly went wrong with FSL's fslmerge.\n\n" \
                  "Merged file: %s\n\nMismatch between merged file volume " \
                  "%d and derivative file %s\n\nEach volume should " \
                  "correspond to the derivative output file for each " \
                  "participant in the model.\n\n" \
                  % (merged_outfile, i, output_file)
            raise Exception(err)

        i += 1


def calculate_measure_mean_in_df(model_df, merge_mask):

    import subprocess
    import pandas as pd

    mm_dict_list = []

    for raw_file in model_df["Raw_Filepath"]:

        mask_string = ["3dmaskave", "-mask", merge_mask, raw_file]

        # calculate
        try:
            retcode = subprocess.check_output(mask_string)
        except Exception as e:
            err = "\n\n[!] AFNI's 3dMaskAve failed for raw output: %s\n" \
                  "Error details: %s\n\n" % (raw_file, e)
            raise Exception(err)

        # if this breaks, 3dmaskave output to STDOUT has changed
        try:
            mean = retcode.split(" ")[0]
        except Exception as e:
            err = "\n\n[!] Something went wrong with parsing the output of " \
                  "AFNI's 3dMaskAve - this is used in calculating the " \
                  "Measure Mean in group analysis.\n\nError details: %s\n\n" \
                  % e
            raise Exception(err)

        mm_dict = {}
        mm_dict["Raw_Filepath"] = raw_file
        mm_dict["Measure_Mean"] = mean
        mm_dict_list.append(mm_dict)

    mm_df = pd.DataFrame(mm_dict_list)

    # demean!
    mm_df["Measure_Mean"] = mm_df["Measure_Mean"].astype(float)
    mm_df["Measure_Mean"] = \
        mm_df["Measure_Mean"].sub(mm_df["Measure_Mean"].mean())

    model_df = pd.merge(model_df, mm_df, how="inner", on=["Raw_Filepath"])

    return model_df


def check_mask_file_resolution(data_file, roi_mask, group_mask, out_dir, \
    output_id=None):

    import os
    import subprocess
    import nibabel as nb

    # let's check if we need to resample the custom ROI mask
    raw_file_img = nb.load(data_file)
    raw_file_hdr = raw_file_img.header
    roi_mask_img = nb.load(roi_mask)
    roi_mask_hdr = roi_mask_img.header

    raw_file_dims = raw_file_hdr.get_zooms()
    roi_mask_dims = roi_mask_hdr.get_zooms()

    if raw_file_dims != roi_mask_dims:
        print("\n\nWARNING: The custom ROI mask file is a different " \
              "resolution than the output data! Resampling the ROI mask " \
              "file to match the original output data!\n\nCustom ROI mask " \
              "file: %s\n\nOutput measure: %s\n\n" % (roi_mask, output_id))

        resampled_outfile = os.path.join(out_dir, \
                                         "resampled_%s" \
                                         % os.path.basename(roi_mask))

        resample_str = ["flirt", "-in", roi_mask, "-ref", group_mask, \
                        "-applyisoxfm", str(raw_file_dims[0]), "-interp", \
                        "nearestneighbour", "-out", resampled_outfile]

        try:
            retcode = subprocess.check_output(resample_str)
        except Exception as e:
            err = "\n\n[!] Something went wrong with running FSL FLIRT for " \
                  "the purpose of resampling the custom ROI mask to match " \
                  "the original output file's resolution.\n\nCustom ROI " \
                  "mask file: %s\n\nError details: %s\n\n" % (roi_mask, e)
            raise Exception(err)

        roi_mask = resampled_outfile

    return roi_mask


def trim_mask(input_mask, ref_mask, output_mask_path):

    import os
    import subprocess

    # mask the mask
    mask_string = ["fslmaths", input_mask, "-mul", ref_mask, output_mask_path]

    try:
        retcode = subprocess.check_output(mask_string)
    except Exception as e:
        err = "\n\n[!] Something went wrong with running FSL's fslmaths for "\
              "the purpose of trimming the custom ROI masks to fit within " \
              "the merged group mask.\n\nCustom ROI mask file: %s\n\nMerged "\
              "group mask file: %s\n\nError details: %s\n\n" \
              % (input_mask, ref_mask, e)
        raise Exception(err)

    return output_mask_path


def calculate_custom_roi_mean_in_df(model_df, roi_mask):

    import os
    import subprocess
    import pandas as pd

    # calculate the ROI means
    roi_dict_list = []

    for raw_file in model_df["Raw_Filepath"]:

        roi_string = ["3dROIstats", "-mask", roi_mask, raw_file]

        try:
            retcode = subprocess.check_output(roi_string)
        except Exception as e:
            err = "\n\n[!] Something went wrong with running AFNI's " \
                  "3dROIstats while calculating the custom ROI means.\n\n" \
                  "Custom ROI mask file: %s\n\nRaw output filepath: %s\n\n" \
                  "Error details: %s\n\n" % (roi_mask, raw_file, e)
            raise Exception(err)

        # process the output string
        roi_means_string = str(retcode.split(raw_file)[1].rstrip("\n"))

        roi_means_string_list = roi_means_string.split("\t")

        # check
        roi_means_list = []
        for roi_mean in roi_means_string_list:
            try:
                roi_mean = float(roi_mean)
            except:
                continue

            roi_means_list.append(roi_mean)

        if len(roi_means_list) != retcode.count("Mean_"):
            err = "\n\n[!] Something went wrong with parsing the output "\
                  "of AFNI's 3dROIstats during the calculation of the " \
                  "custom ROI means.\n\nOutput file: %s\n\nCustom ROI " \
                  "mask file: %s\n\n3dROIstats output: %s\n\n" \
                  % (raw_file, roi_mask, retcode)
            raise Exception(err)

        # add in the custom ROI means!
        roi_dict = {}
        roi_dict["Raw_Filepath"] = raw_file

        i = 1
        for roi_mean in roi_means_list:
            roi_label = "Custom_ROI_Mean_%d" % i
            roi_dict[roi_label] = roi_mean
            i += 1

        roi_dict_list.append(roi_dict)


    roi_df = pd.DataFrame(roi_dict_list)

    # demean!
    i = 1
    for roi_mean in roi_means_list:
        roi_label = "Custom_ROI_Mean_%d" % i
        i += 1
        roi_df[roi_label] = roi_df[roi_label].astype(float)
        roi_df[roi_label] = roi_df[roi_label].sub(roi_df[roi_label].mean())

    model_df = pd.merge(model_df, roi_df, how="inner", on=["Raw_Filepath"])

    return model_df


def parse_out_covariates(design_formula):

    patsy_ops = ["~", "+", "-", "*", "/", ":", "**", ")", "("]

    for op in patsy_ops:
        if op in design_formula:
            design_formula = design_formula.replace(op, " ")

    words = design_formula.split(" ")

    covariates = [x for x in words if x != ""]

    return covariates


def split_groups(pheno_df, group_ev, ev_list, cat_list):

    import pandas as pd

    new_ev_list = []
    new_cat_list = []
    if group_ev not in cat_list:
        err = "\n\n[!] The grouping variable must be one of the categorical "\
              "covariates!\n\n"
        raise Exception(err)

    # map for the .grp file for FLAME
    idx = 1
    keymap = {}
    for val in pheno_df[group_ev]:
        if val not in keymap.keys():
            keymap[val] = idx
            idx += 1
    grp_vector = pheno_df[group_ev].map(keymap)

    # start the split
    pheno_df["subject_key"] = pheno_df["participant_id"]
    join_column = ["subject_key"]

    if "Session" in pheno_df:
        pheno_df["session_key"] = pheno_df["Session"]
        join_column.append("session_key")

    if "Series" in pheno_df:
        pheno_df["series_key"] = pheno_df["Series"]
        join_column.append("series_key")

    group_levels = list(set(pheno_df[group_ev]))

    level_df_list = []
    for level in group_levels:

        level_df = pheno_df[pheno_df[group_ev] == level]
        rename = {}

        for col in level_df.columns:
            if (col != group_ev) and (col not in join_column) and (col in ev_list):
                rename[col] = col + "__FOR_%s_%s" % (group_ev, level)
                if rename[col] not in new_ev_list:
                    new_ev_list.append(rename[col])
                if (col in cat_list) and (rename[col] not in new_cat_list):
                    new_cat_list.append(rename[col])

        for other_lev in group_levels:
            if other_lev != level:
                for col in level_df.columns:
                    if (col != group_ev) and (col not in join_column) and (col in ev_list):
                        newcol = col + "__FOR_%s_%s" % (group_ev, other_lev)
                        level_df[newcol] = 0
                        if newcol not in new_ev_list:
                            new_ev_list.append(newcol)
                        if col in cat_list:
                            if newcol not in new_cat_list:
                                new_cat_list.append(newcol)
        level_df.rename(columns=rename, inplace=True)
        level_df_list.append(level_df)

        # the grouping variable has to be in the categorical list too
        #new_cat_list.append(group_ev)

        # get it back into order!
    pheno_df = pheno_df[join_column].merge(pd.concat(level_df_list), on=join_column)
    pheno_df = pheno_df.drop(join_column,1)

    return pheno_df, grp_vector, new_ev_list, new_cat_list


def patsify_design_formula(formula, categorical_list, encoding="Treatment"):

    closer = ")"
    if encoding == "Treatment":
        closer = ")"
    elif encoding == "Sum":
        closer = ", Sum)"

    # pad with spaces if they aren't present
    formula = formula.replace("+", " + ")
    formula = formula.replace("-", " - ")
    formula = formula.replace("=", " = ")
    formula = formula.replace("(", " ( ").replace(")", " ) ")
    formula = formula.replace("*", " * ")
    formula = formula.replace("/", " / ")

    # pad end with single space so the formula.replace below won't miss the last
    # covariate when relevant
    formula = '{0} '.format(formula)

    for ev in categorical_list:
        if ev in formula:
            new_ev = "C(" + ev + closer
            formula = formula.replace(" {0} ".format(ev), new_ev)

    # remove Intercept - if user wants one, they should add "+ Intercept" when
    # specifying the design formula
    if ("Intercept" not in formula) and ("intercept" not in formula):
        formula = formula + " - 1"
    else:
        # having "Intercept" in the formula is really just a flag to prevent
        # "- 1" from being added to the formula - we don't actually want
        # "Intercept" in the formula
        formula = formula.replace("+ Intercept", "")
        formula = formula.replace("+Intercept", "")
        formula = formula.replace("+ intercept", "")
        formula = formula.replace("+intercept", "")

    return formula


def check_multicollinearity(matrix):

    import numpy as np

    print("\nChecking for multicollinearity in the model..")

    U, s, V = np.linalg.svd(matrix)

    max_singular = np.max(s)
    min_singular = np.min(s)

    print("Max singular: ", max_singular)
    print("Min singular: ", min_singular)
    print("Rank: ", np.linalg.matrix_rank(matrix), "\n")

    if min_singular == 0:
        print('[!] CPAC warns: Detected multicollinearity in the ' \
                  'computed group-level analysis model. Please double-' \
                  'check your model design.\n\n')
    else:
        condition_number = float(max_singular)/float(min_singular)
        print("Condition number: %f" % condition_number)
        if condition_number > 30:
            print('[!] CPAC warns: Detected multicollinearity in the ' \
                      'computed group-level analysis model. Please double-' \
                      'check your model design.\n\n')
        else:
            print('Looks good..\n')


def create_contrasts_dict(dmatrix_obj, contrasts_list, output_measure):

    contrasts_vectors = []

    for con_equation in contrasts_list:

        try:
            lincon = dmatrix_obj.design_info.linear_constraint(str(con_equation))
        except Exception as e:
            err = "\n\n[!] Could not process contrast equation:\n%s\n\n" \
                  "Design matrix EVs/covariates:\n%s\n\nError details:\n%s" \
                  "\n\nNote: If the design matrix EVs are different than " \
                  "what was shown in the model design creator, this may be " \
                  "because missing participants, sessions, or series for " \
                  "this measure (%s) may have altered the group design.\n\n" \
                  % (str(con_equation), dmatrix_obj.design_info.column_names,\
                     e, output_measure)
            raise Exception(err)

        con_vec = lincon.coefs[0]
        contrasts_vectors.append(con_vec)


    return contrasts_vectors


def build_feat_model(model_df, model_name, group_config_file, resource_id,
                     preproc_strat, session_id, series_or_repeated_label):

    #
    # this function runs once per derivative type and preproc strat combo
    # during group analysis
    #

    import os
    import patsy
    import pandas as pd
    import numpy as np

    from CPAC.pipeline import nipype_pipeline_engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.io as nio
    from CPAC.pipeline.cpac_group_runner import load_config_yml

    from CPAC.utils.create_group_analysis_info_files import write_design_matrix_csv, \
        write_blank_contrast_csv

    group_config_obj = load_config_yml(group_config_file)

    pipeline_ID = group_config_obj.pipeline_dir.rstrip('/').split('/')[-1]
    #sublist_txt = group_config_obj.participant_list

    #if sublist_txt == None:
    #    print ("Warning! You have not provided a subject list. CPAC will use all the subjects in pipeline directory")
    #    sublist_txt = group_config_obj.participant_list
    #else:
    #    sublist_txt = group_config_obj.particpant_list

    # remove file names from preproc_strat
    filename = preproc_strat.split("/")[-1]
    preproc_strat = preproc_strat.replace('.nii', '').replace('.gz', '')
    preproc_strat = preproc_strat.lstrip("/").rstrip("/")

    ftest_list = []
    readme_flags = []

    # determine if f-tests are included or not
    custom_confile = group_config_obj.custom_contrasts

    if ((custom_confile is None) or (custom_confile == '') or
            ("None" in custom_confile) or ("none" in custom_confile)):
        custom_confile = None

    #    if (len(group_config_obj.f_tests) == 0) or \
    #            (group_config_obj.f_tests is None):
    #        fTest = False
    #    else:
    #        fTest = True
    #        ftest_list = group_config_obj.f_tests

    #else:
    #    if not os.path.exists(custom_confile):
    #        errmsg = "\n[!] CPAC says: You've specified a custom contrasts " \
    #                 ".CSV file for your group model, but this file cannot " \
    #                 "be found. Please double-check the filepath you have " \
    #                 "entered.\n\nFilepath: %s\n\n" % custom_confile
    #        raise Exception(errmsg)#

    #    with open(custom_confile, "r") as f:
    #        evs = f.readline()

    #    evs = evs.rstrip('\r\n').split(',')
    #    count_ftests = 0

    #    fTest = False

    #    for ev in evs:
    #        if "f_test" in ev:
    #            count_ftests += 1

    # create path for output directory
    model_dir = os.path.join(group_config_obj.output_dir,
                             'cpac_group_analysis',
                             'FSL_FEAT',
                             '{0}'.format(pipeline_ID),
                             'group_model_{0}'.format(model_name))

    out_dir = os.path.join(model_dir,
                           resource_id,
                           session_id,
                           series_or_repeated_label,
                           preproc_strat)

    try:
        preset_contrast = group_config_obj.preset
        preset = True
    except AttributeError:
        preset = False

    if 'sca_roi' in resource_id:
        out_dir = os.path.join(out_dir,
            re.search('sca_ROI_(\d)+', os.path.splitext(\
                os.path.splitext(os.path.basename(\
                    model_df["Filepath"][0]))[0])[0]).group(0))

    if 'dr_tempreg_maps_zstat_files_to_standard_smooth' in resource_id:
        out_dir = os.path.join(out_dir,
            re.search('temp_reg_map_z_(\d)+', os.path.splitext(\
                os.path.splitext(os.path.basename(\
                    model_df["Filepath"][0]))[0])[0]).group(0))

    if 'centrality' in resource_id:
        names = ['degree_centrality_binarize',
                 'degree_centrality_weighted',
                 'eigenvector_centrality_binarize',
                 'eigenvector_centrality_weighted',
                 'lfcd_binarize', 'lfcd_weighted']

        for name in names:
            if name in filename:
                out_dir = os.path.join(out_dir, name)
                break

    if 'tempreg_maps' in resource_id:
        out_dir = os.path.join(out_dir, re.search('\w*[#]*\d+',
            os.path.splitext(os.path.splitext(os.path.basename(\
                model_df["Filepath"][0]))[0])[0]).group(0))

    model_path = os.path.join(out_dir, 'model_files')

    # create the actual directories
    create_dir(model_path, "group analysis output")

    # create new subject list based on which subjects are left after checking
    # for missing outputs
    new_participant_list = []
    for part in model_df["participant_id"]:
        # do this instead of using "set" just in case, to preserve order
        #   only reason there may be duplicates is because of multiple-series
        #   repeated measures runs
        if part not in new_participant_list:
            new_participant_list.append(part)

    if group_config_obj.participant_list == None:
        #participant_list = os.listdir(group_config_obj.pipeline_dir)
        new_sub_file = write_new_sub_file(model_path,
                                          group_config_obj.pipeline_dir,
                                          new_participant_list)
    else:
        new_sub_file = write_new_sub_file(model_path,
                                      group_config_obj.participant_list,
                                      new_participant_list)

    group_config_obj.update('participant_list', new_sub_file)

    num_subjects = len(list(model_df["participant_id"]))

    # start processing the dataframe further
    design_formula = group_config_obj.design_formula

    # demean EVs set for demeaning
    for demean_EV in group_config_obj.ev_selections.get("demean",[]):
        model_df[demean_EV] = model_df[demean_EV].astype(float)
        model_df[demean_EV] = model_df[demean_EV].sub(model_df[demean_EV].mean())

    # demean the motion params
    if ("MeanFD" in design_formula) or ("MeanDVARS" in design_formula):
        params = ["MeanFD_Power", "MeanFD_Jenkinson", "MeanDVARS"]
        for param in params:
            model_df[param] = model_df[param].astype(float)
            model_df[param] = model_df[param].sub(model_df[param].mean())

    # create 4D merged copefile, in the correct order, identical to design
    # matrix
    merge_outfile = model_name + "_" + resource_id + "_merged.nii.gz"
    merge_outfile = os.path.join(model_path, merge_outfile)
    merge_file = create_merged_copefile(model_df["Filepath"].tolist(),
                                        merge_outfile)

    # create merged group mask
    merge_mask_outfile = '_'.join([model_name, resource_id,
                                   "merged_mask.nii.gz"])
    merge_mask_outfile = os.path.join(model_path, merge_mask_outfile)
    merge_mask = create_merge_mask(merge_file, merge_mask_outfile)

    if "Group Mask" in group_config_obj.mean_mask:
        mask_for_means = merge_mask
    else:
        individual_masks_dir = os.path.join(model_path,
                                            "individual_masks")
        create_dir(individual_masks_dir, "individual masks")
        for unique_id, series_id, raw_filepath in zip(
                model_df["participant_id"],
                model_df["Series"], model_df["Raw_Filepath"]):
            mask_for_means_path = os.path.join(individual_masks_dir,
                                               "%s_%s_%s_mask.nii.gz" % (
                                               unique_id, series_id,
                                               resource_id))
            mask_for_means = create_merge_mask(raw_filepath,
                                               mask_for_means_path)
        readme_flags.append("individual_masks")

    # calculate measure means, and demean
    if "Measure_Mean" in design_formula:
        model_df = calculate_measure_mean_in_df(model_df, mask_for_means)

    # calculate custom ROIs, and demean (in workflow?)
    if "Custom_ROI_Mean" in design_formula:

        custom_roi_mask = group_config_obj.custom_roi_mask

        if (custom_roi_mask == None) or (custom_roi_mask == "None") or \
                (custom_roi_mask == "none") or (custom_roi_mask == ""):
            err = "\n\n[!] You included 'Custom_ROI_Mean' in your design " \
                  "formula, but you didn't supply a custom ROI mask file." \
                  "\n\nDesign formula: %s\n\n" % design_formula
            raise Exception(err)

        # make sure the custom ROI mask file is the same resolution as the
        # output files - if not, resample and warn the user
        roi_mask = check_mask_file_resolution(list(model_df["Raw_Filepath"])[0],
                                              custom_roi_mask, mask_for_means,
                                              model_path, resource_id)


        # trim the custom ROI mask to be within mask constraints
        output_mask = os.path.join(model_path, "masked_%s" \
                                   % os.path.basename(roi_mask))
        roi_mask = trim_mask(roi_mask, mask_for_means, output_mask)
        readme_flags.append("custom_roi_mask_trimmed")

        # calculate
        model_df = calculate_custom_roi_mean_in_df(model_df, roi_mask)

        # update the design formula
        new_design_substring = ""

        for col in model_df.columns:
            if "Custom_ROI_Mean_" in str(col):
                if str(col) == "Custom_ROI_Mean_1":
                    new_design_substring = new_design_substring + " %s" % col
                else:
                    new_design_substring = new_design_substring +" + %s" % col
        design_formula = design_formula.replace("Custom_ROI_Mean",
                                                new_design_substring)

    cat_list = []
    if "categorical" in group_config_obj.ev_selections.keys():
        cat_list = group_config_obj.ev_selections["categorical"]

    # prep design for repeated measures, if applicable
    if len(group_config_obj.sessions_list) > 0:
        if "session" in model_df.columns:
            # if these columns were added by the model builder automatically
            design_formula = design_formula + " + session"
            if "session" not in cat_list:
                cat_list.append("session")

    if len(group_config_obj.series_list) > 0:
        design_formula = design_formula + " + Series"
        if "Series" not in cat_list:
            cat_list.append("Series")

    if "session" in model_df.columns:
        # if these columns were added by the model builder automatically
        for col in model_df.columns:
            # should only grab the repeated measures-designed participant_{ID}
            # columns, not the "participant_id" column!
            if "participant_" in col and "_id" not in col:
                design_formula = design_formula + " + %s" % col
                cat_list.append(col)

    # parse out the EVs in the design formula at this point in time
    #   this is essentially a list of the EVs that are to be included
    ev_list = parse_out_covariates(design_formula)

    # SPLIT GROUPS here.
    #   CURRENT PROBLEMS: was creating a few doubled-up new columns
    grp_vector = [1] * num_subjects

    if group_config_obj.group_sep:

        # check if the group_ev parameter is a list instead of a string:
        # this was added to handle the new group-level analysis presets. this
        # is the only modification that was required to the group analysis
        # workflow, and it handles cases where the group variances must be
        # modeled separately, by creating separate groups for the FSL FLAME
        # .grp file.
        #     the group_ev parameter gets sent in as a list if coming from any
        #     of the presets that deal with multiple groups- in these cases,
        #     the pheno_df/design matrix is already set up properly for the
        #     multiple groups, and we need to bypass all of the processing
        #     that usually occurs when the "modeling group variances
        #     separately" option is enabled in the group analysis config YAML
        group_ev = group_config_obj.grouping_var

        if isinstance(group_ev, list) or "," in group_ev:
            grp_vector = []

            if "," in group_ev:
                group_ev = group_ev.split(",")

            if len(group_ev) == 2:
                for x, y in zip(model_df[group_ev[0]], model_df[group_ev[1]]):
                    if x == 1:
                        grp_vector.append(1)
                    elif y == 1:
                        grp_vector.append(2)
                    else:
                        err = "\n\n[!] The two categorical covariates you " \
                              "provided as the two separate groups (in order " \
                              "to model each group's variances separately) " \
                              "either have more than 2 levels (1/0), or are " \
                              "not encoded as 1's and 0's.\n\nCovariates:\n" \
                              "{0}\n{1}\n\n".format(group_ev[0], group_ev[1])
                        raise Exception(err)

            elif len(group_ev) == 3:
                for x, y, z in zip(model_df[group_ev[0]], model_df[group_ev[1]],
                                   model_df[group_ev[2]]):
                    if x == 1:
                        grp_vector.append(1)
                    elif y == 1:
                        grp_vector.append(2)
                    elif z == 1:
                        grp_vector.append(3)
                    else:
                        err = "\n\n[!] The three categorical covariates you " \
                              "provided as the three separate groups (in order " \
                              "to model each group's variances separately) " \
                              "either have more than 2 levels (1/0), or are " \
                              "not encoded as 1's and 0's.\n\nCovariates:\n" \
                              "{0}\n{1}\n{2}\n\n".format(group_ev[0],
                                                         group_ev[1],
                                                         group_ev[2])
                        raise Exception(err)

            else:
                # we're only going to see this if someone plays around with
                # their preset or config file manually
                err = "\n\n[!] If you are seeing this message, it's because:\n" \
                      "1. You are using the group-level analysis presets\n" \
                      "2. You are running a model with multiple groups having " \
                      "their variances modeled separately (i.e. multiple " \
                      "values in the FSL FLAME .grp input file), and\n" \
                      "3. For some reason, the configuration has been set up " \
                      "in a way where CPAC currently thinks you're including " \
                      "only one group, or more than three, neither of which " \
                      "are supported.\n\nGroups provided:\n{0}" \
                      "\n\n".format(str(group_ev))
                raise Exception(err)

        else:
            # model group variances separately
            old_ev_list = ev_list

            model_df, grp_vector, ev_list, cat_list = split_groups(model_df,
                                    group_config_obj.grouping_var,
                                    ev_list, cat_list)

            # make the grouping variable categorical for Patsy (if we try to
            # do this automatically below, it will categorical-ize all of
            # the substrings too)
            design_formula = design_formula.replace(group_config_obj.grouping_var,
                                      "C(" + group_config_obj.grouping_var + ")")
            if group_config_obj.coding_scheme == "Sum":
                design_formula = design_formula.replace(")", ", Sum)")

            # update design formula
            rename = {}
            for old_ev in old_ev_list:
                for new_ev in ev_list:
                    if old_ev + "__FOR" in new_ev:
                        if old_ev not in rename.keys():
                            rename[old_ev] = []
                        rename[old_ev].append(new_ev)

            for old_ev in rename.keys():
                design_formula = design_formula.replace(old_ev,
                                                        " + ".join(rename[old_ev]))

    # prep design formula for Patsy
    design_formula = patsify_design_formula(design_formula, cat_list,
                                            group_config_obj.coding_scheme[0])

    if not preset:
        # send to Patsy
        try:
            dmatrix = patsy.dmatrix(design_formula, model_df)
            dmatrix.design_info.column_names.append(model_df["Filepath"])
            dmatrix_column_names = dmatrix.design_info.column_names
        except Exception as e:
            err = "\n\n[!] Something went wrong with processing the group model "\
                  "design matrix using the Python Patsy package. Patsy might " \
                  "not be properly installed, or there may be an issue with the "\
                  "formatting of the design matrix.\n\nDesign matrix columns: " \
                  "%s\n\nPatsy-formatted design formula: %s\n\nError details: " \
                  "%s\n\n" % (model_df.columns, design_formula, e)
            raise Exception(err)
    else:
        if 'Sessions' in model_df:
            sess_levels = list(set(list(model_df['Sessions'].values)))
            if len(sess_levels) > 1:
                sess_map = {sess_levels[0]: '1', sess_levels[1]: '-1'}
                if len(sess_levels) == 3:
                    sess_map.update({sess_levels[2]: '0'})
                new_sess = [s.replace(s, sess_map[s]) for s in list(model_df['Sessions'].values)]
                model_df['Sessions'] = new_sess
        if 'Series' in model_df:
            sess_levels = list(set(list(model_df['Series'].values)))
            if len(sess_levels) > 1:
                sess_map = {sess_levels[0]: '1', sess_levels[1]: '-1'}
                if len(sess_levels) == 3:
                    sess_map.update({sess_levels[2]: '0'})
                new_sess = [s.replace(s, sess_map[s]) for s in list(model_df['Series'].values)]
                model_df['Series'] = new_sess

        keep_cols = [x for x in model_df.columns if x in design_formula]
        dmatrix = model_df[keep_cols].astype('float')
        dmatrix_column_names = list(dmatrix.columns)

    # check the model for multicollinearity - Patsy takes care of this, but
    # just in case
    check_multicollinearity(np.array(dmatrix))

    dmat_csv_path = os.path.join(model_path, "design_matrix.csv")
    contrast_out_path = os.path.join(out_dir, "contrast.csv")

    # make sure "column_names" is in the same order as the original EV column
    # header ordering in model_df - mainly for repeated measures, to make sure
    # participants_<ID> cols are at end for clarity for users
    dmat_cols = []
    dmat_id_cols = []
    for dmat_col in dmatrix_column_names:
        if 'participant_' in dmat_col:
            dmat_id_cols.append(dmat_col)
        else:
            dmat_cols.append(dmat_col)
    column_names = dmat_cols
    dmat_id_cols = sorted(dmat_id_cols)
    column_names += dmat_id_cols

    # check to make sure there are more time points than EVs!
    if len(column_names) >= num_subjects:
        err = "\n\n################## MODEL NOT GENERATED ##################" \
              "\n\n[!] CPAC says: There are more EVs than there are " \
              "participants currently included in the model for:\n\n" \
              "Derivative: {0}\nSession: {1}\nScan: {2}\nPreproc strategy:" \
              "\n    {3}\n\n" \
              "There must be more participants than EVs in the design.\n\n" \
              "Number of participants: {4}\nNumber of EVs: {5}\n\nEV/" \
              "covariate list: {6}\n\nNote: If you specified to model group " \
              "variances separately, the amount of EVs can nearly double " \
              "once they are split along the grouping variable.\n\nIf the " \
              "number of participants is lower than the number of " \
              "participants in your group analysis inclusion list, this " \
              "may be because not every participant originally included has " \
              "an output for {7} for this scan and preprocessing strategy in " \
              "the individual-level analysis output directory.\n\nDesign " \
              "formula going in: {8}" \
              "\n\n#########################################################" \
              "\n\n".format(resource_id, session_id, series_or_repeated_label,
                            preproc_strat, num_subjects, len(column_names),
                            column_names, resource_id, design_formula)
        print(err)

    # check the merged file's order
    check_merged_file(model_df["Filepath"], merge_file)

    # we must demean the categorical regressors if the Intercept/Grand Mean
    # is included in the model, otherwise FLAME produces blank outputs
    if "Intercept" in column_names:
        cat_indices = []
        col_name_indices = dmatrix.design_info.column_name_indexes
        for col_name in col_name_indices.keys():
            if "C(" in col_name:
                cat_indices.append(int(col_name_indices[col_name]))

        # note: dmat_T is now no longer a DesignMatrix Patsy object, but only
        # an array
        dmat_T = dmatrix.transpose()

        for index in cat_indices:
            new_row = []
            for val in dmat_T[index]:
                new_row.append(val - dmat_T[index].mean())
            dmat_T[index] = new_row

        # we can go back, but we won't be the same
        dmatrix = dmat_T.transpose()
        readme_flags.append("cat_demeaned")

    dmatrix_df = pd.DataFrame(np.array(dmatrix),
                              index=model_df["participant_id"],
                              columns=dmatrix_column_names)
    cols = dmatrix_df.columns.tolist()

    # make sure "column_names" is in the same order as the original EV column
    # header ordering in model_df - mainly for repeated measures, to make sure
    # participants_<ID> cols are at end for clarity for users
    dmat_cols = []
    dmat_id_cols = []
    for dmat_col in cols:
        if 'participant_' in dmat_col:
            dmat_id_cols.append(dmat_col)
        else:
            dmat_cols.append(dmat_col)
    column_names = dmat_cols
    dmat_id_cols = sorted(dmat_id_cols)
    column_names += dmat_id_cols

    dmatrix_df = dmatrix_df[column_names]

    dmat_csv_path = os.path.join(model_path, "design_matrix.csv")

    write_design_matrix_csv(dmatrix_df, model_df["participant_id"],
                            column_names, dmat_csv_path)

    # time for contrasts
    if (group_config_obj.custom_contrasts == None) or (group_config_obj.contrasts == None):
        # if no custom contrasts matrix CSV provided (i.e. the user
        # specified contrasts in the GUI)
        contrasts_columns = column_names
        if group_config_obj.f_tests:
            for i in group_config_obj.f_tests[1:len(group_config_obj.f_tests)-1]:
                contrasts_columns.append('f_test_{0}'.format(i))
    else:
        pass

    contrast_out_path = os.path.join(model_dir, "contrasts.csv")

    if preset:
        cons = pd.read_csv(group_config_obj.custom_contrasts)
        with open(contrast_out_path, "w") as f:
            cons.to_csv(f, index=False)
    else:
        if os.path.isfile(contrast_out_path):
            contrasts_df = pd.read_csv(contrast_out_path)
            if contrasts_df.shape[0] > 1 or np.count_nonzero(contrasts_df.values[0][1:]) > 0:
                msg = "\n\n[!] C-PAC says: It appears you have modified your " \
                      "contrasts CSV file already- back up this file before " \
                      "building your model again to avoid overwriting your " \
                      "changes.\n\nContrasts file:\n{0}" \
                      "\n\n".format(contrast_out_path)
                raise Exception(msg)

        with open(contrast_out_path, "w") as f:
            f.write('Contrasts')
            for col in contrasts_columns:
                f.write(',{0}'.format(col))
            f.write('\ncontrast_1')
            for col in contrasts_columns:
                f.write(',0')

    groups_out_path = os.path.join(model_path, 'groups.txt')
    with open(groups_out_path, 'w') as f:
        for val in grp_vector:
            f.write('{0}\n'.format(val))

    msg = 'Model successfully generated for..\nDerivative: {0}\nSession: {1}' \
          '\nScan: {2}\nPreprocessing strategy:\n    {3}\n\nModel directory:' \
          '\n{4}\n\nGroup configuration file:\n{5}\n\nContrasts template CSV:' \
          '\n{6}\n\nDefine your contrasts in this contrasts template CSV and ' \
          'save your changes, then run FSL-FEAT ' \
          'through the command-line like so:\n\n    cpac group ' \
          'feat run <path to group config.yml>' \
           '\n'.format(resource_id, session_id, series_or_repeated_label,
                       preproc_strat, model_path, group_config_file,
                       contrast_out_path)
    print('-------------------------------------------------------------------')
    print(msg)
    print('-------------------------------------------------------------------')

    return dmat_csv_path, new_sub_file, contrast_out_path
