import re
import os
import sys
import glob

from CPAC.utils.datasource import create_grp_analysis_dataflow
from CPAC.utils import Configuration
from CPAC.utils.utils import prepare_gp_links
from CPAC.group_analysis import create_group_analysis


def write_new_sub_file(current_mod_path, subject_list, new_participant_list):

    # write the new participant list
    new_sub_file = os.path.join(current_mod_path, \
                                    os.path.basename(subject_list))

    try:
        with open(new_sub_file, "w") as f:
            for part_ID in new_participant_list:
                print >>f, part_ID
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

    mask_string = ["fslmaths", merged_file, "-abs", "-Tmin", "-bin", \
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
        test_string = ["3ddot", "-demean",  merged_outfile + "[" + str(i) + "]", \
           output_file]

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
    
    for raw_file in list(model_df["Raw_Filepath"]):
    
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
    raw_file_hdr = raw_file_img.get_header()
    roi_mask_img = nb.load(roi_mask)
    roi_mask_hdr = roi_mask_img.get_header()

    raw_file_dims = raw_file_hdr.get_zooms()
    roi_mask_dims = roi_mask_hdr.get_zooms()

    if raw_file_dims != roi_mask_dims:
        print "\n\nWARNING: The custom ROI mask file is a different " \
              "resolution than the output data! Resampling the ROI mask " \
              "file to match the original output data!\n\nCustom ROI mask " \
              "file: %s\n\nOutput measure: %s\n\n" % (roi_mask, output_id)

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
    
    for raw_file in list(model_df["Raw_Filepath"]):
        
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

    patsy_ops = ["~","+","-","*","/",":","**",")","("]

    for op in patsy_ops:
        if op in design_formula:
            design_formula = design_formula.replace(op," ")

    words = design_formula.split(" ")

    covariates = [x for x in words if x != ""]

    return covariates



def split_groups(pheno_df, group_ev, ev_list, cat_list):   
            
    import pandas as pd
        
    new_ev_list = []
    new_cat_list = []
    print group_ev
    print cat_list
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
    grp_vector = list(pheno_df[group_ev].map(keymap))
            
    # start the split
    pheno_df["subject_key"] = pheno_df["Participant"]
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

    for ev in categorical_list:
        if ev in formula:
            new_ev = "C(" + ev + closer
            formula = formula.replace(ev, new_ev)

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

    print "\nChecking for multicollinearity in the model.."

    U, s, V = np.linalg.svd(matrix)

    max_singular = np.max(s)
    min_singular = np.min(s)

    print "Max singular: ", max_singular
    print "Min singular: ", min_singular
    print "Rank: ", np.linalg.matrix_rank(matrix), "\n"

    if min_singular == 0:

        print '[!] CPAC warns: Detected multicollinearity in the ' \
                  'computed group-level analysis model. Please double-' \
                  'check your model design.\n\n'

    else:

        condition_number = float(max_singular)/float(min_singular)
        print "Condition number: %f\n\n" % condition_number
        if condition_number > 30:

            print '[!] CPAC warns: Detected multicollinearity in the ' \
                      'computed group-level analysis model. Please double-' \
                      'check your model design.\n\n'



def create_contrasts_dict(dmatrix_obj, contrasts_list, output_measure):

    contrasts_dict = {}

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
        contrasts_dict[con_equation] = con_vec

    return contrasts_dict



def prep_group_analysis_workflow(model_df, pipeline_config_path, \
    model_name, group_config_path, resource_id, preproc_strat, \
    series_or_repeated_label):
    
    #
    # this function runs once per derivative type and preproc strat combo
    # during group analysis
    #

    import os
    import patsy
    import numpy as np

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.io as nio

    from CPAC.pipeline.cpac_group_runner import load_config_yml
    from CPAC.utils.create_flame_model_files import create_flame_model_files
    from CPAC.utils.create_group_analysis_info_files import write_design_matrix_csv

    pipeline_config_obj = load_config_yml(pipeline_config_path)
    group_config_obj = load_config_yml(group_config_path)

    pipeline_ID = pipeline_config_obj.pipelineName

    # remove file names from preproc_strat
    filename = preproc_strat.split("/")[-1]
    preproc_strat = preproc_strat.replace(filename,"")
    preproc_strat = preproc_strat.lstrip("/").rstrip("/")

    # get thresholds
    z_threshold = float(group_config_obj.z_threshold[0])

    p_threshold = float(group_config_obj.p_threshold[0])

    sub_id_label = group_config_obj.participant_id_label

    ftest_list = []
    readme_flags = []

    # determine if f-tests are included or not
    custom_confile = group_config_obj.custom_contrasts

    if ((custom_confile == None) or (custom_confile == '') or \
            ("None" in custom_confile) or ("none" in custom_confile)):

        custom_confile = None

        if (len(group_config_obj.f_tests) == 0) or \
            (group_config_obj.f_tests == None):
            fTest = False
        else:
            fTest = True
            ftest_list = group_config_obj.f_tests

    else:

        if not os.path.exists(custom_confile):
            errmsg = "\n[!] CPAC says: You've specified a custom contrasts " \
                     ".CSV file for your group model, but this file cannot " \
                     "be found. Please double-check the filepath you have " \
                     "entered.\n\nFilepath: %s\n\n" % custom_confile
            raise Exception(errmsg)

        with open(custom_confile,"r") as f:
            evs = f.readline()

        evs = evs.rstrip('\r\n').split(',')
        count_ftests = 0

        fTest = False

        for ev in evs:
            if "f_test" in ev:
                count_ftests += 1

        if count_ftests > 0:
            fTest = True


    # create path for output directory
    out_dir = os.path.join(group_config_obj.output_dir, \
        "group_analysis_results_%s" % pipeline_ID, \
        "group_model_%s" % model_name, resource_id, \
        series_or_repeated_label, preproc_strat)

    if 'sca_roi' in resource_id:
        out_dir = os.path.join(out_dir, \
            re.search('sca_roi_(\d)+',os.path.splitext(\
                os.path.splitext(os.path.basename(\
                    model_df["Filepath"][0]))[0])[0]).group(0))
            
    if 'dr_tempreg_maps_zstat_files_to_standard_smooth' in resource_id:
        out_dir = os.path.join(out_dir, \
            re.search('temp_reg_map_z_(\d)+',os.path.splitext(\
                os.path.splitext(os.path.basename(\
                    model_df["Filepath"][0]))[0])[0]).group(0))
            
    if 'centrality' in resource_id:
        names = ['degree_centrality_binarize', 'degree_centrality_weighted', \
                 'eigenvector_centrality_binarize', \
                 'eigenvector_centrality_weighted', \
                 'lfcd_binarize', 'lfcd_weighted']

        for name in names:
            if name in filename:
                out_dir = os.path.join(out_dir, name)
                break

    if 'tempreg_maps' in resource_id:
        out_dir = os.path.join(out_dir, re.search('\w*[#]*\d+', \
            os.path.splitext(os.path.splitext(os.path.basename(\
                model_df["Filepath"][0]))[0])[0]).group(0))

    model_path = os.path.join(out_dir, 'model_files')

    second_half_out = \
        out_dir.split("group_analysis_results_%s" % pipeline_ID)[1]

    # generate working directory for this output's group analysis run
    work_dir = os.path.join(pipeline_config_obj.workingDirectory, \
        "group_analysis", second_half_out.lstrip("/"))

    log_dir = os.path.join(pipeline_config_obj.logDirectory, \
        "group_analysis", second_half_out.lstrip("/"))

    # create the actual directories
    create_dir(model_path, "group analysis output")
    create_dir(work_dir, "group analysis working")
    create_dir(log_dir, "group analysis logfile")


    # create new subject list based on which subjects are left after checking
    # for missing outputs
    new_participant_list = []
    for part in list(model_df["Participant"]):
        # do this instead of using "set" just in case, to preserve order
        #   only reason there may be duplicates is because of multiple-series
        #   repeated measures runs
        if part not in new_participant_list:
            new_participant_list.append(part)

    new_sub_file = write_new_sub_file(model_path, \
                                      group_config_obj.participant_list, \
                                      new_participant_list)

    group_config_obj.update('participant_list',new_sub_file)

    num_subjects = len(list(model_df["Participant"]))


    # start processing the dataframe further
    design_formula = group_config_obj.design_formula

    # demean EVs set for demeaning
    for demean_EV in group_config_obj.ev_selections["demean"]:
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

    merge_file = create_merged_copefile(list(model_df["Filepath"]), \
                                        merge_outfile)

    # create merged group mask
    merge_mask_outfile = model_name + "_" + resource_id + \
                             "_merged_mask.nii.gz"
    merge_mask_outfile = os.path.join(model_path, merge_mask_outfile)
    merge_mask = create_merge_mask(merge_file, merge_mask_outfile)

    if "Group Mask" in group_config_obj.mean_mask:
        mask_for_means = merge_mask
    else:
        individual_masks_dir = os.path.join(model_path, "individual_masks")
        create_dir(individual_masks_dir, "individual masks")
        for unique_id, series_id, raw_filepath in zip(model_df["Participant"],
            model_df["Series"], model_df["Raw_Filepath"]):
            
            mask_for_means_path = os.path.join(individual_masks_dir,
                "%s_%s_%s_mask.nii.gz" % (unique_id, series_id, resource_id))
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
        roi_mask = check_mask_file_resolution(list(model_df["Raw_Filepath"])[0], \
                                              custom_roi_mask, mask_for_means, \
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
        design_formula = design_formula.replace("Custom_ROI_Mean", \
                                                new_design_substring)


    cat_list = []
    if "categorical" in group_config_obj.ev_selections.keys():
        cat_list = group_config_obj.ev_selections["categorical"]


    # prep design for repeated measures, if applicable
    if len(group_config_obj.sessions_list) > 0:
        design_formula = design_formula + " + Session"
        if "Session" not in cat_list:
            cat_list.append("Session")
    if len(group_config_obj.series_list) > 0:
        design_formula = design_formula + " + Series"
        if "Series" not in cat_list:
            cat_list.append("Series")
    for col in list(model_df.columns):
        if "participant_" in col:
            design_formula = design_formula + " + %s" % col
            cat_list.append(col)


    # parse out the EVs in the design formula at this point in time
    #   this is essentially a list of the EVs that are to be included
    ev_list = parse_out_covariates(design_formula)


    # SPLIT GROUPS here.
    #   CURRENT PROBLEMS: was creating a few doubled-up new columns
    grp_vector = [1] * num_subjects

    if group_config_obj.group_sep:

        # model group variances separately
        old_ev_list = ev_list

        model_df, grp_vector, ev_list, cat_list = split_groups(model_df, \
                                group_config_obj.grouping_var, \
                                ev_list, cat_list)

        # make the grouping variable categorical for Patsy (if we try to
        # do this automatically below, it will categorical-ize all of 
        # the substrings too)
        design_formula = design_formula.replace(group_config_obj.grouping_var, \
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
            design_formula = design_formula.replace(old_ev, \
                                                   " + ".join(rename[old_ev]))


    # prep design formula for Patsy
    design_formula = patsify_design_formula(design_formula, cat_list, \
                         group_config_obj.coding_scheme[0])
    print design_formula
    # send to Patsy
    try:
        dmatrix = patsy.dmatrix(design_formula, model_df)
    except Exception as e:
        err = "\n\n[!] Something went wrong with processing the group model "\
              "design matrix using the Python Patsy package. Patsy might " \
              "not be properly installed, or there may be an issue with the "\
              "formatting of the design matrix.\n\nPatsy-formatted design " \
              "formula: %s\n\nError details: %s\n\n" \
              % (model_df.columns, design_formula, e)
        raise Exception(err)

    print dmatrix.design_info.column_names
    print dmatrix

    # check the model for multicollinearity - Patsy takes care of this, but
    # just in case
    check_multicollinearity(np.array(dmatrix))

    # prepare for final stages
    column_names = dmatrix.design_info.column_names

    # what is this for?
    design_matrix = np.array(dmatrix, dtype=np.float16)
    
        
    # check to make sure there are more time points than EVs!
    if len(column_names) >= num_subjects:
        err = "\n\n[!] CPAC says: There are more EVs than there are " \
              "participants currently included in the model for %s. There " \
              "must be more participants than EVs in the design.\n\nNumber " \
              "of participants: %d\nNumber of EVs: %d\n\nEV/covariate list: "\
              "%s\n\nNote: If you specified to model group " \
              "variances separately, the amount of EVs can nearly double " \
              "once they are split along the grouping variable.\n\n" \
              "If the number of subjects is lower than the number of " \
              "subjects in your group analysis subject list, this may be " \
              "because not every subject in the subject list has an output " \
              "for %s in the individual-level analysis output directory.\n\n"\
              % (resource_id, num_subjects, len(column_names), column_names, \
                 resource_id)
        raise Exception(err)

    # time for contrasts
    contrasts_dict = None

    if ((custom_confile == None) or (custom_confile == '') or \
            ("None" in custom_confile) or ("none" in custom_confile)):

        # if no custom contrasts matrix CSV provided (i.e. the user
        # specified contrasts in the GUI)
        contrasts_list = group_config_obj.contrasts
        contrasts_dict = create_contrasts_dict(dmatrix, contrasts_list,
            resource_id)

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

    # send off the info so the FLAME input model files can be generated!
    mat_file, grp_file, con_file, fts_file = create_flame_model_files(dmatrix, \
        column_names, contrasts_dict, custom_confile, ftest_list, \
        group_config_obj.group_sep, grp_vector, group_config_obj.coding_scheme[0], \
        model_name, resource_id, model_path)

    dmat_csv_path = os.path.join(model_path, "design_matrix.csv")
    write_design_matrix_csv(dmatrix, model_df["Participant"], column_names, \
        dmat_csv_path)

    # workflow time
    wf_name = "%s_%s" % (resource_id, series_or_repeated_label)
    wf = pe.Workflow(name=wf_name)

    wf.base_dir = work_dir
    crash_dir = os.path.join(pipeline_config_obj.crashLogDirectory, \
                             "group_analysis", model_name)

    wf.config['execution'] = {'hash_method': 'timestamp', \
                              'crashdump_dir': crash_dir} 

    # gpa_wf
    # Creates the actual group analysis workflow
    gpa_wf = create_group_analysis(fTest, "gp_analysis_%s" % wf_name)

    gpa_wf.inputs.inputspec.merged_file = merge_file
    gpa_wf.inputs.inputspec.merge_mask = merge_mask

    gpa_wf.inputs.inputspec.z_threshold = z_threshold
    gpa_wf.inputs.inputspec.p_threshold = p_threshold
    gpa_wf.inputs.inputspec.parameters = (pipeline_config_obj.FSLDIR, \
                                          'MNI152')

    gpa_wf.inputs.inputspec.mat_file = mat_file
    gpa_wf.inputs.inputspec.con_file = con_file
    gpa_wf.inputs.inputspec.grp_file = grp_file

    if fTest:
        gpa_wf.inputs.inputspec.fts_file = fts_file      

    # ds
    # Creates the datasink node for group analysis
    ds = pe.Node(nio.DataSink(), name='gpa_sink')
     
    #     if c.mixedScanAnalysis == True:
    #         out_dir = re.sub(r'(\w)*scan_(\w)*(\d)*(\w)*[/]', '', out_dir)
              
    ds.inputs.base_directory = str(out_dir)
    ds.inputs.container = ''
        
    ds.inputs.regexp_substitutions = [(r'(?<=rendered)(.)*[/]','/'),
                                      (r'(?<=model_files)(.)*[/]','/'),
                                      (r'(?<=merged)(.)*[/]','/'),
                                      (r'(?<=stats/clusterMap)(.)*[/]','/'),
                                      (r'(?<=stats/unthreshold)(.)*[/]','/'),
                                      (r'(?<=stats/threshold)(.)*[/]','/'),
                                      (r'_cluster(.)*[/]',''),
                                      (r'_slicer(.)*[/]',''),
                                      (r'_overlay(.)*[/]','')]
   

    ########datasink connections#########
    #if fTest:
    #    wf.connect(gp_flow, 'outputspec.fts',
    #               ds, 'model_files.@0') 
        
    #wf.connect(gp_flow, 'outputspec.mat',
    #           ds, 'model_files.@1' )
    #wf.connect(gp_flow, 'outputspec.con',
    #           ds, 'model_files.@2')
    #wf.connect(gp_flow, 'outputspec.grp',
    #           ds, 'model_files.@3')
    wf.connect(gpa_wf, 'outputspec.merged',
               ds, 'merged')
    wf.connect(gpa_wf, 'outputspec.zstats',
               ds, 'stats.unthreshold')
    wf.connect(gpa_wf, 'outputspec.zfstats',
               ds,'stats.unthreshold.@01')
    wf.connect(gpa_wf, 'outputspec.fstats',
               ds,'stats.unthreshold.@02')
    wf.connect(gpa_wf, 'outputspec.cluster_threshold_zf',
               ds, 'stats.threshold')
    wf.connect(gpa_wf, 'outputspec.cluster_index_zf',
               ds,'stats.clusterMap')
    wf.connect(gpa_wf, 'outputspec.cluster_localmax_txt_zf',
               ds, 'stats.clusterMap.@01')
    wf.connect(gpa_wf, 'outputspec.overlay_threshold_zf',
               ds, 'rendered')
    wf.connect(gpa_wf, 'outputspec.rendered_image_zf',
               ds, 'rendered.@01')
    wf.connect(gpa_wf, 'outputspec.cluster_threshold',
               ds,  'stats.threshold.@01')
    wf.connect(gpa_wf, 'outputspec.cluster_index',
               ds, 'stats.clusterMap.@02')
    wf.connect(gpa_wf, 'outputspec.cluster_localmax_txt',
               ds, 'stats.clusterMap.@03')
    wf.connect(gpa_wf, 'outputspec.overlay_threshold',
               ds, 'rendered.@02')
    wf.connect(gpa_wf, 'outputspec.rendered_image',
               ds, 'rendered.@03')
       
    ######################################

    # Run the actual group analysis workflow
    wf.run()

    print "\n\nWorkflow finished for model %s\n\n" % wf_name



def run(config, subject_infos, resource):
    import re
    import commands
    commands.getoutput('source ~/.bashrc')
    import os
    import sys
    import pickle
    import yaml
    
    c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))
    
    prep_group_analysis_workflow(c, pickle.load(open(resource, 'r') ), \
        pickle.load(open(subject_infos, 'r')))



