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
        err = "\n\n[!] CPAC says: Could not write new participant list for "\
              "current model and derivative in group-level analysis. Ensure "\
              "you have write access to the directory:\n%s\n\nError " \
              "details: %s\n\n" \
              % (current_mod_path, e)
        raise Exception(err)

    return new_sub_file



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
    mm_df["Measure_Mean"] = mm_df["Measure_Mean"].sub(mm_df["Measure_Mean"].mean())
    
    model_df = pd.merge(model_df, mm_df, how="inner", on=["Raw_Filepath"])
    
    return model_df



def check_mask_file_resolution(data_file, roi_mask, out_dir, output_id=None):

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

        resample_str = ["flirt", "-in", roi_mask, "-ref", roi_mask, \
                        "-applyisoxfm", raw_file_dims[0], "-out", \
                        resampled_outfile]

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
              % (roi_mask, merge_mask, e)
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
        roi_means_string = retcode.split(raw_file)[1].rstrip("\n")

        if "\t[0]?\t" in roi_means_string:
            roi_means_string = roi_means_string.replace("\t[0]?\t","")

        roi_means_string_list = roi_means_string.split("\t")

        # check
        roi_means_list = []
        for roi_mean in roi_means_string_list:
            try:
                roi_means_list.append(float(roi_mean))
            except:
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



def prep_group_analysis_workflow(model_df, pipeline_config_obj, \
    model_name, group_config_obj, resource_id, preproc_strat, \
    series_or_repeated_label):
    
    #
    # this function runs once per derivative type and preproc strat combo
    # during group analysis
    #

    import os

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.io as nio

    pipeline_ID = pipeline_config_obj.pipeline_name

    # get thresholds
    z_threshold = float(group_config_obj.z_threshold[0])

    p_threshold = float(group_config_obj.p_threshold[0])

    sub_id_label = group_config_obj.subject_id_label

    # determine if f-tests are included or not
    custom_confile = group_config_obj.custom_contrasts

    if ((custom_confile == None) or (custom_confile == '') or \
            ("None" in custom_confile) or ("none" in custom_confile)):

        if (len(group_config_obj.f_tests) == 0) or \
            (group_config_obj.f_tests == None):
            fTest = False
        else:
            fTest = True

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

    model_path = os.path.join(out_dir, 'model_files')

    # generate working directory for this output's group analysis run
    work_dir = os.path.join(c.workingDirectory, "group_analysis", model_name,\
        resource_id, series_or_repeated_label, preproc_strat)

    log_dir = os.path.join(out_dir, 'logs', resource_id, \
        'model_%s' % model_name)

    # create the actual directories
    if not os.path.isdir(model_path):
        try:
            os.makedirs(model_path)
        except Exception as e:
            err = "\n\n[!] Could not create the group analysis output " \
                  "directories.\n\nAttempted directory creation: %s\n\n" \
                  "Error details: %s\n\n" % (model_path, e)
            raise Exception(err)

    if not os.path.isdir(work_dir):
        try:
            os.makedirs(work_dir)
        except Exception as e:
            err = "\n\n[!] Could not create the group analysis working " \
                  "directories.\n\nAttempted directory creation: %s\n\n" \
                  "Error details: %s\n\n" % (model_path, e)
            raise Exception(err)

    if not os.path.isdir(log_dir):
        try:
            os.makedirs(log_dir)
        except Exception as e:
            err = "\n\n[!] Could not create the group analysis logfile " \
                  "directories.\n\nAttempted directory creation: %s\n\n" \
                  "Error details: %s\n\n" % (model_path, e)
            raise Exception(err)


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

    group_conf.update('participant_list',new_sub_file)


    # start processing the dataframe further
    design_formula = group_config_obj.design_formula

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
    if group_config_obj.mean_mask[0] == "Group Mask":
        merge_mask_outfile = os.path.basename(merge_file) + "_mask.nii.gz"
        merge_mask = create_merged_mask(merge_file, merge_mask_outfile)

    # calculate measure means, and demean
    if "Measure_Mean" in design_formula:
        model_df = calculate_measure_mean_in_df(model_df, merge_mask)

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
                                              custom_roi_mask, model_path, \
                                              resource_id)

        # if using group merged mask, trim the custom ROI mask to be within
        # its constraints
        if merge_mask:
            output_mask = os.path.join(model_path, "group_masked_%s" \
                                       % os.path.basename(input_mask))
            roi_mask = trim_mask(roi_mask, merge_mask, output_mask)

        # calculate
        model_df = calculate_custom_roi_mean_in_df(model_df, roi_mask)   

    


    # modeling group variances separately

    # add repeated measures 1's matrices

    # patsify model DF, drop columns not in design formula

    # process contrasts


        
    wf = pe.Workflow(name=resource_id)

    wf.base_dir = work_dir
    crash_dir = os.path.join(pipeline_config_obj.crashLogDirectory, \
                             "group_analysis", model_name)

    wf.config['execution'] = {'hash_method': 'timestamp', \
                              'crashdump_dir': crash_dir}








    if "Measure_Mean" in design_formula:
        measure_mean = pe.Node(util.Function(input_names=['model_df',
                                                          'merge_mask'],
                                       output_names=['model_df'],
                                       function=calculate_measure_mean_in_df),
                                       name='measure_mean')
        measure_mean.inputs.model_df = model_df

        wf.connect(merge_mask, "out_file", measure_mean, "merge_mask")


    if "Custom_ROI_Mean" in design_formula:
        roi_mean = pe.Node(util.Function())


    group_config_obj.custom_roi_mask
    






    #----------------

    import yaml
    import pandas as pd


    # load group analysis model configuration file
    try:
        with open(os.path.realpath(group_config_file),"r") as f:
            group_conf = Configuration(yaml.load(f))
    except Exception as e:
        err_string = "\n\n[!] CPAC says: Could not read group model " \
                     "configuration YML file. Ensure you have read access " \
                     "for the file and that it is formatted properly.\n\n" \
                     "Configuration file: %s\n\nError details: %s" \
                     % (group_config_file, e)
        raise Exception(err_string)


    # gather all of the information
    # - lists of all the participant unique IDs (participant_site_session) and
    # of all of the series IDs present in output_file_list
    # - also returns the pipeline ID
    new_participant_list, all_series_names, pipeline_ID = \
        gather_new_participant_list(output_path_file, output_file_list)

     

      

    # create the path string for the group analysis output
    #    replicate the directory path of one of the participant's output
    #    folder path to the derivative's file, but replace the participant ID
    #    with the group model name
    #        this is to ensure nothing gets overwritten between strategies
    #        or thresholds, etc.
    out_dir = os.path.dirname(output_file_list[0]).split(pipeline_ID + '/')
    out_dir = out_dir[1].split(out_dir[1].split("/")[-1])[0]
    out_dir = os.path.join(group_conf.output_dir, out_dir)
    out_dir = out_dir.replace(new_participant_list[0], \
                  'group_analysis_results_%s/_grp_model_%s' \
                  % (pipeline_ID, group_conf.model_name))

    # !!!!!!!!!!
    if (group_conf.repeated_measures == True) and (series_ids[0] != None):
        out_dir = out_dir.replace(series_ids[0] + "/", "multiple_series")

    # create model file output directories
    model_out_dir = os.path.join(group_conf.output_dir, \
        'group_analysis_results_%s/_grp_model_%s' \
        %(pipeline_ID, group_conf.model_name))

    mod_path = os.path.join(model_out_dir, 'model_files')

    if not os.path.isdir(mod_path):
        os.makedirs(mod_path)

    # current_mod_path = folder under
    #   "/gpa_output/_grp_model_{model name}/model_files/{current derivative}"
    current_mod_path = os.path.join(mod_path, resource)

    if not os.path.isdir(current_mod_path):
        os.makedirs(current_mod_path)

        
    # create new subject list based on which subjects are left after checking
    # for missing outputs
    new_sub_file = write_new_sub_file(current_mod_path, \
                       group_conf.subject_list, new_participant_list)

    group_conf.update('subject_list',new_sub_file)


    # create new design matrix with only the subjects that are left






    # Run 'create_fsl_model' script to extract phenotypic data from
    # the phenotypic file for each of the subjects in the subject list

    # get the motion statistics parameter file, if present
    # get the parameter file so it can be passed to create_fsl_model.py
    # so MeanFD or other measures can be included in the design matrix


    ''' okay, here we go... how are we handling series? because here it needs to take in '''
    ''' the appropriate series to get the appropriate parameter file ! ! ! '''

    ''' MAY HAVE TO GO BACK ON THIS, and just have one series sent in per this function...'''

    power_params_files = {}

    measure_list = ['MeanFD_Power', 'MeanFD_Jenkinson', 'MeanDVARS']

    for measure in measure_list:
    
        if measure in group_conf.design_formula:

            for series_id in all_series_names:

                parameter_file = os.path.join(c.outputDirectory, \
                                              pipeline_ID, \
                                              '%s%s_all_params.csv' % \
                                              (series_id.strip('_'), \
                                              threshold_val))

                if not os.path.exists(parameter_file):
                    err = "\n\n[!] CPAC says: Could not find or open the motion "\
                          "parameter file. This is necessary if you have " \
                          "included any of the MeanFD measures in your group " \
                          "model.\n\nThis file can usually be found in the " \
                          "output directory of your individual-level analysis " \
                          "runs. If it is not there, double-check to see if " \
                          "individual-level analysis had completed successfully."\
                          "\n\nPath not found: %s\n\n" % parameter_file
                    raise Exception(err)


                power_params_files[series_id] = parameter_file
                

            break
            
    else:
    
        power_params_files = None



    # path to the pipeline folder to be passed to create_fsl_model.py
    # so that certain files like output_means.csv can be accessed
    pipeline_path = os.path.join(c.outputDirectory, pipeline_ID)

    # generate working directory for this output's group analysis run
    workDir = '%s/group_analysis/%s/%s' % (c.workingDirectory, \
                                               group_conf.model_name, \
                                               resource)
            
    # this makes strgy_path basically the directory path of the folders after
    # the resource/derivative folder level         
    strgy_path = os.path.dirname(output_file_list[0]).split(resource)[1]

    # get rid of periods in the path
    for ch in ['.']:
        if ch in strgy_path:
            strgy_path = strgy_path.replace(ch, "")
                
    # create nipype-workflow-name-friendly strgy_path
    # (remove special characters)
    strgy_path_name = strgy_path.replace('/', "_")

    workDir = workDir + '/' + strgy_path_name



    # merge the subjects for this current output
    # then, take the group mask, and iterate over the list of subjects
    # to extract the mean of each subject using the group mask
    merge_output, merge_mask_output, merge_output_dir = \
        create_merged_files(workDir, resource, output_file_list)

    
    # CALCULATE THE MEANS of each output using the group mask
    derivative_means_dict, roi_means_dict = \
        calculate_output_means(resource, output_file_list, \
                               group_conf.mean_mask, \
                               group_conf.design_formula, \
                               group_conf.custom_roi_mask, pipeline_path, \
                               merge_output_dir, c.identityMatrix)


    measure_dict = {}

    # extract motion measures from CPAC-generated power params file
    if power_params_files != None:
        for param_file in power_params_files.values():
            new_measure_dict = get_measure_dict(param_file)
            measure_dict.update(new_measure_dict)


    # combine the motion measures dictionary with the measure_mean
    # dictionary (if it exists)
    if derivative_means_dict:
        measure_dict["Measure_Mean"] = derivative_means_dict

    # run create_fsl_model.py to generate the group analysis models
    
    from CPAC.utils import create_fsl_model, kill_me
    create_fsl_model.run(group_conf, resource, parameter_file, \
                             derivative_means_dict, roi_means_dict, \
                                 current_mod_path, True)


    # begin GA workflow setup

    if not os.path.exists(new_sub_file):
        raise Exception("path to input subject list %s is invalid" % new_sub_file)
        
    #if c.mixedScanAnalysis == True:
    #    wf = pe.Workflow(name = 'group_analysis/%s/grp_model_%s'%(resource, os.path.basename(model)))
    #else:

    wf = pe.Workflow(name = resource)

    wf.base_dir = workDir
    wf.config['execution'] = {'hash_method': 'timestamp', 'crashdump_dir': os.path.abspath(c.crashLogDirectory)}
    log_dir = os.path.join(group_conf.output_dir, 'logs', 'group_analysis', resource, 'model_%s' % (group_conf.model_name))
        

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    else:
        pass


    # gp_flow
    # Extracts the model files (.con, .grp, .mat, .fts) from the model
    # directory and sends them to the create_group_analysis workflow gpa_wf

    gp_flow = create_grp_analysis_dataflow("gp_dataflow_%s" % resource)
    gp_flow.inputs.inputspec.grp_model = os.path.join(mod_path, resource)
    gp_flow.inputs.inputspec.model_name = group_conf.model_name
    gp_flow.inputs.inputspec.ftest = fTest
  

    # gpa_wf
    # Creates the actual group analysis workflow

    gpa_wf = create_group_analysis(fTest, "gp_analysis_%s" % resource)

    gpa_wf.inputs.inputspec.merged_file = merge_output
    gpa_wf.inputs.inputspec.merge_mask = merge_mask_output

    gpa_wf.inputs.inputspec.z_threshold = z_threshold
    gpa_wf.inputs.inputspec.p_threshold = p_threshold
    gpa_wf.inputs.inputspec.parameters = (c.FSLDIR, 'MNI152')
    
   
    wf.connect(gp_flow, 'outputspec.mat',
               gpa_wf, 'inputspec.mat_file')
    wf.connect(gp_flow, 'outputspec.con',
               gpa_wf, 'inputspec.con_file')
    wf.connect(gp_flow, 'outputspec.grp',
                gpa_wf, 'inputspec.grp_file')
           
    if fTest:
        wf.connect(gp_flow, 'outputspec.fts',
                   gpa_wf, 'inputspec.fts_file')
        

    # ds
    # Creates the datasink node for group analysis
       
    ds = pe.Node(nio.DataSink(), name='gpa_sink')
     
    if 'sca_roi' in resource:
        out_dir = os.path.join(out_dir, \
            re.search('sca_roi_(\d)+',os.path.splitext(os.path.splitext(os.path.basename(output_file_list[0]))[0])[0]).group(0))
            
            
    if 'dr_tempreg_maps_zstat_files_to_standard_smooth' in resource:
        out_dir = os.path.join(out_dir, \
            re.search('temp_reg_map_z_(\d)+',os.path.splitext(os.path.splitext(os.path.basename(output_file_list[0]))[0])[0]).group(0))
            
            
    if 'centrality' in resource:
        names = ['degree_centrality_binarize', 'degree_centrality_weighted', \
                 'eigenvector_centrality_binarize', 'eigenvector_centrality_weighted', \
                 'lfcd_binarize', 'lfcd_weighted']

        for name in names:
            if name in os.path.basename(output_file_list[0]):
                out_dir = os.path.join(out_dir, name)
                break

    if 'tempreg_maps' in resource:
        out_dir = os.path.join(out_dir, \
            re.search('\w*[#]*\d+', os.path.splitext(os.path.splitext(os.path.basename(output_file_list[0]))[0])[0]).group(0))
        
#     if c.mixedScanAnalysis == True:
#         out_dir = re.sub(r'(\w)*scan_(\w)*(\d)*(\w)*[/]', '', out_dir)
              
    ds.inputs.base_directory = out_dir
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
    if fTest:
        wf.connect(gp_flow, 'outputspec.fts',
                   ds, 'model_files.@0') 
        
    wf.connect(gp_flow, 'outputspec.mat',
               ds, 'model_files.@1' )
    wf.connect(gp_flow, 'outputspec.con',
               ds, 'model_files.@2')
    wf.connect(gp_flow, 'outputspec.grp',
               ds, 'model_files.@3')
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

    
    print "\n\nWorkflow finished for model %s and resource %s\n\n" \
          % (os.path.basename(group_conf.output_dir), resource)



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



