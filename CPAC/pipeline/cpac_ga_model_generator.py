import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys
import glob

from CPAC.utils.datasource import create_grp_analysis_dataflow
from CPAC.utils import Configuration
from CPAC.utils.utils import prepare_gp_links
from CPAC.group_analysis import create_group_analysis



def prep_group_analysis_workflow(c, resource, subject_infos):
    
    #
    # this function runs once per output file during group analysis
    #

    # p_id = a list of pipeline IDs, i.e. the name of the output folder for
    #        the strat
    
    # s_ids = a list of all the subject IDs

    # scan_ids = a list of scan IDs

    # s_paths = a list of all of the filepaths of this particular output
    #           file that prep_group_analysis_workflow is being called for

    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))

    # set this to False for now
    fTest = False


    for config in c.modelConfigs:

       
        try:
            group_conf = Configuration(yaml.load(open(os.path.realpath(config), 'r')))
        except:
            raise Exception("Error in reading %s configuration file" % config)


        # only run model if the derivative was selected for it - this change
        # will allow different models to have different derivative lists
        if resource in group_conf.derivative_list:
        
            group_sublist_file = open(group_conf.subject_list, 'r')

            group_sublist_items = group_sublist_file.readlines()

            group_sublist = [line.rstrip('\n') for line in group_sublist_items \
                                  if not (line == '\n') and not line.startswith('#')]

            # list of subjects for which paths which DO exist
            exist_paths = []

            # paths to the actual derivatives for those subjects
            derivative_paths = []



            ''' begin iteration through group subject list for processing '''

            for ga_sub in group_sublist:

                # ga_sub = subject ID taken off the group analysis subject list

                # let's check to make sure the subject list is formatted for
                # repeated measures properly if repeated measures is enabled
                # and vice versa
                if (group_conf.repeated_measures == True) and (',' not in ga_sub):
                    print '\n\n'
                    print '[!] CPAC says: The group analysis subject list ' \
                            'is not in the appropriate format for repeated ' \
                            'measures.\n'
                    print 'Please use the appropriate format as described in ' \
                            'the CPAC User Guide or turn off Repeated Measures ' \
                            'in the CPAC pipeline configuration editor, found ' \
                            'in the \'Group Analysis Settings\' tab of the ' \
                            'pipeline configuration editor.\n'
                    print 'NOTE: CPAC generates a properly-formatted group ' \
                            'analysis subject list meant for running repeated ' \
                            'measures when you create your original subject ' \
                            'list. Look for \'subject_list_group_analysis_' \
                            'repeated_measures.txt\' in the directory where ' \
                            'you created your subject list.\n\n'
                    raise Exception

                elif (group_conf.repeated_measures == False) and (',' in ga_sub):
                    print '\n\n'
                    print '[!] CPAC says: It looks like your group analysis ' \
                            'subject list is formatted for running repeated ' \
                            'measures, but \'Run Repeated Measures\' is not ' \
                            'enabled in the pipeline configuration, found in ' \
                            'the \'Group Analysis Settings\' tab of the ' \
                            'pipeline configuration editor.\n'
                    print 'Double-check your pipeline configuration?\n\n'
                    raise Exception



                ''' process subject ids for repeated measures, if it is on '''
                # if repeated measures is being run and the subject list
                # is a list of subject IDs and scan IDs concatenated
                if (group_conf.repeated_measures == True):

                    # sub.count(',') equals 1 when there is either multiple scans
                    # or multiple sessions but not both, for repeated measures

                    # sub.count(',') equals 2 when there are multiple sessions
                    # AND scans, for repeated measures

                    if ga_sub.count(',') == 1:
                        sub_id = ga_sub.split(',',1)[0]
                        other_id = ga_sub.split(',',1)[1]

                    elif ga_sub.count(',') == 2:
                        sub_id = ga_sub.split(',',2)[0]
                        scan_id = ga_sub.split(',',2)[1]
                        session_id = ga_sub.split(',',2)[2]



                ''' drop subjects from the group subject list '''
                # check the path files in path_files_here folder in the
                # subject's output folder - and drop any subjects from the
                # group analysis subject list which do not exist in the paths
                # to the output files

                '''
                REVISIT THIS LATER to establish a potentially better way to
                pull output paths (instead of path_files_here)
                '''

                for path in s_paths:

                    if (group_conf.repeated_measures == True):

                        if ga_sub.count(',') == 1:
                            if (sub_id in path) and (other_id in path):
                                exist_paths.append(ga_sub)
                                derivative_paths.append(path)

                        elif ga_sub.count(',') == 2:
                            if (sub_id in path) and (scan_id in path) and \
                                    (session_id in path):
                                exist_paths.append(ga_sub)
                                derivative_paths.append(path)

                    else:

                        if ga_sub in path:
                            exist_paths.append(ga_sub)


                # END subject-dropping!

            if len(derivative_paths) == 0:
                print '\n\n\n[!] CPAC says: None of the subjects listed in the ' \
                      'group analysis subject list were found to have outputs ' \
                      'produced by individual-level analysis.\n\nEnsure that ' \
                      'the subjects listed in your group analysis subject list ' \
                      'are the same as the ones included in the individual-' \
                      'level analysis you are running group-level analysis for.' \
                      '\n\n\n'
                raise Exception

            ''' END subject list iteration '''
 


            # check to see if any derivatives of subjects are missing
            if len(list(set(subject_list) - set(exist_paths))) >0:
                print "List of outputs missing for subjects:"
                print list(set(subject_list) - set(exist_paths))
                print "..for derivatives:"
                print resource
                print "..at paths:"
                print os.path.dirname(s_paths[0]).replace(s_ids[0], '*')

        

            # create the path string for the group analysis output
            out_dir = os.path.dirname(s_paths[0]).split(p_id[0] + '/')
            out_dir = os.path.join(conf.output_dir, out_dir[1])
            out_dir = out_dir.replace(s_ids[0], 'group_analysis_results_%s/_grp_model_%s'%(p_id[0],group_conf.model_name))

            mod_path = os.path.join(out_dir, 'model_files')


            if not os.path.isdir(mod_path):
                os.makedirs(mod_path)

        
            ''' write the new subject list '''
            new_sub_file = os.path.join(mod_path, os.path.basename(group_conf.subject_list))

            try:

                f = open(new_sub_file, 'w')
         
                for sub in exist_paths:
                    print >>f, sub
        
                f.close()

            except:

                print "Error: Could not open subject list file: ", new_sub_file
                raise Exception


            group_conf.update('subject_list',new_sub_file)

            sub_id_label = group_conf.subject_id_label


            # Run 'create_fsl_model' script to extract phenotypic data from
            # the phenotypic file for each of the subjects in the subject list



            ''' get the motion statistics parameter file, if present '''
            # get the parameter file so it can be passed to create_fsl_model.py
            # so MeanFD or other measures can be included in the design matrix
            parameter_file = os.path.join(c.outputDirectory, p_id[0], '%s_threshold_%s_all_params.csv'%(scan_ids[0].strip('_'),threshold_val))

            if 1 in c.runGenerateMotionStatistics:

                if not os.path.exists(parameter_file):
                    print '\n\n[!] CPAC says: Could not open the parameter file. ' \
                          'If Generate Motion Statistics is enabled, this can ' \
                          'usually be found in the output directory of your ' \
                          'individual-level analysis runs.\n'
                    print 'Path not found: ', parameter_file, '\n\n'
                    raise Exception

            elif (1 not in c.runGenerateMotionStatistics) and (os.path.exists(parameter_file)):

                if not os.path.exists(parameter_file):
                    print '\n\n[!] CPAC says: Could not open the parameter file. ' \
                          'If Generate Motion Statistics is enabled, this can ' \
                          'usually be found in the output directory of your ' \
                          'individual-level analysis runs.\n'
                    print 'Path not found: ', parameter_file, '\n\n'
                    raise Exception

            else:

                def no_measures_error(measure):
                    print '\n\n[!] CPAC says: The measure %s was included in ' \
                          'your group analysis design matrix formula, but ' \
                          'Generate Motion Statistics was not run during ' \
                          'individual-level analysis.\n' % measure
                    print 'Please run Generate Motion Statistics if you wish ' \
                          'to include this measure in your model.\n'
                    print 'If you HAVE completed a run with this option ' \
                          'enabled, then you are seeing this error because ' \
                          'the motion parameter file normally created by this ' \
                          'option is missing.\n\n'
                    raise Exception

                for measure in measure_list:
                    if (measure in group_conf.design_formula):
                        no_measures_error(measure)

                parameter_file = None




            # path to the pipeline folder to be passed to create_fsl_model.py
            # so that certain files like output_means.csv can be accessed
            pipeline_path = os.path.join(c.outputDirectory, p_id[0])

            # the current output that cpac_group_analysis_pipeline.py and
            # create_fsl_model.py is currently being run for
            current_output = s_paths[0].replace(pipeline_path, '').split('/')[2]


            ''' merge the remaining subjects for this current output '''
            # then, take the group mask, and iterate over the list of subjects
            # remaining to extract the mean of each subject using the group
            # mask

            merge_input = " "
            merge_output = modpath + "/" + current_output + "_merged.nii.gz"
            merge_mask_output = modpath + "/" + current_output + "_merged_mask.nii.gz"

            for derivative_path in derivative_paths:
                merge_input = merge_input + derivative_path

            merge_string = "fslmerge -t %s %s" % (merge_output, merge_input)

            # MERGE the remaining outputs
            try:
                commands.getoutput(merge_string)
            except Exception as e:
                print "FSL Merge failed for output: %s\n" % current_output
                print "error: %s\n\n" % e
                raise

            merge_mask_string = "fslmaths %s -abs -Tmin -bin %s" % (merge_output, merge_mask_output)

            # CREATE A MASK of the merged file
            try:
                commands.getoutput(merge_mask_string)
            except Exception as e:
                print "FSL Mask failed for output: %s\n" % current_output
                print "error: %s\n\n" % e
                raise


            derivative_mean = {}

            # CALCULATE THE MEANS of each remaining output using the group mask
            for derivative_path in derivative_paths:

                try:
                    3dmaskave_output = commands.getoutput("3dmaskave -mask %s %s" % (merge_mask_output, derivative_path))
                except Exception as e:
                    print "AFNI 3dmaskave failed for output: %s\n" % current_output
                    print "error: %s\n\n" % e
                    raise

                # get the subject ID of the current derivative path reliably
                derivative_path_subID = derivative_path.replace(pipeline_path,"").strip("/").split("/")[0]

                # this crazy-looking command simply extracts the mean from the
                # verbose AFNI 3dmaskave output string
                derivative_mean[derivative_path_subID] = 3dmaskave_output.split("\n")[-1].split(" ")[0]

                # derivative_mean is now something like this:
                # { 'sub001': 0.3124, 'sub002': 0.2981, .. }



            ''' run create_fsl_model.py to generate the group analysis models '''
            try:

                from CPAC.utils import create_fsl_model
                create_fsl_model.run(conf, fTest, parameter_file, derivative_mean, pipeline_path, current_output, True)

            except Exception, e:

                print "FSL Group Analysis model not successfully created - error in create_fsl_model script"
                #print "Error ->", e
                raise

