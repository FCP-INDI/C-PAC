import re
import os

import glob
#from CPAC.utils import Configuration as c 
from CPAC.utils.datasource import create_grp_analysis_dataflow
from CPAC.utils.utils import prepare_gp_links
from CPAC.pipeline.cpac_group_runner import load_config_yml
from CPAC.group_analysis import create_group_analysis
from CPAC.pipeline.cpac_group_runner import prep_feat_inputs

def main():

    from multiprocessing import Process
    import argparse
    import sys 
    parser = argparse.ArgumentParser()
    parser.add_argument("group_config_file", type=str, help='provide the path top the group config file')
    parser.add_argument("pipeline_output_folder", type=str,help='provide the path to the output folder of your group config files')

    args = parser.parse_args()

    # let's get the show on the road
    procss = []
    #group_config_file = args.group_config_file_path
    #pipeline_output_folder = args.pipeline_output_folder_path


        # create the analysis DF dictionary
    analysis_dict = prep_feat_inputs(args.group_config_file,args.pipeline_output_folder)
    print analysis_dict
    for unique_resource_id in analysis_dict.keys():
        # unique_resource_id is a 5-long tuple:
        #    ( model name, group model config file, output measure name,
        #          preprocessing strategy string,
        #          series_id or "repeated_measures" )

        model_name = unique_resource_id[0]
        group_config_file = unique_resource_id[1]
        resource_id = unique_resource_id[2]
        preproc_strat = unique_resource_id[3]
        series_or_repeated = unique_resource_id[4]

        model_df = analysis_dict[unique_resource_id]

        if not c.runOnGrid:
            from CPAC.pipeline.cpac_ga_model_generator import \
                prep_group_analysis_workflow

            procss.append(Process(target=prep_group_analysis_workflow,
                                  args=(model_df, model_name,
                                        group_config_file, resource_id,
                                        preproc_strat,
                                        series_or_repeated)))
        else:
            print "\n\n[!] CPAC says: Group-level analysis has not yet " \
                  "been implemented to handle runs on a cluster or " \
                  "grid.\n\nPlease turn off 'Run CPAC On A Cluster/" \
                  "Grid' in order to continue with group-level " \
                  "analysis. This will submit the job to only one " \
                  "node, however.\n\nWe will update users on when this " \
                  "feature will be available through release note " \
                  "announcements.\n\n"

    manage_processes(procss, c.outputDirectory, c.numGPAModelsAtOnce)
    

if __name__ == "__main__":
  main()