import re
import os
import glob
#from CPAC.utils import Configuration as c 
from CPAC.utils.datasource import create_grp_analysis_dataflow
from CPAC.utils.utils import prepare_gp_links
from CPAC.pipeline.cpac_group_runner import load_config_yml
from CPAC.group_analysis import create_group_analysis
from CPAC.pipeline.cpac_group_runner import prep_feat_inputs
from CPAC.pipeline.cpac_ga_model_generator import prep_group_analysis_workflow


def prep():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("group_config_file", type=str,
                        help='provide the path to the group config file')
    args = parser.parse_args()

    group_config_obj = load_config_yml(args.group_config_file)
    pipeline_output_folder = group_config_obj.pipeline_dir

    analysis_dict = prep_feat_inputs(args.group_config_file,
                                     pipeline_output_folder)
    dmat_csv_path = ""
    new_sub_file = ""

    for unique_resource_id in analysis_dict.keys():

        # unique_resource_id is a 6-long tuple:
        #    ( model name, group model config file, output measure name,
        #          preprocessing strategy string, session_id,
        #          series_id or "repeated_measures" )

        model_name = unique_resource_id[0]

        group_config_file = unique_resource_id[1]
        resource_id = unique_resource_id[2]
        preproc_strat = unique_resource_id[3]
        session = unique_resource_id[4]
        series_or_repeated = unique_resource_id[5]
        model_df = analysis_dict[unique_resource_id]

        dmat_csv_path, new_sub_file, empty_csv = prep_group_analysis_workflow(model_df,
                                                                              model_name,
                                                                              args.group_config_file,
                                                                              resource_id,
                                                                              preproc_strat,
                                                                              session,
                                                                              series_or_repeated)

    return dmat_csv_path, new_sub_file, empty_csv

  
if __name__ == "__main__":
    prep()
