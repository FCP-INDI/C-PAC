

def write_design_matrix_csv(group_model_dataframe, outfile_path):

    import os
    import pandas as pd

    try:
    	group_model_dataframe.to_csv(outfile_path)
    except Exception as e:
    	err = "\n\n[!] Could not write the design matrix dataframe to the " \
    	      "CSV file: %s\n\nError details: %s\n\n" % (outfile_path, e)



def write_custom_readme_file():
    pass
	