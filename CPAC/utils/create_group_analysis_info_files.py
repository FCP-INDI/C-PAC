

def write_design_matrix_csv(patsy_dmatrix, participant_column, column_names,
	outfile_path):

    import os
    import patsy
    import pandas as pd

    try:
    	group_model_dataframe = pd.DataFrame(data=patsy_dmatrix, \
    		                                 index=participant_column, \
    		                                 columns=column_names)
    	group_model_dataframe.to_csv(outfile_path)
    except Exception as e:
    	err = "\n\n[!] Could not write the design matrix dataframe to the " \
    	      "CSV file: %s\n\nError details: %s\n\n" % (outfile_path, e)



def write_custom_readme_file():
    pass
	