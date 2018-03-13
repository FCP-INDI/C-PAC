

def write_design_matrix_csv(patsy_dmatrix, participant_column, column_names,
                            outfile_path):

    import pandas as pd

    try:
        # TODO
        # patsy dmatrix might actually already be a dataframe- rename this
        # variable please
        patsy_dmatrix.to_csv(outfile_path)
    except:
        group_model_dataframe = pd.DataFrame(data=patsy_dmatrix,
                                             index=participant_column,
                                             columns=column_names)
        group_model_dataframe.to_csv(outfile_path)

    #except Exception as e:
    #   err = "\n\n[!] Could not write the design matrix dataframe to the " \
    #          "CSV file: %s\n\nError details: %s\n\n" % (outfile_path, e)
    #    raise Exception(err)


def write_custom_readme_file():
    pass

