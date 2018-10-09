

def run_gather_outputs_func(pipeline_out_dir):
    from CPAC.pipeline import cpac_group_runner as cgr
    df_dct = cgr.gather_outputs(pipeline_out_dir, ["functional_to_standard"],
                                None, False, False, get_func=True)
    print(df_dct)

