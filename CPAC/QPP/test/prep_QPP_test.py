from CPAC.QPP.prep_QPP import split_subdfs

#from CPAC.pipeline.cpac_group_runner import gather_outputs
#df_dct = gather_outputs(pipeline_folder, ['functional_nuisance_residuals'], None, False, False, get_func=True)

import pickle
import pandas as pd
with open('/home/nrajamani/C-PAC/CPAC/QPP/op_dict.pkl', 'r') as f:
    df_dct = pickle.load(f)

print('\n\nInput dataframe:')
print(df_dct)
print('\n\n')

def test_basic_df_split():
    qpp_dict = split_subdfs(df_dct)
    print('with all sessions and all scans:')
    print(qpp_dict)
    print('\n\n')

def test_scan_split():
    qpp_dict = split_subdfs(df_dct,None,['rest_run-1'],['None'])
    print('with only scan rest_run-1:')
    print(qpp_dict)
    print('\n\n')

def test_sess_split():
    qpp_dict = split_subdfs(df_dct, sess_inclusion=None)
    print('with session 2 only:')
    print(qpp_dict)
    print('\n\n')

def test_sess_scan_split():
    qpp_dict = split_subdfs(df_dct,['2'],['rest_run-2'],['None'])
    print('with session 2 only; and with scan rest_run-2 only:')
    print(qpp_dict)
    print('\n\n')


#test_basic_df_split()
#test_scan_split()
#test_sess_split()
test_sess_scan_split()
