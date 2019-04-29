

from CPAC.QPP.prep_QPP import split_subdfs

#from CPAC.pipeline.cpac_group_runner import gather_outputs
#df_dct = gather_outputs(pipeline_folder, ['functional_nuisance_residuals'], None, False, False, get_func=True)

import pickle
with open('/cpac_team/cpac/resources/pipeline_folder_sess-scan.pkl', 'r') as f:
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
    qpp_dict = split_subdfs(df_dct, scan_inclusion=['rest_run-1'])
    print('with only scan rest_run-1:')
    print(qpp_dict)
    print('\n\n')

def test_sess_split():
    qpp_dict = split_subdfs(df_dct, sess_inclusion=['2'])
    print('with session 2 only:')
    print(qpp_dict)
    print('\n\n')

def test_sess_scan_split():
    qpp_dict = split_subdfs(df_dct, sess_inclusion=['2'],
                            scan_inclusion=['rest_run-2'])
    print('with session 2 only; and with scan rest_run-2 only:')
    print(qpp_dict)
    print('\n\n')


test_basic_df_split()
test_scan_split()
test_sess_split()
test_sess_scan_split()

