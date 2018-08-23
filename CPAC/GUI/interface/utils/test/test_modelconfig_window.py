
import wx
from CPAC.GUI.interface.utils import modelconfig_window


def run_get_pheno_header(new_modelconfig, ga_settings):
    with open(ga_settings['pheno_file'], 'r') as f:
        new_modelconfig.get_pheno_header(f)
    return new_modelconfig.phenoHeaderItems

def run_read_phenotypic(new_modelconfig, ga_settings):
    pheno_dct = new_modelconfig.read_phenotypic(ga_settings['pheno_file'],
                                                ga_settings['ev_selections'])
    return pheno_dct


ga_settings = {'pheno_file': '/Users/steven.giavasis/run/cpac/v120/configs/participants.tsv',
               'participant_id_label': 'participant_id'}

app = wx.PySimpleApp()
frame = wx.Frame(None)
new_modelconfig = modelconfig_window.ModelConfig(frame,
                                                 gpa_settings=ga_settings)

pheno_header_items = run_get_pheno_header(new_modelconfig, ga_settings)
print('Output:\n{0}'.format(pheno_header_items))

pheno_dct = run_read_phenotypic(new_modelconfig, ga_settings)
print('Output:\n{0}'.format(pheno_dct))

ga_settings = {'pheno_file': '/Users/steven.giavasis/run/cpac/pheno_test_benchmark_12.csv',
               'participant_id_label': 'Participant',
               'ev_selections': {'categorical': ['Diagnosis', 'Sex']}}

new_modelconfig = modelconfig_window.ModelConfig(frame,
                                                 gpa_settings=ga_settings)

pheno_header_items = run_get_pheno_header(new_modelconfig, ga_settings)
print('Output:\n{0}'.format(pheno_header_items))

pheno_dct = run_read_phenotypic(new_modelconfig, ga_settings)
print('Output:\n{0}'.format(pheno_dct))