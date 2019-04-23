import yaml
from CPAC.nuisance import NuisanceRegressor


def test_representations():

    assert NuisanceRegressor.encode({
        'CerebrospinalFluid': {
            'erode_mask': True,
            'extraction_resolution': 2,
            'summary': {
                'method': 'PC',
                'components': 5,
            }
        }
    }) == 'CSF-2mmE-PC5'

    selector_test = yaml.load("""

        tCompCor:
            summary:
                method: PC
                components: 5
            threshold: 1.5SD
            by_slice: true
        
        aCompCor:
            summary:
                method: PC
                components: 5
            tissues:
                - WhiteMatter
                - CerebrospinalFluid
            extraction_resolution: 2
        
    """)

    NuisanceRegressor.encode({
        'tCompCor': selector_test['tCompCor']
    }) == 'tC-S1.5SD-PC5'

    assert NuisanceRegressor.encode({
        'aCompCor': selector_test['aCompCor']
    }) == 'aC-WM+CSF-PC5'

