# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Miscellaneous utilities
"""
# This functionality is adapted from poldracklab/niworkflows:
# https://github.com/poldracklab/niworkflows/blob/master/niworkflows/utils/misc.py        
# https://fmriprep.readthedocs.io/
# https://poldracklab.stanford.edu/
# We are temporarily maintaining our own copy for more granular control.


__all__ = ['get_template_specs']


def get_template_specs(in_template, template_spec=None, default_resolution=1):
    """
    Parse template specifications

    >>> get_template_specs('MNI152NLin2009cAsym', {'suffix': 'T1w'})[1]  # doctest: +SKIP
    {'resolution': 1}

    >>> get_template_specs('MNI152NLin2009cAsym', {'res': '2', 'suffix': 'T1w'})[1]  # doctest: +SKIP
    {'resolution': '2'}

    >>> get_template_specs('MNIInfant', {'res': '2', 'cohort': '10', 'suffix': 'T1w'})[1]  # doctest: +SKIP
    {'resolution': '2', 'cohort': '10'}

    >>> get_template_specs('MNI152NLin2009cAsym',
    ...                    {'suffix': 'T1w', 'cohort': 1})[1] # doctest: +IGNORE_EXCEPTION_DETAIL +SKIP
    Traceback (most recent call last):
    RuntimeError:
    ...

    >>> get_template_specs('MNI152NLin2009cAsym',
    ...                    {'suffix': 'T1w', 'res': '1|2'})[1] # doctest: +IGNORE_EXCEPTION_DETAIL +SKIP
    Traceback (most recent call last):
    RuntimeError:
    ...

    """
    # from templateflow.api import get as get_template
    # Massage spec (start creating if None)
    template_spec = template_spec or {}
    template_spec['desc'] = template_spec.get('desc', None)
    template_spec['atlas'] = template_spec.get('atlas', None)
    template_spec['resolution'] = template_spec.pop(
        'res', template_spec.get('resolution', default_resolution))

    common_spec = {'resolution': template_spec['resolution']}
    if 'cohort' in template_spec:
        common_spec['cohort'] = template_spec['cohort']

if __name__ == '__main__':
    pass
