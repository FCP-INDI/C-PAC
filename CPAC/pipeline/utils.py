# Copyright (C) 2021-2022  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""C-PAC pipeline engine utilities"""
from typing import Union


def source_set(sources: Union[str, list, set]) -> set:
    """Given a CpacProvenance, return a set of {resource}:{source} strings
  
    Parameters
    ----------
    sources: str, list, or set

    Returns
    -------
    set

    Examples
    --------
    >>> source_set([[[['bold:func_ingress',
    ...       'desc-preproc_bold:func_reorient',
    ...       'desc-preproc_bold:func_truncate'],
    ...      ['TR:func_metadata_ingress'],
    ...      ['tpattern:func_metadata_ingress'],
    ...      'desc-preproc_bold:func_slice_time'],
    ...     [['bold:func_ingress',
    ...       'desc-preproc_bold:func_reorient',
    ...       'desc-preproc_bold:func_truncate'],
    ...      ['bold:func_ingress', 'desc-reorient_bold:func_reorient'],
    ...      'motion-basefile:get_motion_ref_fmriprep_reference'],
    ...     'desc-preproc_bold:motion_correction_only_mcflirt'],
    ...         [[['bold:func_ingress',
    ...           'desc-preproc_bold:func_reorient',
    ...           'desc-preproc_bold:func_truncate'],
    ...      ['bold:func_ingress', 'desc-reorient_bold:func_reorient'],
    ...      'motion-basefile:get_motion_ref_fmriprep_reference'],
    ...     [[['bold:func_ingress',
    ...        'desc-preproc_bold:func_reorient',
    ...        'desc-preproc_bold:func_truncate'],
    ...       ['TR:func_metadata_ingress'],
    ...       ['tpattern:func_metadata_ingress'],
    ...       'desc-preproc_bold:func_slice_time'],
    ...      [['bold:func_ingress',
    ...        'desc-preproc_bold:func_reorient',
    ...        'desc-preproc_bold:func_truncate'],
    ...       ['bold:func_ingress', 'desc-reorient_bold:func_reorient'],
    ...       'motion-basefile:get_motion_ref_fmriprep_reference'],
    ...      'desc-preproc_bold:motion_correction_only_mcflirt'],
    ...     ['FSL-AFNI-bold-ref:template_resample'],
    ...     ['FSL-AFNI-brain-mask:template_resample'],
    ...     ['FSL-AFNI-brain-probseg:template_resample'],
    ...     'space-bold_desc-brain_mask:bold_mask_fsl_afni'],
    ...     'desc-preproc_bold:bold_masking']) == set({
    ...     'FSL-AFNI-bold-ref:template_resample',
    ...     'FSL-AFNI-brain-mask:template_resample',
    ...     'FSL-AFNI-brain-probseg:template_resample',
    ...     'TR:func_metadata_ingress',
    ...     'bold:func_ingress',
    ...     'desc-preproc_bold:bold_masking',
    ...     'desc-preproc_bold:func_reorient',
    ...     'desc-preproc_bold:func_slice_time',
    ...     'desc-preproc_bold:func_truncate',
    ...     'desc-preproc_bold:motion_correction_only_mcflirt',
    ...     'desc-reorient_bold:func_reorient',
    ...     'motion-basefile:get_motion_ref_fmriprep_reference',
    ...     'space-bold_desc-brain_mask:bold_mask_fsl_afni',
    ...     'tpattern:func_metadata_ingress'})
    True
    """
    _set = set()
    if isinstance(sources, str):
        _set.add(sources)
    if isinstance(sources, (set, list)):
        for item in sources:
            _set.update(source_set(item))
    return _set
