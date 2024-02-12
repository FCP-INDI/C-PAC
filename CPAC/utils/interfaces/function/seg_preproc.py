"""Function interfaces for seg_preproc"""
from nipype.interfaces import utility as util


def pick_tissue_from_labels_file_interface(input_names=None):
    """Function to create a Function interface for
    CPAC.seg_preproc.utils.pick_tissue_from_labels_file

    Parameters
    ----------
    input_names : list, optional

    Returns
    -------
    nipype.interfaces.base.core.Interface
    """
    # pylint: disable=import-outside-toplevel
    from CPAC.seg_preproc.utils import pick_tissue_from_labels_file

    if input_names is None:
        input_names = ['multiatlas_Labels', 'csf_label', 'gm_label',
                       'wm_label']
    return util.Function(
        input_names=input_names,
        output_names=['csf_mask', 'gm_mask', 'wm_mask'],
        function=pick_tissue_from_labels_file)
