# -*- coding: utf-8 -*-

import os
import nipype.pipeline.engine as pe
from nipype.interfaces.afni import preprocess
import nipype.interfaces.utility as util
from CPAC.alff.utils import get_opt_string


def create_alff(wf_name='alff_workflow'):
    """
    Calculate Amplitude of low frequency oscillations (ALFF) and fractional ALFF maps

    Parameters
    ----------
    wf_name : string
        Workflow name

    Returns
    -------
    alff_workflow : workflow object
        ALFF workflow

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/alff/alff.py>`_

    Workflow Inputs::

        hp_input.hp : list of float
            high pass frequencies

        lp_input.lp : list of float
            low pass frequencies

        inputspec.rest_res : string
            Path to existing Nifti file. Nuisance signal regressed functional image.

        inputspec.rest_mask : string
            Path to existing Nifti file. A mask volume(derived by dilating the motion corrected functional volume) in native space


    Workflow Outputs::

        outputspec.alff_img : string
            Path to Nifti file. Image containing the sum of the amplitudes in the low frequency band

        outputspec.falff_img : string
            Path to Nifti file. Image containing the sum of the amplitudes in the low frequency band divided by the amplitude of the total frequency

        outputspec.alff_Z_img : string
            Path to Nifti file. Image containing Normalized ALFF Z scores across full brain in native space

        outputspec.falff_Z_img : string
            Path to Nifti file. Image containing Normalized fALFF Z scores across full brain in native space


    Order of Commands:

    - Filter the input file rest file( slice-time, motion corrected and nuisance regressed) ::
        3dBandpass -prefix residual_filtered.nii.gz
                    0.009 0.08 residual.nii.gz

    - Calculate ALFF by taking the standard deviation of the filtered file ::
        3dTstat -stdev
                -mask rest_mask.nii.gz
                -prefix residual_filtered_3dT.nii.gz
                residual_filtered.nii.gz

    - Calculate the standard deviation of the unfiltered file ::
        3dTstat -stdev
                -mask rest_mask.nii.gz
                -prefix residual_3dT.nii.gz
                residual.nii.gz

    - Calculate fALFF ::
        3dcalc -a rest_mask.nii.gz
               -b residual_filtered_3dT.nii.gz
               -c residual_3dT.nii.gz
               -expr '(1.0*bool(a))*((1.0*b)/(1.0*c))' -float

    - Normalize ALFF/fALFF to Z-score across full brain ::

        fslstats
        ALFF.nii.gz
        -k rest_mask.nii.gz
        -m > mean_ALFF.txt ; mean=$( cat mean_ALFF.txt )

        fslstats
        ALFF.nii.gz
        -k rest_mask.nii.gz
        -s > std_ALFF.txt ; std=$( cat std_ALFF.txt )

        fslmaths
        ALFF.nii.gz
        -sub ${mean}
        -div ${std}
        -mas rest_mask.nii.gz ALFF_Z.nii.gz

        fslstats
        fALFF.nii.gz
        -k rest_mask.nii.gz
        -m > mean_fALFF.txt ; mean=$( cat mean_fALFF.txt )

        fslstats
        fALFF.nii.gz
        -k rest_mask.nii.gz
        -s > std_fALFF.txt
        std=$( cat std_fALFF.txt )

        fslmaths
        fALFF.nii.gz
        -sub ${mean}
        -div ${std}
        -mas rest_mask.nii.gz
        fALFF_Z.nii.gz

    .. exec::
        from CPAC.alff import create_alff
        wf = create_alff()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/alff.dot'
        )

    High Level Workflow Graph:

    .. image:: ../../images/generated/alff.png
        :width: 500

    Detailed Workflow Graph:

    .. image:: ../../images/generated/alff_detailed.png
        :width: 500


    References
    ----------

    .. [1] Zou, Q.-H., Zhu, C.-Z., Yang, Y., Zuo, X.-N., Long, X.-Y., Cao, Q.-J., Wang, Y.-F., et al. (2008). An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: fractional ALFF. Journal of neuroscience methods, 172(1), 137-41. doi:10.10

    Examples
    --------

    >>> alff_w = create_alff()
    >>> alff_w.inputs.hp_input.hp = [0.01]
    >>> alff_w.inputs.lp_input.lp = [0.1]
    >>> alff_w.get_node('hp_input').iterables = ('hp', [0.01])
    >>> alff_w.get_node('lp_input').iterables = ('lp', [0.1])
    >>> alff_w.inputs.inputspec.rest_res = '/home/data/subject/func/rest_bandpassed.nii.gz'
    >>> alff_w.inputs.inputspec.rest_mask= '/home/data/subject/func/rest_mask.nii.gz'
    >>> alff_w.run() # doctest: +SKIP
    """

    wf = pe.Workflow(name=wf_name)
    input_node = pe.Node(util.IdentityInterface(fields=['rest_res',
                                                        'rest_mask']),
                         name='inputspec')

    input_node_hp = pe.Node(util.IdentityInterface(fields=['hp']),
                            name='hp_input')

    input_node_lp = pe.Node(util.IdentityInterface(fields=['lp']),
                            name='lp_input')

    output_node = pe.Node(util.IdentityInterface(fields=['alff_img',
                                                         'falff_img']),
                          name='outputspec')

    # filtering
    bandpass = pe.Node(interface=preprocess.Bandpass(),
                       name='bandpass_filtering')
    bandpass.inputs.outputtype = 'NIFTI_GZ'
    bandpass.inputs.out_file = os.path.join(os.path.curdir,
                                            'residual_filtered.nii.gz')

    wf.connect(input_node_hp, 'hp', bandpass, 'highpass')
    wf.connect(input_node_lp, 'lp', bandpass, 'lowpass')
    wf.connect(input_node, 'rest_res', bandpass, 'in_file')

    get_option_string = pe.Node(util.Function(input_names=['mask'],
                                              output_names=['option_string'],
                                              function=get_opt_string),
                                name='get_option_string')

    wf.connect(input_node, 'rest_mask', get_option_string, 'mask')

    # standard deviation over frequency
    try:
        from nipype.interfaces.afni import utils as afni_utils
        stddev_filtered = pe.Node(interface=afni_utils.TStat(),
                                  name='stddev_filtered')
    except ImportError:
        stddev_filtered = pe.Node(interface=preprocess.TStat(),
                                  name='stddev_filtered')

    stddev_filtered.inputs.outputtype = 'NIFTI_GZ'
    stddev_filtered.inputs.out_file = os.path.join(os.path.curdir,
                                                   'alff.nii.gz')
                                                
    wf.connect(bandpass, 'out_file', stddev_filtered, 'in_file')
    wf.connect(get_option_string, 'option_string', stddev_filtered, 'options')

    wf.connect(stddev_filtered, 'out_file', output_node, 'alff_img')

    # standard deviation of the unfiltered nuisance corrected image
    try:
        stddev_unfiltered = pe.Node(interface=afni_utils.TStat(),
                                    name='stddev_unfiltered')
    except UnboundLocalError:
        stddev_unfiltered = pe.Node(interface=preprocess.TStat(),
                                    name='stddev_unfiltered')

    stddev_unfiltered.inputs.outputtype = 'NIFTI_GZ'
    stddev_unfiltered.inputs.out_file = os.path.join(os.path.curdir,
                                                     'residual_3dT.nii.gz')

    wf.connect(input_node, 'rest_res', stddev_unfiltered, 'in_file')
    wf.connect(get_option_string, 'option_string', stddev_unfiltered,
               'options')

    # falff calculations
    try:
        falff = pe.Node(interface=afni_utils.Calc(), name='falff')
    except UnboundLocalError:
        falff = pe.Node(interface=preprocess.Calc(), name='falff')

    falff.inputs.args = '-float'
    falff.inputs.expr = '(1.0*bool(a))*((1.0*b)/(1.0*c))'
    falff.inputs.outputtype = 'NIFTI_GZ'
    falff.inputs.out_file = os.path.join(os.path.curdir, 'falff.nii.gz')

    wf.connect(input_node, 'rest_mask', falff, 'in_file_a')
    wf.connect(stddev_filtered, 'out_file', falff, 'in_file_b')
    wf.connect(stddev_unfiltered, 'out_file', falff, 'in_file_c')

    wf.connect(falff, 'out_file', output_node, 'falff_img')

    return wf


def run_alff(input_fmri, func_brain_mask, hp=0.01, lp=0.1, out_dir=None,
             run=True):
    """Runner function for the create_alff workflow builder."""

    import os
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = 'alff'

    workflow = pe.Workflow(name='{0}_workflow'.format(output))

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    num_cores_per_subject = 1

    alff = create_alff('alff_falff')

    alff.inputs.inputspec.rest_res = os.path.abspath(input_fmri)
    alff.inputs.inputspec.rest_mask = os.path.abspath(func_brain_mask)
    alff.inputs.hp_input.hp = float(hp)
    alff.inputs.lp_input.lp = float(lp)

    ds = pe.Node(nio.DataSink(), name='datasink_{0}'.format(output))
    ds.inputs.base_directory = workflow_dir

    workflow.connect(alff, 'outputspec.alff_img', ds, 'alff')
    workflow.connect(alff, 'outputspec.falff_img', ds, 'falff')

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))
        return outpath

    else:
        return workflow, workflow.base_dir
