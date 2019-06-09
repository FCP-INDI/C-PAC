import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.afni import preprocess

from CPAC.func_preproc.utils import add_afni_prefix


def slice_timing_wf(name='slice_timing'):

    # allocate a workflow object
    wf = pe.Workflow(name=name)

    # configure the workflow's input spec
    inputNode = pe.Node(util.IdentityInterface(fields=['func_ts',
                                                       'tr',
                                                       'tpattern']),
                        name='inputspec')

    # configure the workflow's output spec
    outputNode = pe.Node(util.IdentityInterface(fields=['edited_func']),
                         name='outputspec')

    # create TShift AFNI node
    func_slice_timing_correction = pe.Node(interface=preprocess.TShift(),
                                           name='slice_timing')
    func_slice_timing_correction.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(inputNode, 'func_ts', func_slice_timing_correction, 'in_file')
    wf.connect(inputNode, 'tr', func_slice_timing_correction, 'tr')

    # if not "Use NIFTI Header" in c.slice_timing_pattern:

    # add the @ prefix to the tpattern file going into
    # AFNI 3dTshift - needed this so the tpattern file
    # output from get_scan_params would be tied downstream
    # via a connection (to avoid poofing)
    add_prefix = pe.Node(util.Function(input_names=['tpattern'],
                                       output_names=['afni_prefix'],
                                       function=add_afni_prefix),
                         name='slice_timing_add_afni_prefix')
    wf.connect(inputNode, 'tpattern', add_prefix, 'tpattern')
    wf.connect(add_prefix, 'afni_prefix',
               func_slice_timing_correction, 'tpattern')

    wf.connect(func_slice_timing_correction, 'out_file',
               outputNode, 'slice_time_corrected')

    return wf