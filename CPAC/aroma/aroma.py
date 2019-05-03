from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.fsl.aroma import ICA_AROMA
 

def create_aroma(tr=None, wf_name='create_aroma'):
    """ICA-AROMA takes in a functional file, along with the movement
    parameters, warp file and a mat file to denoise the artifacts using
    ICA-based methods and produces an output directory which contains the
    denoised file, classification features, and a melodic.ica directory.

    Order of commands and inputs:

    -- FSL-McFlirt: in_file : denoise_file
                    out_file : par_file
                    Parameters :
                    - save plots

    -- FSL-BET:     Brain extraction of the denoised file to provide as the
                    mask file for the ICA-AROMA interface
                    in_file: denoise_file
                    out_file: mask_file
                    Parameters: 
                     -f : 0.3

    --ICA-AROMA  :  Takes in the denoise_file, par_file (FSL-McFlirt),
                    mask_file (FSL-BET), fnirt_file, mat_file
                    in_file : denoise_file
                              mat_file
                              fnirt_warp_file
                              out_dir
                    out_file: out_dir
                    Parameters :
                     - denoise_type : 'aggr', 'nonaggr'
                     - TR : s
                     - dim:
     """

    preproc = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['denoise_file',
                                                       'mat_file',
                                                       'fnirt_warp_file']),
                        name='inputspec')

    inputNode_params = pe.Node(util.IdentityInterface(fields=['denoise_type',
                                                              'dim']),
                               name='params')

    outputNode = pe.Node(util.IdentityInterface(fields=['aggr_denoised_file',
                                                        'nonaggr_denoised_file']),
                         name='outputspec')

    par_mcflirt = pe.Node(interface = fsl.MCFLIRT(),name='par_mcflirt')
    par_mcflirt.inputs.save_plots = True
    preproc.connect(inputNode,'denoise_file', par_mcflirt,'in_file')
    preproc.connect(par_mcflirt,'par_file', outputNode,'par_file')

    bet_aroma = pe.Node(interface=fsl.BET(),name='bet_aroma')
    bet_aroma.inputs.frac = 0.3
    bet_aroma.inputs.mask = True
    preproc.connect(inputNode,'denoise_file', bet_aroma,'in_file')
    preproc.connect(bet_aroma,'mask_file', outputNode,'mask_aroma')
    
    aroma = pe.Node(ICA_AROMA(), name='aroma_wf')
    aroma.inputs.out_dir = '.'
    if tr:
        aroma.inputs.TR = tr

    preproc.connect(inputNode,'denoise_file', aroma,'in_file')
    preproc.connect(inputNode,'mat_file', aroma,'mat_file')
    preproc.connect(inputNode,'fnirt_warp_file', aroma,'fnirt_warp_file')
    preproc.connect(par_mcflirt,'par_file', aroma,'motion_parameters')
    preproc.connect(bet_aroma,'mask_file', aroma,'mask')
    preproc.connect(inputNode_params,'denoise_type', aroma,'denoise_type')
    preproc.connect(inputNode_params,'dim', aroma,'dim')
    preproc.connect(aroma,'nonaggr_denoised_file', outputNode,'nonaggr_denoised_file')
    preproc.connect(aroma,'aggr_denoised_file', outputNode,'aggr_denoised_file')
	
    return preproc
