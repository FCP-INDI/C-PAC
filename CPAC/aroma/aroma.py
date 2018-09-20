from nipype.interfaces import afni 
from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
#from .aroma import ICA_AROMA
import os
from nipype.interfaces.fsl.aroma import ICA_AROMA
#from CPAC.aroma import (ICA_AROMA, ICA_AROMA_functions) 



def create_aroma(wf_name='create_aroma'):
    preproc = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['denoise_file','mat_file','out_dir','fnirt_warp_file']),name='inputspec')

    inputNode_params = pe.Node(util.IdentityInterface(fields=['denoise_type','TR','dim', 'warp_file_boolean']),name='params')

    outputNode = pe.Node(util.IdentityInterface(fields=['aggr_denoised_file','nonaggr_denoised_file']),name = 'outputspec')
 
   
    par_mcflirt = pe.Node(interface = fsl.MCFLIRT(),name='par_mcflirt')
    par_mcflirt.inputs.save_plots = True
    preproc.connect(inputNode,'denoise_file',par_mcflirt,'in_file')
    preproc.connect(par_mcflirt,'par_file',outputNode,'par_file')
	
	 
    bet_aroma = pe.Node(interface=fsl.BET(),name='bet_aroma')
    bet_aroma.inputs.frac=0.3
    bet_aroma.inputs.mask=True
    preproc.connect(inputNode,'denoise_file',bet_aroma,'in_file')
    preproc.connect(bet_aroma,'mask_file',outputNode,'mask_aroma')
    
    
    #wf.connect(inputnode, 'underlay', resample_u, 'file_')
    #wf.connect(resample_u, 'new_fname', outputnode,'resampled_underlay')
    
    aroma_wf = pe.Node(ICA_AROMA(),name='aroma_wf')
    preproc.connect(inputNode,'out_dir',aroma_wf,'out_dir')
    preproc.connect(inputNode,'denoise_file',aroma_wf,'in_file')
    preproc.connect(inputNode,'mat_file',aroma_wf,'mat_file')
    preproc.connect(inputNode,'fnirt_warp_file',aroma_wf,'fnirt_warp_file')
    preproc.connect(par_mcflirt,'par_file',aroma_wf,'motion_parameters')
    preproc.connect(bet_aroma,'mask_file',aroma_wf,'mask')
    preproc.connect(inputNode_params,'denoise_type',aroma_wf,'denoise_type')
    preproc.connect(inputNode_params,'TR',aroma_wf,'TR')
    preproc.connect(inputNode_params,'dim',aroma_wf,'dim')
   # preproc.connect(aroma_wf,'out_dir',outputNode,'aroma_dir')
    preproc.connect(aroma_wf,'aggr_denoised_file',outputNode,'aggr_denoised_file')
    preproc.connect(aroma_wf,'nonaggr_denoised_file',outputNode,'nonaggr_denoised_file')
	
    return preproc
