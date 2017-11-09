
import nipype.pipeline.engine as pe
from nipype.interfaces import afni,fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype
import os,glob,sys

"""
  Fieldmap correction takes in an input magnitude image which is Skull Stripped (Tight).
  The magnitude images are obtained from each echo series. It also requires a phase image
  as an input, the phase image is a subtraction of the two phase images from each echo.
  
  Order of commands and inputs:
    -- Extraction of brain volumes: using first three columes from the fMRI
                                    
    -- SkullStrip: SkullStrip (or BET) is used to strip the non-brain (tissue) regions
                   from the fMRI
                   
    -- PrepareFieldMap: Preparing the fieldmap. 
                           For more details, check:https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide
                           in_file = phase_file, magnitude_file
                           deltaTE = Delay in the time taken by the scanners to record the two echos.
                                     (default == 2.46 ms)
                           scanner = SIEMENS, default.
    -- FUGUE           : One of the steps in EPI-DistCorrection toolbox, it unwarps the fieldmaps 
                            in_file = field map which is a 4D image (containing 2 unwarpped image)
                            mask_file = epi_mask.nii
                            dwell_to_asymm ratio = (0.77e-3 * 3)/(2.46e-3)
  """                          
                           
inputNode = pe.Node(util.IdentityInterface(fields=['distcorr']),
                        name='inputspec')

outputNode = pe.Node(util.IdentityInterface(fields=['roi_file','fieldmap','epireg']),name='outputspec')
preproc = pe.Workflow(name='preprocflow')
## Specify commands to be run

# Extract first three volumes from fmri
fslroi = pe.Node(interface=fsl.ExtractROI(),name='fslroi')
fslroi.inputs.t_min=0
fslroi.inputs.t_size=3

preproc.connect(inputNode,'distcorr',fslroi,'in_file') 

preproc.connect(fslroi,'roi_file',outputNode,'fslroi_file')
# Skullstrip

skullstrip = pe.Node(interface=afni.preprocess.SkullStrip(),name='skullstrip')
skullstrip.inputs.outputtype='NIFTI_GZ'

preproc.connect(inputNode,'fmap_mag',skullstrip,'in_file')


#SkullStrip the anatomy file
skullstrip_anat = pe.Node(interface=afni.preprocess.SkullStrip(),name='skullstrip_anat')
skullstrip_anat.inputs.outputtype='NIFTI_GZ'

preproc.connect(inputNode,'anat',skullstrip_anat,'in_file')
preproc.connect(skullstrip_anat,'out_file',outputNode,'stripped_anat')

# Prepare Fieldmap
prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(),name='prepare')
prepare.inputs.output_type = "NIFTI_GZ"
prepare.inputs.delta_TE = 2.46

preproc.connect(inputNode,'fmap_pha',prepare,'in_phase')
preproc.connect(inputNode,'fmap_mag',prepare,'in_magnitude')
preproc.connect(prepare,'out_fieldmap',outputNode,'fieldmap')

#perform fugue
fugue = pe.Node(interface = fsl.FUGUE(),name='fugue')
fugue.unwarp_direction='y'
fugue.save_unmasked_fmap=True
fugue.shift_out_file = 'opname_shiftout_file'
fugue.dwell_to_asym_ratio=(0.77e-3*3)/(2.46e-3)

preproc.connect(inputNode,'mask',fugue,'mask_file')
preproc.connect(prepare,'out_fieldmap',fugue,'in_file')


# Co-Register EPI and Correct field inhomogeniety distortions
epireg = pe.Node(interface=fsl.epi.EpiReg(), name='epireg')
epireg.inputs.echospacing=0.00046
epireg.inputs.pedir='-y'
epireg.inputs.output_type='NIFTI_GZ'
preproc.connect(skullstrip_anat,'out_file', epireg,'t1_brain')
preproc.connect(skullstrip,'out_file',epireg,'fmapmagbrain')
preproc.connect(prepare,'out_fieldmap',epireg,'fmap')
preproc.connect(epireg,'out_file',outputNode,'epireg')


## Connect nodes and run workflow


#preproc.connect(skullstrip,'out_file',outputNode,'mag_ss')                 
#preproc.connect(skullstrip_anat,'out_file',outputNode,'anat_ss')
#preproc.connect(skullstrip,'out_file',prepare,'in_magnitude')
#preproc.connect(prepare,'out_fieldmap',outputNode,'fmap_prep')
#preproc.connect(skullstrip_anat,'out_file',epireg,'t1_brain')
#preproc.connect(skullstrip,'out_file',epireg,'fmapmagbrain')
#preproc.connect(prepare,'out_fieldmap',epireg,'fmap')
#preproc.connect(epireg,'out_file',datasink,'epireg')


#preproc.run('MultiProc',plugin_args={'n_procs':4})
