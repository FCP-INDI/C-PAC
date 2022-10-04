import os
import numpy as np
import subprocess
import re
import scipy.io
import numpy as np
import re

import os
import nipype
from nipype.pipeline import engine as pe
from nipype.interfaces import afni
from nipype.interfaces import ants
from nipype.interfaces import fsl
from nipype.interfaces import freesurfer
from nipype.interfaces import utility
import nipype.interfaces.utility as util
from nipype.interfaces.fsl import utils as fsl_utils
from nipype.interfaces.fsl import maths as fsl_maths
from nipype.interfaces.io import DataSink
from nipype.interfaces.base import CommandLine

InflateExtraScale = 1

#wb_command -volume-fill-holes "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
def run_volume_fill_holes(in_file):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), 'T1wImageBrainMask_1mm.nii.gz')
        cmd=  ['wb_command', '-volume-fill-holes', in_file , out_file]
        retcode = subprocess.check_output(cmd)
        return out_file
        
def run_surface_apply_affine(in_file,outfile_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), outfile_name)
        cmd=  ['wb_command','-surface-apply-affine', in_file,'/data3/cnl/tgeorge/pfreesurfer-bash-runs/python-c_ras.mat', out_file]
        retcode = subprocess.check_output(cmd)
        return out_file


def surface_apply_warpfield(in_file,outfile_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), outfile_name)
        cmd=  ['wb_command','-surface-apply-warpfield', in_file, '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/merge_template_to_t1_xfms_82/einv1_merged.nii.gz',out_file, '-fnirt', '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/merge_t1_to_template_xfms_82/e1_merged.nii.gz']
        retcode = subprocess.check_output(cmd)
        return out_file

          
def run_volume_label_1(in_file,file_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), file_name)
        cmd=  ['wb_command', '-volume-label-import', in_file, "/data3/cnl/tgeorge/pfreesurfer-bash-runs/FreeSurferAllLut.txt", out_file, "-drop-unused-labels"]
        retcode = subprocess.check_output(cmd)
        return out_file

def run_volume_label_2(in_file):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), '2.nii.gz')
        cmd=  ['wb_command', '-volume-label-import', in_file,'/data3/cnl/tgeorge/pfreesurfer-bash-runs/FreeSurferSubcorticalLabelTableLut.txt' , out_file, "-discard-others"]
        retcode = subprocess.check_output(cmd)
        return out_file
        

     
def run_surface_avg(in_file_1,in_file_2):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), 'midthickness.native.surf.gii')
        cmd=  ['wb_command', '-surface-average', out_file ,'-surf', in_file_1,'-surf', in_file_2]
        retcode = subprocess.check_output(cmd)
        return out_file

def run_generate_sur_inflated(in_file_1, in_file_2):

        import os
        import subprocess
        out_file_1 = os.path.join(os.getcwd(), 'inflated.native.surf.gii')
        out_file_2 = os.path.join(os.getcwd(), 'very_inflated.native.surf.gii')
        cmd=  ['wb_command','-surface-generate-inflated', in_file_1, out_file_1, out_file_2, '-iterations-scale', in_file_2]
        retcode = subprocess.check_output(cmd)
        return out_file_1, out_file_2

def add_spec_file(spec_file,in_file,structure):

        import os
        import subprocess
        cmd=  ['wb_command','-add-to-spec-file',spec_file,structure, in_file]
        retcode = subprocess.check_output(cmd)
        return spec_file

def get_nativeinflationscale(InflateExtraScale, in_val):

        import os
        import subprocess
        
        NativeInflationScale = os.path.join(os.getcwd(),'NativeInflationScale')
        NativeInflationScale = format(InflateExtraScale*0.75*(float(in_val)/32492),'.4f')
        
        return NativeInflationScale

def get_native_verts(in_file):

     import os
     import subprocess
     import numpy as np

     NativeVerts = os.path.join(os.getcwd(), 'NativeVerts')
     out_file = os.path.join(os.getcwd(), 'verts_file_info.txt')

     cmd = ['wb_command', '-file-information', in_file]
     file_info = subprocess.check_output(cmd)



     with open('verts_file_info.txt', 'w') as f:
            for v in file_info:

                f.write(str(file_info))



    

     for line in open(out_file, 'r'):


         NativeVerts = line.strip().split('Number of Vertices:')
         NativeVerts = NativeVerts[1].strip().split('\\n')
         NativeVerts = NativeVerts[0]

    
   
     return NativeVerts



def run_set_structure_1(in_file,structure):

        import os
        import subprocess
        cmd=  ['wb_command','-set-structure', in_file, structure ]
        retcode = subprocess.check_output(cmd)
        return in_file
        
def run_set_structure_2(in_file,structure):

        import os
        import subprocess
        cmd=  ['wb_command','-set-structure', in_file, structure, '-surface-type', 'SPHERICAL' ]
        retcode = subprocess.check_output(cmd)
        return in_file

        
def run_set_structure_3(in_file,structure, surface_secondary_type):

        import os
        import subprocess
        cmd=  ['wb_command','-set-structure', in_file, structure, '-surface-type', 'ANATOMICAL', '-surface-secondary-type', surface_secondary_type ]
        retcode = subprocess.check_output(cmd)
        return in_file
        



def run_set_maps(in_file, map_name):

        import os
        import subprocess
        cmd=  ['wb_command','-set-map-names', in_file, '-map', '1', map_name]
        retcode = subprocess.check_output(cmd)
        return in_file
        
def run_metric_palette(in_file):

        import os
        import subprocess
        cmd=  ['wb_command', '-metric-palette', in_file , 'MODE_AUTO_SCALE_PERCENTAGE', '-pos-percent', '2', '98', '-palette-name', 'Gray_Interp', '-disp-pos', 'true', '-disp-neg','true', '-disp-zero', 'true']
        retcode = subprocess.check_output(cmd)
        return in_file



def run_set_mapnames(in_file,map_name):

        import os
        import subprocess
        cmd=  ['wb_command', '-set-map-names', in_file ,'-map', '1', map_name]
        retcode = subprocess.check_output(cmd)
        return in_file
        
def run_set_giftilabels(in_file,prefix,file_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), file_name)
        cmd=  ['wb_command', '-gifti-label-add-prefix',in_file, prefix, out_file]
        retcode = subprocess.check_output(cmd)
        return out_file

def run_metric_maths_1(in_file, file_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), file_name)
        cmd=  ['wb_command','-metric-math', "thickness > 0" , out_file, '-var',  'thickness' , in_file]
        retcode = subprocess.check_output(cmd)
        return out_file



def run_metric_maths_2(in_file,file_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), file_name)
        cmd=  ['wb_command', '-metric-math', "abs(thickness)", out_file, '-var', 'thickness', in_file]
        retcode = subprocess.check_output(cmd)
        return out_file
        
def run_metric_maths_3(in_file,outfile_name):

        import os
        import subprocess
        cmd=  ['wb_command','-metric-math', "var * -1", in_file, '-var' , 'var', in_file ]
        retcode = subprocess.check_output(cmd)
        return in_file
                

        

def run_metric_palette(in_file):

        import os
        import subprocess
        cmd=  ['wb_command', '-metric-palette', in_file, 'MODE_AUTO_SCALE_PERCENTAGE', '-pos-percent', '4', '96', '-interpolate', 'true', '-palette-name', 'videen_style', '-disp-pos', 'true', '-disp-neg','false', '-disp-zero', 'false']
        retcode = subprocess.check_output(cmd)
        return in_file
        

        
def run_metric_fillholes(surface_file,metric_infile,outfile_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), outfile_name)
        cmd= ["wb_command", "-metric-fill-holes",surface_file, metric_infile, out_file]
        retcode = subprocess.check_output(cmd)
        return out_file
        
def run_metric_removeislands(surface_file,metric_infile,outfile_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), outfile_name)
        cmd=  ["wb_command", "-metric-remove-islands", surface_file, metric_infile, out_file]
        retcode = subprocess.check_output(cmd)
        return out_file
        
def run_setmaps(in_file, map_name):

        import os
        import subprocess
        cmd=  ["wb_command", "-set-map-names", in_file, "-map", "1",map_name ]
        retcode = subprocess.check_output(cmd)
        return in_file
        
def run_metric_dilate(in_file,surface_file,outfile_name):

        import os
        import subprocess
        out_file = os.path.join(os.getcwd(), outfile_name)
        cmd=  ["wb_command", "-metric-dilate", in_file, surface_file, "10", out_file, "-nearest"]
        retcode = subprocess.check_output(cmd)
        return out_file
        
def create_res_mat(mri_info):


     cmd = ['mri_info', mri_info]
     mri_info_file = subprocess.check_output(cmd)



     with open('mri_info.txt', 'w') as f:
            for v in mri_info_file:
				
                f.write(str(mri_info_file))
                
            


     out_file = os.path.join(os.getcwd(), 'mri_info.txt')

     for line in open(out_file, 'r'):
    
       
         cr_val = line.strip().split('c_r =')
         cr = cr_val[1].strip().split('\\n')
        
         ca_val = line.strip().split('c_a =')
         ca = ca_val[1].strip().split('\\n')
     
         cs_val = line.strip().split('c_s =')
         cs = cs_val[1].strip().split('\\n')
    
    

     cr = re.sub(r"[\n\t\s]*", "", cr[0])
     cr = float(cr)

     ca = re.sub(r"[\n\t\s]*", "", ca[0])
     ca = float(ca)

     cs = re.sub(r"[\n\t\s]*", "", cs[0])
     cs = float(cs)


     cr_matrix = np.array([1, 0, 0, cr])
     ca_matrix = np.array([0, 1, 0, ca])
     cs_matrix = np.array([0, 0, 1, cs])
     id_mastrix = np.array([0, 0, 0, 1])

     final_mat = np.concatenate((cr_matrix, ca_matrix, cs_matrix, id_mastrix), axis=0)
     final_mat = final_mat.reshape(4,4)
     file_path = '/data3/cnl/tgeorge/pfreesurfer-bash-runs/python-c_ras.mat'

     np.savetxt('python-c_ras.mat', final_mat,fmt='%1.0f')
create_res_mat("/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/mri/wmparc.mgz")
     
def run_workflow():

    wf = pe.Workflow(name="dcan_postfreesurfer_combined")
    wf_dir = os.path.join(os.getcwd(), "dcan_postfreesurfer_wf")
    wf.base_dir = wf_dir
    wf.config["execution"] = {"hash_method": "timestamp", "crashdump_dir": wf_dir}


    #mri_convert -rt nearest -rl "$T1wFolder"/"$T1wImage".nii.gz "$FreeSurferFolder"/mri/"$Image".mgz "$T1wFolder"/"$Image"_1mm.nii.gz

    mri_convert_wmparc = pe.Node(interface=freesurfer.preprocess.MRIConvert(),name='mri_convert_wmparc')
    mri_convert_wmparc.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/mri/wmparc.mgz"
    mri_convert_wmparc.out_file = "wmparc_1mm.nii.gz"
    mri_convert_wmparc.inputs.resample_type = 'nearest'
    mri_convert_wmparc.inputs.reslice_like = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/fast_bias_field_correction_52/get_anat_restore/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths.nii.gz"


    #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$Image".nii.gz

    apply_warp_1_wmparc = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_1_wmparc')
    apply_warp_1_wmparc.inputs.interp = 'nn'
    apply_warp_1_wmparc.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_1_wmparc.out_file = "wmparc.nii.gz"
    apply_warp_1_wmparc.relwarp = 'TRUE'
    apply_warp_1_wmparc.premat = '/usr/local/fsl/etc/flirtsch/ident.mat'
    wf.connect(mri_convert_wmparc ,'out_file',apply_warp_1_wmparc,'in_file')

    #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$Image".nii.gz
    apply_warp_2_wmparc = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_2_wmparc')
    apply_warp_2_wmparc.inputs.interp = 'nn'
    apply_warp_2_wmparc.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_2_wmparc.out_file = "wmparc.nii.gz"
    apply_warp_2_wmparc.relwarp = True
    apply_warp_2_wmparc.inputs.field_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/merge_t1_to_template_xfms_82/e1_merged.nii.gz"
    wf.connect(mri_convert_wmparc ,'out_file',apply_warp_2_wmparc,'in_file')



    #wb_command -volume-label-import "$T1wFolder"/"$Image".nii.gz "$FreeSurferLabels" "$T1wFolder"/"$Image".nii.gz -drop-unused-labels
    volume_label_imports_1_wmparc = ['import os', 'import subprocess']
    volume_label_1_wmparc = pe.Node(util.Function(input_names=['in_file', 'file_name'],
                                               output_names=['out_file'],
                                               function=run_volume_label_1,
                                               imports=volume_label_imports_1_wmparc),
                                 name='volume_label_1_wmparc')

    wf.connect(apply_warp_1_wmparc,'out_file',volume_label_1_wmparc,'in_file')
    volume_label_1_wmparc.inputs.file_name = 'wmparc.nii.gz'



    #wb_command -volume-label-import "$AtlasSpaceFolder"/"$Image".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/"$Image".nii.gz -drop-unused-labels
    volume_label_imports_2_wmparc = ['import os', 'import subprocess']
    volume_label_2_wmparc = pe.Node(util.Function(input_names=['in_file'],
                                               output_names=['out_file'],
                                               function=run_volume_label_2,
                                               imports=volume_label_imports_2_wmparc),
                                 name='volume_label_2_wmparc')

    wf.connect(apply_warp_2_wmparc,'out_file',volume_label_2_wmparc,'in_file')
    
    ##########################################
    
    #mri_convert -rt nearest -rl "$T1wFolder"/"$T1wImage".nii.gz "$FreeSurferFolder"/mri/"$Image".mgz "$T1wFolder"/"$Image"_1mm.nii.gz

    mri_convert_aparca2009saseg = pe.Node(interface=freesurfer.preprocess.MRIConvert(),name='mri_convert_aparca2009saseg')
    mri_convert_aparca2009saseg.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/mri/aparc.a2009s+aseg.mgz"
    mri_convert_aparca2009saseg.out_file = "aparca2009saseg_1mm.nii.gz"
    mri_convert_aparca2009saseg.inputs.resample_type = 'nearest'
    mri_convert_aparca2009saseg.inputs.reslice_like = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/fast_bias_field_correction_52/get_anat_restore/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths.nii.gz"


    #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$Image".nii.gz

    apply_warp_1_aparca2009saseg = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_1_aparca2009saseg')
    apply_warp_1_aparca2009saseg.inputs.interp = 'nn'
    apply_warp_1_aparca2009saseg.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_1_aparca2009saseg.out_file = "aparc.a2009saseg.nii.gz"
    apply_warp_1_aparca2009saseg.relwarp = 'TRUE'
    apply_warp_1_aparca2009saseg.premat = '/usr/local/fsl/etc/flirtsch/ident.mat'
    wf.connect(mri_convert_aparca2009saseg ,'out_file',apply_warp_1_aparca2009saseg,'in_file')

    #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$Image".nii.gz
    apply_warp_2_aparca2009saseg = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_2_aparca2009saseg')
    apply_warp_2_aparca2009saseg.inputs.interp = 'nn'
    apply_warp_2_aparca2009saseg.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_2_aparca2009saseg.out_file = "aparca2009saseg.nii.gz"
    apply_warp_2_aparca2009saseg.relwarp = True
    apply_warp_2_aparca2009saseg.inputs.field_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/merge_t1_to_template_xfms_82/e1_merged.nii.gz"
    wf.connect(mri_convert_aparca2009saseg ,'out_file',apply_warp_2_aparca2009saseg,'in_file')



    #wb_command -volume-label-import "$T1wFolder"/"$Image".nii.gz "$FreeSurferLabels" "$T1wFolder"/"$Image".nii.gz -drop-unused-labels
    volume_label_imports_1_aparca2009saseg = ['import os', 'import subprocess']
    volume_label_1_aparca2009saseg = pe.Node(util.Function(input_names=['in_file', 'file_name'],
                                               output_names=['out_file'],
                                               function=run_volume_label_1,
                                               imports=volume_label_imports_1_aparca2009saseg),
                                 name='volume_label_1_aparca2009saseg')

    volume_label_1_aparca2009saseg.inputs.file_name = 'aparca2009saseg.nii.gz'
    wf.connect(apply_warp_1_aparca2009saseg,'out_file',volume_label_1_aparca2009saseg,'in_file')




    #wb_command -volume-label-import "$AtlasSpaceFolder"/"$Image".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/"$Image".nii.gz -drop-unused-labels
    volume_label_imports_2_aparca2009saseg = ['import os', 'import subprocess']
    volume_label_2_aparca2009saseg = pe.Node(util.Function(input_names=['in_file'],
                                               output_names=['out_file'],
                                               function=run_volume_label_2,
                                               imports=volume_label_imports_2_aparca2009saseg),
                                 name='volume_label_2_aparca2009saseg')

    wf.connect(apply_warp_2_aparca2009saseg,'out_file',volume_label_2_aparca2009saseg,'in_file')
    
    ###########
    
    #mri_convert -rt nearest -rl "$T1wFolder"/"$T1wImage".nii.gz "$FreeSurferFolder"/mri/"$Image".mgz "$T1wFolder"/"$Image"_1mm.nii.gz

    mri_convert_aparcaseg = pe.Node(interface=freesurfer.preprocess.MRIConvert(),name='mri_convert_aparcaseg')
    mri_convert_aparcaseg.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/mri/aparc+aseg.mgz"
    mri_convert_aparcaseg.out_file = "aparc+aseg_1mm.nii.gz"
    mri_convert_aparcaseg.inputs.resample_type = 'nearest'
    mri_convert_aparcaseg.inputs.reslice_like = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/fast_bias_field_correction_52/get_anat_restore/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths.nii.gz"


    #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$Image".nii.gz

    apply_warp_1_aparcaseg = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_1_aparcaseg')
    apply_warp_1_aparcaseg.inputs.interp = 'nn'
    apply_warp_1_aparcaseg.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_1_aparcaseg.out_file = "wmparc.nii.gz"
    apply_warp_1_aparcaseg.relwarp = 'TRUE'
    apply_warp_1_aparcaseg.premat = '/usr/local/fsl/etc/flirtsch/ident.mat'
    wf.connect(mri_convert_aparcaseg ,'out_file',apply_warp_1_aparcaseg,'in_file')

    #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$Image".nii.gz
    apply_warp_2_aparcaseg = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_2_aparcaseg')
    apply_warp_2_aparcaseg.inputs.interp = 'nn'
    apply_warp_2_aparcaseg.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_2_aparcaseg.out_file = "wmparc.nii.gz"
    apply_warp_2_aparcaseg.relwarp = True
    apply_warp_2_aparcaseg.inputs.field_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/merge_t1_to_template_xfms_82/e1_merged.nii.gz"
    wf.connect(mri_convert_aparcaseg ,'out_file',apply_warp_2_aparcaseg,'in_file')



    #wb_command -volume-label-import "$T1wFolder"/"$Image".nii.gz "$FreeSurferLabels" "$T1wFolder"/"$Image".nii.gz -drop-unused-labels
    volume_label_imports_1_aparcaseg = ['import os', 'import subprocess']
    volume_label_1_aparcaseg = pe.Node(util.Function(input_names=['in_file','file_name'],
                                               output_names=['out_file'],
                                               function=run_volume_label_1,
                                               imports=volume_label_imports_1_aparcaseg),
                                 name='volume_label_1_aparcaseg')

    volume_label_1_aparcaseg.inputs.file_name = 'aparcaseg.nii.gz'
    wf.connect(apply_warp_1_aparcaseg,'out_file',volume_label_1_aparcaseg,'in_file')




    #wb_command -volume-label-import "$AtlasSpaceFolder"/"$Image".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/"$Image".nii.gz -drop-unused-labels
    volume_label_imports_2_aparcaseg = ['import os', 'import subprocess']
    volume_label_2_aparcaseg = pe.Node(util.Function(input_names=['in_file'],
                                               output_names=['out_file'],
                                               function=run_volume_label_2,
                                               imports=volume_label_imports_2_aparcaseg),
                                 name='volume_label_2_aparcaseg')

    wf.connect(apply_warp_2_aparcaseg,'out_file',volume_label_2_aparcaseg,'in_file')
    
    
    #fslmaths "$T1wFolder"/wmparc_1mm.nii.gz -bin -dilD -dilD -dilD -ero -ero "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz

    fsl_maths_1 = pe.Node(interface=fsl.maths.UnaryMaths(),name='fsl_maths_1')
    fsl_maths_1.inputs.operation = "bin"
    fsl_maths_1.inputs.args = "-dilD -dilD -dilD -ero -ero"
    fsl_maths_1.out_file = "T1wImageBrainMask_1mm.nii.gz"
    wf.connect(mri_convert_wmparc, 'out_file',fsl_maths_1, 'in_file')

    #wb_command -volume-fill-holes "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
    volume_fill_holes_imports = ['import os', 'import subprocess']
    volume_fill_holes = pe.Node(util.Function(input_names=['in_file'],
                                               output_names=['out_file'],
                                               function=run_volume_fill_holes,
                                               imports=volume_fill_holes_imports),
                                 name='volume_fill-holes')

    wf.connect(fsl_maths_1,'out_file',volume_fill_holes,'in_file')

    #fslmaths "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -bin "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz

    fsl_maths_2 = pe.Node(interface=fsl.maths.UnaryMaths(),name='fsl_maths_2')
    fsl_maths_2.inputs.operation = "bin"
    fsl_maths_2.out_file = "T1wImageBrainMask_1mm.nii.gz"
    wf.connect(volume_fill_holes,'out_file',fsl_maths_2,'in_file')
    
   #applywarp --rel --interp=nn -i "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$T1wImageBrainMask".nii.gz
   
    apply_warp_1_brainmask = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_1_brainmask')
    apply_warp_1_brainmask.inputs.interp = 'nn'
    apply_warp_1_brainmask.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_1_brainmask.out_file = "T1wImageBrainMask.nii.gz"
    apply_warp_1_brainmask.relwarp = 'TRUE'
    apply_warp_1_brainmask.premat = '/usr/local/fsl/etc/flirtsch/ident.mat'
    wf.connect(fsl_maths_2 ,'out_file',apply_warp_1_brainmask,'in_file')

   
   #applywarp --rel --interp=nn -i "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$T1wImageBrainMask".nii.gz

    apply_warp_2_brainmask = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_2_brainmask')
    apply_warp_2_brainmask.inputs.interp = 'nn'
    apply_warp_2_brainmask.inputs.ref_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    apply_warp_2_brainmask.out_file = "T1wImageBrainMask.nii.gz"
    apply_warp_2_brainmask.relwarp = True
    apply_warp_2_brainmask.inputs.field_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/merge_t1_to_template_xfms_82/e1_merged.nii.gz"
    wf.connect(fsl_maths_2 ,'out_file',apply_warp_2_brainmask,'in_file')

    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$T1wFolder"/"$T1wImage".nii.gz
    
    T1w_add_spec_file_imports = ['import os', 'import subprocess']
    T1w_add_spec_file = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=T1w_add_spec_file_imports),
                                 name='T1w_add_spec_file')
    T1w_add_spec_file.inputs.spec_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/T1w_sub-0050952_ses-1.native.wb.spec"
    T1w_add_spec_file.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/fast_bias_field_correction_52/get_anat_restore/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths.nii.gz"
    T1w_add_spec_file.inputs.structure = "INVALID"
    

    

    #  wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz
    
    Atlasspace_add_spec_file_imports = ['import os', 'import subprocess']
    Atlasspace_add_spec_file = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=Atlasspace_add_spec_file_imports),
                                 name='Atlasspace_add_spec_file')
    Atlasspace_add_spec_file.inputs.spec_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/MNINonLinear_sub-0050952_ses-1.native.wb.spec"
    Atlasspace_add_spec_file.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    Atlasspace_add_spec_file.inputs.structure = "INVALID"

    
   
    #wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz

    Atlasspace_HighResMesh_add_spec_file_imports = ['import os', 'import subprocess']
    Atlasspace_HighResMesh_add_spec_file = pe.Node(util.Function(input_names=['spec_file', 'in_file', 'structure'],
                                               output_names=['out_file'],
                                               function=add_spec_file,
                                               imports=Atlasspace_HighResMesh_add_spec_file_imports),
                                 name='Atlasspace_HighResMesh_add_spec_file')
    Atlasspace_HighResMesh_add_spec_file.inputs.spec_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/sub-0050952_ses-1.164k_fs_LR.wb.spec"
    Atlasspace_HighResMesh_add_spec_file.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    Atlasspace_HighResMesh_add_spec_file.inputs.structure = 'INVALID'
    wf.add_nodes([Atlasspace_HighResMesh_add_spec_file])
    
    
    Atlasspace_LowResMesh_add_spec_file_imports = ['import os', 'import subprocess']
    Atlasspace_LowResMesh_add_spec_file = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                               output_names=['out_file'],
                                               function=add_spec_file,
                                               imports=Atlasspace_LowResMesh_add_spec_file_imports),
                                 name='Atlasspace_LowResMesh_add_spec_file')
    Atlasspace_LowResMesh_add_spec_file.inputs.spec_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/MNINonLinear_sub-0050952_ses-1.32k_fs_LR.wb.spec"
    Atlasspace_LowResMesh_add_spec_file.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/FSL-ABCD_T1_to_template_82/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths_warp.nii.gz"
    Atlasspace_LowResMesh_add_spec_file.inputs.structure = "INVALID"

    wf.add_nodes([Atlasspace_LowResMesh_add_spec_file])
   
    #wb_command -add-to-spec-file "$T1wFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec INVALID "$T1wFolder"/"$T1wImage".nii.gz
    T1w_LowResMesh_add_spec_file_imports = ['import os', 'import subprocess']
    T1w_LowResMesh_add_spec_file = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                                   output_names=['out_file'],
                                                   function=add_spec_file,
                                                   imports=T1w_LowResMesh_add_spec_file_imports),
                                     name='T1w_LowResMesh_add_spec_file')
    T1w_LowResMesh_add_spec_file.inputs.spec_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/T1w_sub-0050952_ses-1.32k_fs_LR.wb.spec"
    T1w_LowResMesh_add_spec_file.inputs.in_file = "/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/fast_bias_field_correction_52/get_anat_restore/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths.nii.gz"
    T1w_LowResMesh_add_spec_file.inputs.structure = "INVALID"
   
    
    ####### T2 image ###
    
    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$T1wFolder"/"$T1wImage".nii.gz
    
    #T2w_add_spec_file_imports = ['import os', 'import subprocess']
    #T2w_add_spec_file = pe.Node(util.Function(input_names=['in_file'],
    #                                          output_names=['out_file'],
    #                                          function=run_add_to_specfile_native,
    #                                         imports=T2w_add_spec_file_imports),
    #                            name='T2w_add_spec_file')
    #T2w_add_spec_file.inputs.in_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/T1w_sub-0050952_ses-1.native.wb.spec"
    
    #wf.add_nodes([T2w_add_spec_file])
    

    #  wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz
    
    #T2w_Atlasspace_add_spec_file_imports = ['import os', 'import subprocess']
    #T2w_Atlasspace_add_spec_file = pe.Node(util.Function(input_names=['in_file'],
    #                                         output_names=['out_file'],
    #                                          function=run_add_to_specfile_atlas,
    #                                         imports=T2w_Atlasspace_add_spec_file_imports),
    #                             name='T2w_Atlasspace_add_spec_file')
    #T2w_Atlasspace_add_spec_file.inputs.in_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/MNINonLinear_sub-0050952_ses-1.native.wb.spec"
    #wf.add_nodes([T2w_Atlasspace_add_spec_file])

    
   
    #wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec INVALID "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz

    #T2w_Atlasspace_HighResMesh_add_spec_file_imports = ['import os', 'import subprocess']
    #T2w_Atlasspace_HighResMesh_add_spec_file = pe.Node(util.Function(input_names=['in_file'],
    #                                           output_names=['out_file'],
    #                                           function=run_add_to_specfile_atlas,
    #                                          imports=T2w_Atlasspace_HighResMesh_add_spec_file_imports),
    #                             name='Atlasspace_HighResMesh_add_spec_file')
    #T2w_Atlasspace_HighResMesh_add_spec_file.inputs.in_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/sub-0050952_ses-1.164k_fs_LR.wb.spec"
    
    #wf.add_nodes(T2w_[Atlasspace_HighResMesh_add_spec_file])
    
    #T2w_Atlasspace_LowResMesh_add_spec_file_imports = ['import os', 'import subprocess']
    #T2w_Atlasspace_LowResMesh_add_spec_file = pe.Node(util.Function(input_names=['in_file'],
    #                                          output_names=['out_file'],
    #                                          function=run_add_to_specfile_atlas,
    #                                           imports=T2w_Atlasspace_LowResMesh_add_spec_file_imports),
    #                             name='T2w_Atlasspace_LowResMesh_add_spec_file')
    #T2w_Atlasspace_LowResMesh_add_spec_file.inputs.in_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/MNINonLinear_sub-0050952_ses-1.32k_fs_LR.wb.spec"


    #wf.add_nodes([T2w_Atlasspace_LowResMesh_add_spec_file])
   
    #wb_command -add-to-spec-file "$T1wFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec INVALID "$T1wFolder"/"$T1wImage".nii.gz
    #T2w_LowResMesh_add_spec_file_imports = ['import os', 'import subprocess']
    #T2w_LowResMesh_add_spec_file = pe.Node(util.Function(input_names=['in_file'],
    #                                              output_names=['out_file'],
    #                                               function=run_add_to_specfile_native,
    #                                               imports=T2w_LowResMesh_add_spec_file_imports),
    #                                name='T2w_LowResMesh_add_spec_file')
    #T2w_LowResMesh_add_spec_file.inputs.in_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec-files/T1w_sub-0050952_ses-1.32k_fs_LR.wb.spec"
    
    
     #applywarp --interp=nn -i "$AtlasSpaceFolder"/wmparc.nii.gz -r "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz
    
    subcortical_apply_warp_1 = pe.Node(interface=fsl.ApplyWarp(),name='subcortical_apply_warp_1')
    subcortical_apply_warp_1.inputs.interp = 'nn'
    subcortical_apply_warp_1.inputs.ref_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/Greyordinates/Atlas_ROIs.2.nii.gz"
    subcortical_apply_warp_1.out_file = "wmparc.2.nii.gz"
    wf.connect(apply_warp_2_wmparc, 'out_file', subcortical_apply_warp_1, 'in_file')
    
    #wb_command -volume-label-import "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz -drop-unused-labels
    subcortical_volume_label_1_imports = ['import os', 'import subprocess']
    subcortical_volume_label_1 = pe.Node(util.Function(input_names=['in_file','file_name'],
                                               output_names=['out_file'],
                                               function=run_volume_label_1,
                                               imports=subcortical_volume_label_1_imports),
                                 name='subcortical_volume_label_1')
                                 
    subcortical_volume_label_1.inputs.file_name = "wmparc.2.nii.gz"
    wf.connect(subcortical_apply_warp_1,'out_file',subcortical_volume_label_1,'in_file')

    
    #applywarp --interp=nn -i "$SurfaceAtlasDIR"/Avgwmparc.nii.gz -r "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$AtlasSpaceFolder"/ROIs/Atlas_wmparc."$GrayordinatesResolution".nii.gz
    subcortical_apply_warp_2 = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_2')
    subcortical_apply_warp_2.inputs.interp = 'nn'
    subcortical_apply_warp_2.inputs.in_file = '/data3/cnl/tgeorge/pfreesurfer-bash-runs/standard_mesh_atlases/Avgwmparc.nii.gz'
    subcortical_apply_warp_2.inputs.ref_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/Greyordinates/Atlas_ROIs.2.nii.gz"
    subcortical_apply_warp_2.out_file = "Atlas_wmparc.2.nii.gz"
    
    
    
  #wb_command -volume-label-import "$AtlasSpaceFolder"/ROIs/Atlas_wmparc."$GrayordinatesResolution".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/ROIs/Atlas_wmparc."$GrayordinatesResolution".nii.gz -drop-unused-labels
  
    subcortical_volume_label_2_imports = ['import os', 'import subprocess']
    subcortical_volume_label_2 = pe.Node(util.Function(input_names=['in_file','file_name'],
                                               output_names=['out_file'],
                                               function=run_volume_label_1,
                                               imports=subcortical_volume_label_2_imports),
                                 name='volume_label_2')
    subcortical_volume_label_2.inputs.file_name = "Atlas_wmparc.2.nii.gz"
    wf.connect(subcortical_apply_warp_2,'out_file',subcortical_volume_label_2,'in_file')
    
  #wb_command -volume-label-import "$AtlasSpaceFolder"/ROIs/wmparc."$GrayordinatesResolution".nii.gz ${SubcorticalGrayLabels} "$AtlasSpaceFolder"/ROIs/ROIs."$GrayordinatesResolution".nii.gz -discard-others
    subcortical_volume_label_3_imports = ['import os', 'import subprocess']
    subcortical_volume_label_3 = pe.Node(util.Function(input_names=['in_file'],
                                               output_names=['out_file'],
                                               function=run_volume_label_2,
                                               imports=subcortical_volume_label_3_imports),
                                 name='subcortical_volume_label_3')
                                 
    wf.connect(subcortical_volume_label_1,'out_file',subcortical_volume_label_3,'in_file')
   
 # applywarp --interp=spline -i "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage".nii.gz -r "$AtlasSpaceFolder"/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage"."$GrayordinatesResolution".nii.gz

    subcortical_apply_warp_3 = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_3')
    subcortical_apply_warp_3.inputs.interp = 'spline'
    subcortical_apply_warp_3.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/fast_bias_field_correction_52/get_anat_restore/sub-0050952_T1w_resample_noise_corrected_corrected_warp_maths_restore_maths.nii.gz'
    subcortical_apply_warp_3.inputs.ref_file = "/data3/cnl/tgeorge/pfreesurfer-bash-runs/Greyordinates/Atlas_ROIs.2.nii.gz"
    subcortical_apply_warp_3.out_file = "2.nii.gz"
 
     # Left White surface

    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii

    l_white_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_white_mris_convert')
    l_white_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.white'
    l_white_mris_convert.inputs.out_file = 'L.white.native.surf.gii'

     #wb_command -set-structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type $Type$Secondary

    l_white_set_structure_imports = ['import os', 'import subprocess']
    l_white_set_structure = pe.Node(util.Function(input_names=['in_file','structure', 'surface_secondary_type'],
                                              output_names=['out_file'],
                                              function=run_set_structure_3,
                                              imports=l_white_set_structure_imports),
                                name='l_white_set_structure')
    l_white_set_structure.inputs.structure = "CORTEX_LEFT"
    l_white_set_structure.inputs.surface_secondary_type = "GRAY_WHITE"
    wf.connect(l_white_mris_convert,'converted',l_white_set_structure,'in_file')

    #wb_command -surface-apply-affine "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$FreeSurferFolder"/mri/c_ras.mat "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    l_white_apply_affine_import = ['import os', 'import subprocess']
    l_white_apply_affine= pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_surface_apply_affine,
                                              imports=l_white_apply_affine_import),
                                name='l_white_apply_affine')
    l_white_apply_affine.inputs.outfile_name = "L.white.native.surf.gii"
    wf.connect(l_white_set_structure,'out_file',l_white_apply_affine,'in_file')


    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    l_white_add_spec_file_T1w_import = ['import os', 'import subprocess']
    l_white_add_spec_file_T1w = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=l_white_add_spec_file_T1w_import),
                                name='l_white_add_spec_file_T1w')

    
    l_white_add_spec_file_T1w.inputs.structure = 'CORTEX_LEFT'
    wf.connect(T1w_add_spec_file, 'out_file', l_white_add_spec_file_T1w, 'spec_file')
    wf.connect(l_white_apply_affine,'out_file',l_white_add_spec_file_T1w,'in_file')



    #wb_command -surface-apply-warpfield "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$InverseAtlasTransform".nii.gz "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii -fnirt "$AtlasTransform".nii.gz
    l_white_apply_warpfield_import = ['import os', 'import subprocess']
    l_white_apply_warpfield = pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=surface_apply_warpfield,
                                              imports=l_white_apply_warpfield_import),
                                name='l_white_apply_warpfield')

    l_white_apply_warpfield.inputs.outfile_name = 'L.white.native.surf.gii'
    wf.connect(l_white_apply_affine,'out_file',l_white_apply_warpfield,'in_file')


    #wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    l_white_add_spec_file_Atlasspace_import = ['import os', 'import subprocess']
    l_white_add_spec_file_Atlasspace = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=l_white_add_spec_file_Atlasspace_import),
                                name='l_white_add_spec_file_Atlasspace')

   
    l_white_add_spec_file_Atlasspace.inputs.structure = 'CORTEX_LEFT'
    wf.connect(l_white_apply_warpfield,'out_file',l_white_add_spec_file_Atlasspace,'in_file')
    wf.connect(Atlasspace_add_spec_file,'out_file',l_white_add_spec_file_Atlasspace,'spec_file')
   

    #  Left Pial surface

    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii


    l_pial_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_pial_mris_convert')
    l_pial_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.pial'
    l_pial_mris_convert.inputs.out_file = 'L.pial.native.surf.gii'

     #wb_command -set-structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type $Type$Secondary

    l_pial_set_structure_imports = ['import os', 'import subprocess']
    l_pial_set_structure = pe.Node(util.Function(input_names=['in_file','structure', 'surface_secondary_type'],
                                              output_names=['out_file'],
                                              function=run_set_structure_3,
                                              imports=l_pial_set_structure_imports),
                                name='l_pial_set_structure')
    l_pial_set_structure.inputs.structure = "CORTEX_LEFT"
    l_pial_set_structure.inputs.surface_secondary_type = "PIAL"
    wf.connect(l_pial_mris_convert,'converted',l_pial_set_structure,'in_file')

    #wb_command -surface-apply-affine "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$FreeSurferFolder"/mri/c_ras.mat "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    l_pial_apply_affine_import = ['import os', 'import subprocess']
    l_pial_apply_affine= pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_surface_apply_affine,
                                              imports=l_pial_apply_affine_import),
                                name='l_pial_apply_affine')
    l_pial_apply_affine.inputs.outfile_name = "L.pial.native.surf.gii"
    wf.connect(l_pial_set_structure,'out_file',l_pial_apply_affine,'in_file')


    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    l_pial_add_spec_file_T1w_import = ['import os', 'import subprocess']
    l_pial_add_spec_file_T1w = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=l_pial_add_spec_file_T1w_import),
                                name='l_pial_add_spec_file_T1w')

    
    l_pial_add_spec_file_T1w.inputs.structure = 'CORTEX_LEFT'
    wf.connect(l_pial_apply_affine,'out_file',l_pial_add_spec_file_T1w,'in_file')
    wf.connect(l_white_add_spec_file_T1w,'out_file',l_pial_add_spec_file_T1w,'spec_file')


    #wb_command -surface-apply-warpfield "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$InverseAtlasTransform".nii.gz "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii -fnirt "$AtlasTransform".nii.gz
    l_pial_apply_warpfield_import = ['import os', 'import subprocess']
    l_pial_apply_warpfield = pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=surface_apply_warpfield,
                                              imports=l_pial_apply_warpfield_import),
                                name='l_pial_apply_warpfield')

    l_pial_apply_warpfield.inputs.outfile_name = 'L.white.native.surf.gii'
    wf.connect(l_pial_apply_affine,'out_file',l_pial_apply_warpfield,'in_file')


    #wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    l_pial_add_spec_file_Atlasspace_import = ['import os', 'import subprocess']
    l_pial_add_spec_file_Atlasspace = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=l_pial_add_spec_file_Atlasspace_import),
                                name='l_pial_add_spec_file_Atlasspace')

    
    l_pial_add_spec_file_Atlasspace.inputs.structure = 'CORTEX_LEFT'
    wf.connect(l_pial_apply_warpfield,'out_file',l_pial_add_spec_file_Atlasspace,'in_file')
    wf.connect(l_white_add_spec_file_Atlasspace,'out_file',l_pial_add_spec_file_Atlasspace,'spec_file')



    # Right White surface

    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii

    r_white_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_white_mris_convert')
    r_white_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.white'
    r_white_mris_convert.inputs.out_file = 'R.white.native.surf.gii'


    #wb_command -set-structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type $Type$Secondary

    r_white_set_structure_imports = ['import os', 'import subprocess']
    r_white_set_structure = pe.Node(util.Function(input_names=['in_file','structure', 'surface_secondary_type'],
                                            output_names=['out_file'],
                                            function=run_set_structure_3,
                                            imports=r_white_set_structure_imports),
                              name='r_white_set_structure')
    r_white_set_structure.inputs.structure = "CORTEX_RIGHT"
    r_white_set_structure.inputs.surface_secondary_type = "GRAY_WHITE"
    wf.connect(r_white_mris_convert,'converted',r_white_set_structure,'in_file')

    #wb_command -surface-apply-affine "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$FreeSurferFolder"/mri/c_ras.mat "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    r_white_apply_affine_import = ['import os', 'import subprocess']
    r_white_apply_affine= pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                            output_names=['out_file'],
                                            function=run_surface_apply_affine,
                                            imports=r_white_apply_affine_import),
                              name='r_white_apply_affine')
    r_white_apply_affine.inputs.outfile_name = "R.white.native.surf.gii"
    wf.connect(r_white_set_structure,'out_file',r_white_apply_affine,'in_file')


    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    r_white_add_spec_file_T1w_import = ['import os', 'import subprocess']
    r_white_add_spec_file_T1w = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                            output_names=['out_file'],
                                            function=add_spec_file,
                                            imports=r_white_add_spec_file_T1w_import),
                              name='r_white_add_spec_file_T1w')

  
    r_white_add_spec_file_T1w.inputs.structure = 'CORTEX_RIGHT'
    wf.connect(r_white_apply_affine,'out_file',r_white_add_spec_file_T1w,'in_file')
    wf.connect(l_pial_add_spec_file_T1w,'out_file',r_white_add_spec_file_T1w,'spec_file')


    #wb_command -surface-apply-warpfield "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$InverseAtlasTransform".nii.gz "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii -fnirt "$AtlasTransform".nii.gz
    r_white_apply_warpfield_import = ['import os', 'import subprocess']
    r_white_apply_warpfield = pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                            output_names=['out_file'],
                                            function=surface_apply_warpfield,
                                            imports=r_white_apply_warpfield_import),
                              name='r_white_apply_warpfield')

    r_white_apply_warpfield.inputs.outfile_name = 'R.white.native.surf.gii'
    wf.connect(r_white_apply_affine,'out_file',r_white_apply_warpfield,'in_file')


    #wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    r_white_add_spec_file_Atlasspace_import = ['import os', 'import subprocess']
    r_white_add_spec_file_Atlasspace = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                            output_names=['out_file'],
                                            function=add_spec_file,
                                            imports=r_white_add_spec_file_Atlasspace_import),
                              name='r_white_add_spec_file_Atlasspace')


    r_white_add_spec_file_Atlasspace.inputs.structure = 'CORTEX_RIGHT'
    wf.connect(r_white_apply_warpfield,'out_file',r_white_add_spec_file_Atlasspace,'in_file')
    wf.connect(l_pial_add_spec_file_Atlasspace,'out_file',r_white_add_spec_file_Atlasspace,'spec_file')

  #  Right Pial surface

    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii


    r_pial_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_pial_mris_convert')
    r_pial_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.pial'
    r_pial_mris_convert.inputs.out_file = 'R.pial.native.surf.gii'


    #wb_command -set-structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type $Type$Secondary

    r_pial_set_structure_imports = ['import os', 'import subprocess']
    r_pial_set_structure = pe.Node(util.Function(input_names=['in_file','structure', 'surface_secondary_type'],
                                            output_names=['out_file'],
                                            function=run_set_structure_3,
                                            imports=r_pial_set_structure_imports),
                              name='r_pial_set_structure')
    r_pial_set_structure.inputs.structure = "CORTEX_RIGHT"
    r_pial_set_structure.inputs.surface_secondary_type = "PIAL"
    wf.connect(r_pial_mris_convert,'converted',r_pial_set_structure,'in_file')

    #wb_command -surface-apply-affine "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$FreeSurferFolder"/mri/c_ras.mat "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    r_pial_apply_affine_import = ['import os', 'import subprocess']
    r_pial_apply_affine = pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                            output_names=['out_file'],
                                            function=run_surface_apply_affine,
                                            imports=r_pial_apply_affine_import),
                              name='r_pial_apply_affine')
    r_pial_apply_affine.inputs.outfile_name = "R.pial.native.surf.gii"
    wf.connect(r_pial_set_structure,'out_file',r_pial_apply_affine,'in_file')


    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    r_pial_add_spec_file_T1w_import = ['import os', 'import subprocess']
    r_pial_add_spec_file_T1w = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                            output_names=['out_file'],
                                            function=add_spec_file,
                                            imports=r_pial_add_spec_file_T1w_import),
                              name='r_pial_add_spec_file_T1w')

   
    r_pial_add_spec_file_T1w.inputs.structure = 'CORTEX_RIGHT'
    wf.connect(r_pial_apply_affine,'out_file',r_pial_add_spec_file_T1w,'in_file')
    wf.connect(r_white_add_spec_file_T1w,'out_file',r_pial_add_spec_file_T1w,'spec_file')


    #wb_command -surface-apply-warpfield "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$InverseAtlasTransform".nii.gz "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii -fnirt "$AtlasTransform".nii.gz
    r_pial_apply_warpfield_import = ['import os', 'import subprocess']
    r_pial_apply_warpfield = pe.Node(util.Function(input_names=['in_file','outfile_name'],
                                            output_names=['out_file'],
                                            function=surface_apply_warpfield,
                                            imports=r_pial_apply_warpfield_import),
                              name='r_pial_apply_warpfield')

    r_pial_apply_warpfield.inputs.outfile_name = 'R.pial.native.surf.gii'
    wf.connect(r_pial_apply_affine,'out_file',r_pial_apply_warpfield,'in_file')


    #wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    r_pial_add_spec_file_Atlasspace_import = ['import os', 'import subprocess']
    r_pial_add_spec_file_Atlasspace = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                            output_names=['out_file'],
                                            function=add_spec_file,
                                            imports=r_pial_add_spec_file_Atlasspace_import),
                              name='r_pial_add_spec_file_Atlasspace')

    
    r_pial_add_spec_file_Atlasspace.inputs.structure = 'CORTEX_RIGHT'
    wf.connect(r_pial_apply_warpfield,'out_file',r_pial_add_spec_file_Atlasspace,'in_file')
    wf.connect(r_white_add_spec_file_Atlasspace,'out_file',r_pial_add_spec_file_Atlasspace,'spec_file')
   
  
   ####### Block 7: #####
   # T1W left
   
    #wb_command -surface-average "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii


    T1w_l_surface_avg_imports = ['import os', 'import subprocess']
    T1w_l_surface_avg = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                              output_names=['out_file'],
                                              function=run_surface_avg,
                                              imports=T1w_l_surface_avg_imports),
                                name='T1w_l_surface_avg')

    
    wf.connect(l_white_apply_affine,'out_file',T1w_l_surface_avg, 'in_file_1')
    wf.connect(l_pial_apply_affine,'out_file', T1w_l_surface_avg, 'in_file_2')

    #wb_command -set-structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii ${Structure} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS

    T1w_l_set_structure_imports = ['import os', 'import subprocess']
    T1w_l_set_structure = pe.Node(util.Function(input_names=['in_file','structure', 'surface_secondary_type'],
                                              output_names=['out_file'],
                                              function=run_set_structure_3,
                                              imports=T1w_l_set_structure_imports),
                                name='T1w_l_set_structure')


    T1w_l_set_structure.inputs.structure = "CORTEX_LEFT"
    T1w_l_set_structure.inputs.surface_secondary_type = "MIDTHICKNESS"
    wf.connect(T1w_l_surface_avg,'out_file',T1w_l_set_structure,'in_file')


    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii


    T1w_l_spec_file_imports_1 = ['import os', 'import subprocess']
    T1w_l_spec_file_1 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=T1w_l_spec_file_imports_1),
                               name='T1w_l_spec_file_1')

    T1w_l_spec_file_1.inputs.structure = "CORTEX_LEFT"
    T1w_l_spec_file_1.inputs.spec_file = '/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/spec_files_blk7/T1w_sub-0050952_ses-1.native.wb.spec'
    wf.connect(T1w_l_set_structure,'out_file',T1w_l_spec_file_1,'in_file')
    wf.connect(r_white_add_spec_file_T1w,'out_file',T1w_l_spec_file_1,'spec_file')


    #get number of vertices from native file
    #NativeVerts=$(wb_command -file-information "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii | grep 'Number of Vertices:' | cut -f2 -d: | tr -d '[:space:]')

    T1w_l_nativeverts_imports = ['import os', 'import subprocess']
    T1w_l_nativeverts = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=get_native_verts,
                                              imports=T1w_l_nativeverts_imports),
                                name='T1w_l_nativeverts')

   
    wf.connect(T1w_l_set_structure,'out_file',T1w_l_nativeverts,'in_file')
    
    #HCP fsaverage_LR32k used -iterations-scale 0.75. Compute new param value for native mesh density
    #NativeInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $NativeVerts / 32492" | bc -l)
    
    T1w_l_nativeinflationscale_imports = ['import os', 'import subprocess']
    T1w_l_nativeinflationscale = pe.Node(util.Function(input_names=['InflateExtraScale','in_val'],
                                              output_names=['out_file'],
                                              function=get_nativeinflationscale,
                                              imports=T1w_l_nativeinflationscale_imports),
                                name='T1w_l_nativeinflationscale')
    
    
    T1w_l_nativeinflationscale.inputs.InflateExtraScale = 1
    wf.connect(T1w_l_nativeverts,'out_file',T1w_l_nativeinflationscale,'in_val')
    
    

    #wb_command -surface-generate-inflated "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii       #"$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii -iterations-scale $NativeInflationScale



    T1w_l_sur_inflated_imports = ['import os', 'import subprocess']
    T1w_l_sur_inflated = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                            output_names=['out_file_1','out_file_2'],
                                            function=run_generate_sur_inflated,
                                            imports=T1w_l_sur_inflated_imports),
                               name='T1w_l_sur_inflated')


    wf.connect(T1w_l_set_structure,'out_file',T1w_l_sur_inflated,'in_file_1')
    wf.connect(T1w_l_nativeinflationscale,'out_file',T1w_l_sur_inflated,'in_file_2')

    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii

    T1w_l_spec_file_imports_2 = ['import os', 'import subprocess']
    T1w_l_spec_file_2 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=T1w_l_spec_file_imports_2),
                                name='T1w_l_spec_file_2')

    
    T1w_l_spec_file_2.inputs.structure = "CORTEX_LEFT"
    wf.connect(T1w_l_spec_file_1,'out_file',T1w_l_spec_file_2,'spec_file')
    wf.connect(T1w_l_sur_inflated,'out_file_1',T1w_l_spec_file_2,'in_file')



   #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii


    T1w_l_spec_file_imports_3 = ['import os', 'import subprocess']
    T1w_l_spec_file_3 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=T1w_l_spec_file_imports_3),
                                    name='T1w_l_spec_file_3')

    wf.connect(T1w_l_spec_file_2,'out_file',T1w_l_spec_file_3,'spec_file')
    T1w_l_spec_file_3.inputs.structure = "CORTEX_LEFT"

    wf.connect(T1w_l_sur_inflated,'out_file_2',T1w_l_spec_file_3,'in_file')


   

    ################################################ Atlasspace l hemi ############################################



    #wb_command -surface-average "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii


    Atlasspace_l_surface_avg_imports = ['import os', 'import subprocess']
    Atlasspace_l_surface_avg = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                               output_names=['out_file'],
                                               function=run_surface_avg,
                                               imports=Atlasspace_l_surface_avg_imports),
                              name='Atlasspace_l_surface_avg')

    
    wf.connect(l_white_apply_warpfield, 'out_file',Atlasspace_l_surface_avg, 'in_file_1')
    wf.connect(l_pial_apply_warpfield, 'out_file',Atlasspace_l_surface_avg, 'in_file_2')

    #wb_command -set-structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii ${Structure} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS

    Atlasspace_l_set_structure_imports = ['import os', 'import subprocess']
    Atlasspace_l_set_structure = pe.Node(util.Function(input_names=['in_file','structure','surface_secondary_type'],
                                              output_names=['out_file'],
                                              function=run_set_structure_3,
                                              imports=Atlasspace_l_set_structure_imports),
                                name='Atlasspace_l_set_structure')


    Atlasspace_l_set_structure.inputs.structure = "CORTEX_LEFT"
    Atlasspace_l_set_structure.inputs.surface_secondary_type = "MIDTHICKNESS"
    wf.connect(Atlasspace_l_surface_avg,'out_file',Atlasspace_l_set_structure,'in_file')


    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii


    Atlasspace_l_spec_file_1_imports = ['import os', 'import subprocess']
    Atlasspace_l_spec_file_1 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=Atlasspace_l_spec_file_1_imports),
                               name='Atlasspace_l_spec_file_1')

    Atlasspace_l_spec_file_1.inputs.structure = "CORTEX_LEFT"
    wf.connect(Atlasspace_l_set_structure,'out_file',Atlasspace_l_spec_file_1,'in_file')
    wf.connect(r_pial_add_spec_file_Atlasspace, 'out_file',Atlasspace_l_spec_file_1, 'spec_file')


    #get number of vertices from native file
    #NativeVerts=$(wb_command -file-information "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii | grep 'Number of Vertices:' | cut -f2 -d: | tr -d '[:space:]')

    Atlasspace_l_nativeverts_imports = ['import os', 'import subprocess']
    Atlasspace_l_nativeverts = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=get_native_verts,
                                             imports=Atlasspace_l_nativeverts_imports),
                                name='Atlasspace_l_nativeverts')
                                
    wf.connect(Atlasspace_l_set_structure,'out_file', Atlasspace_l_nativeverts, 'in_file')



    #HCP fsaverage_LR32k used -iterations-scale 0.75. Compute new param value for native mesh density
    #NativeInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $NativeVerts / 32492" | bc -l)
    
    Atlasspace_l_nativeinflationscale_imports = ['import os', 'import subprocess']
    Atlasspace_l_nativeinflationscale = pe.Node(util.Function(input_names=['InflateExtraScale','in_val'],
                                              output_names=['out_file'],
                                              function=get_nativeinflationscale,
                                              imports=Atlasspace_l_nativeinflationscale_imports),
                                name='Atlasspace_l_nativeinflationscale')
    
    
    Atlasspace_l_nativeinflationscale.inputs.InflateExtraScale = 1
    wf.connect(Atlasspace_l_nativeverts,'out_file',Atlasspace_l_nativeinflationscale,'in_val')
    
    
    

	 #wb_command -surface-generate-inflated "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii       #"$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii -iterations-scale $NativeInflationScale



    Atlasspace_l_sur_inflated_imports = ['import os', 'import subprocess']
    Atlasspace_l_sur_inflated = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                            output_names=['out_file_1','out_file_2'],
                                            function=run_generate_sur_inflated,
                                            imports=Atlasspace_l_sur_inflated_imports),
                               name='Atlasspace_l_sur_inflated')


    wf.connect(Atlasspace_l_set_structure,'out_file',Atlasspace_l_sur_inflated,'in_file_1')
    wf.connect(Atlasspace_l_nativeinflationscale,'out_file',Atlasspace_l_sur_inflated,'in_file_2')

    
    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii
    Atlasspace_l_spec_file_imports_2 = ['import os', 'import subprocess']
    Atlasspace_l_spec_file_2 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=Atlasspace_l_spec_file_imports_2),
                                name='Atlasspace_l_spec_file_2')

    Atlasspace_l_spec_file_2.inputs.structure = "CORTEX_LEFT"

    wf.connect(Atlasspace_l_sur_inflated,'out_file_1',Atlasspace_l_spec_file_2,'in_file')
    wf.connect(Atlasspace_l_spec_file_1,'out_file',Atlasspace_l_spec_file_2,'spec_file')

    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii

    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii
    Atlasspace_l_spec_file_imports_3 = ['import os', 'import subprocess']
    Atlasspace_l_spec_file_3 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=Atlasspace_l_spec_file_imports_3),
                                name='Atlasspace_l_spec_file_3')

    Atlasspace_l_spec_file_3.inputs.structure = "CORTEX_LEFT"

    wf.connect(Atlasspace_l_sur_inflated,'out_file_2',Atlasspace_l_spec_file_3,'in_file')
    wf.connect(Atlasspace_l_spec_file_2,'out_file',Atlasspace_l_spec_file_3,'spec_file')



############################################################## T1w right hemi ##############################

   #wb_command -surface-average "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii


    T1w_r_surface_avg_imports = ['import os', 'import subprocess']
    T1w_r_surface_avg = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                              output_names=['out_file'],
                                              function=run_surface_avg,
                                              imports=T1w_r_surface_avg_imports),
                                name='T1w_r_surface_avg')

   
    wf.connect(r_white_apply_affine, 'out_file',T1w_r_surface_avg, 'in_file_1')
    wf.connect(r_pial_apply_affine, 'out_file',T1w_r_surface_avg, 'in_file_2')

    #wb_command -set-structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii ${Structure} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS

    T1w_r_set_structure_imports = ['import os', 'import subprocess']
    T1w_r_set_structure = pe.Node(util.Function(input_names=['in_file','structure','surface_secondary_type'],
                                              output_names=['out_file'],
                                              function=run_set_structure_3,
                                              imports=T1w_r_set_structure_imports),
                                name='T1w_r_set_structure')


    T1w_r_set_structure.inputs.structure = "CORTEX_RIGHT"
    T1w_r_set_structure.inputs.surface_secondary_type = "MIDTHICKNESS"
    wf.connect(T1w_r_surface_avg,'out_file',T1w_r_set_structure,'in_file')


    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii


    T1w_r_spec_file_imports_1 = ['import os', 'import subprocess']
    T1w_r_spec_file_1 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=T1w_r_spec_file_imports_1),
                               name='T1w_r_spec_file_1')

    T1w_r_spec_file_1.inputs.structure = "CORTEX_RIGHT"
    wf.connect(T1w_l_spec_file_3,'out_file',T1w_r_spec_file_1,'spec_file')
    wf.connect(T1w_r_set_structure,'out_file',T1w_r_spec_file_1,'in_file')



    #get number of vertices from native file
    #NativeVerts=$(wb_command -file-information "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii | grep 'Number of Vertices:' | cut -f2 -d: | tr -d '[:space:]')

    T1w_r_nativeverts_imports = ['import os', 'import subprocess']
    T1w_r_nativeverts = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=get_native_verts,
                                              imports=T1w_r_nativeverts_imports),
                                name='T1w_r_nativeverts')

   
    wf.connect(T1w_r_set_structure,'out_file',T1w_r_nativeverts,'in_file')
    
    #HCP fsaverage_LR32k used -iterations-scale 0.75. Compute new param value for native mesh density
    #NativeInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $NativeVerts / 32492" | bc -l)
    
    T1w_r_nativeinflationscale_imports = ['import os', 'import subprocess']
    T1w_r_nativeinflationscale = pe.Node(util.Function(input_names=['InflateExtraScale','in_val'],
                                              output_names=['out_file'],
                                              function=get_nativeinflationscale,
                                              imports=T1w_r_nativeinflationscale_imports),
                                name='T1w_r_nativeinflationscale')
    
    
    T1w_r_nativeinflationscale.inputs.InflateExtraScale = 1
    wf.connect(T1w_r_nativeverts,'out_file',T1w_r_nativeinflationscale,'in_val')
    
    

    #wb_command -surface-generate-inflated "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii       #"$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii -iterations-scale $NativeInflationScale



    T1w_r_sur_inflated_imports = ['import os', 'import subprocess']
    T1w_r_sur_inflated = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                            output_names=['out_file_1','out_file_2'],
                                            function=run_generate_sur_inflated,
                                            imports=T1w_r_sur_inflated_imports),
                               name='T1w_r_sur_inflated')


    wf.connect(T1w_r_set_structure,'out_file',T1w_r_sur_inflated,'in_file_1')
    wf.connect(T1w_r_nativeinflationscale,'out_file',T1w_r_sur_inflated,'in_file_2')

    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii

    T1w_r_spec_file_imports_2 = ['import os', 'import subprocess']
    T1w_r_spec_file_2 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=T1w_r_spec_file_imports_2),
                                name='T1w_r_spec_file_2')

    wf.connect(T1w_r_spec_file_1,'out_file', T1w_r_spec_file_2,'spec_file')
    T1w_r_spec_file_2.inputs.structure = "CORTEX_LEFT"

    wf.connect(T1w_r_sur_inflated,'out_file_1',T1w_r_spec_file_2,'in_file')



   #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii


    T1w_r_spec_file_imports_3 = ['import os', 'import subprocess']
    T1w_r_spec_file_3 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=T1w_r_spec_file_imports_3),
                                    name='T1w_r_spec_file_3')

    wf.connect(T1w_r_spec_file_2,'out_file', T1w_r_spec_file_3,'spec_file')
    T1w_r_spec_file_3.inputs.structure = "CORTEX_LEFT"

    wf.connect(T1w_r_sur_inflated,'out_file_2',T1w_r_spec_file_3,'in_file')
    
    
    ######################################### Atlasspace R hemi #######################################




    #wb_command -surface-average "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii -surf "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii


    Atlasspace_r_surface_avg_imports = ['import os', 'import subprocess']
    Atlasspace_r_surface_avg = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                               output_names=['out_file'],
                                               function=run_surface_avg,
                                               imports=Atlasspace_r_surface_avg_imports),
                              name='Atlasspace_r_surface_avg')


    wf.connect(r_pial_apply_warpfield, 'out_file',Atlasspace_r_surface_avg, 'in_file_1')
    wf.connect(r_white_apply_warpfield,'out_file',Atlasspace_r_surface_avg, 'in_file_2')

    #wb_command -set-structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii ${Structure} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS

    Atlasspace_r_set_structure_imports = ['import os', 'import subprocess']
    Atlasspace_r_set_structure = pe.Node(util.Function(input_names=['in_file','structure','surface_secondary_type'],
                                              output_names=['out_file'],
                                              function=run_set_structure_3,
                                              imports=Atlasspace_r_set_structure_imports),
                                name='Atlasspace_r_set_structure')

    Atlasspace_r_set_structure.inputs.surface_secondary_type = "MIDTHICKNESS"
    Atlasspace_r_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(Atlasspace_r_surface_avg,'out_file',Atlasspace_r_set_structure,'in_file')


    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii


    Atlasspace_r_spec_file_1_imports = ['import os', 'import subprocess']
    Atlasspace_r_spec_file_1 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=Atlasspace_r_spec_file_1_imports),
                               name='Atlasspace_r_spec_file_1')

    Atlasspace_r_spec_file_1.inputs.structure = "CORTEX_RIGHT"
    wf.connect(Atlasspace_l_spec_file_3,'out_file', Atlasspace_r_spec_file_1, 'spec_file')
    wf.connect(Atlasspace_r_set_structure,'out_file',Atlasspace_r_spec_file_1,'in_file')



    #get number of vertices from native file
    #NativeVerts=$(wb_command -file-information "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii | grep 'Number of Vertices:' | cut -f2 -d: | tr -d '[:space:]')

    Atlasspace_r_nativeverts_imports = ['import os', 'import subprocess']
    Atlasspace_r_nativeverts = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=get_native_verts,
                                             imports=Atlasspace_r_nativeverts_imports),
                                name='Atlasspace_r_nativeverts')
                                
    wf.connect(Atlasspace_r_set_structure,'out_file', Atlasspace_r_nativeverts, 'in_file')



    #HCP fsaverage_LR32k used -iterations-scale 0.75. Compute new param value for native mesh density
    #NativeInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $NativeVerts / 32492" | bc -l)
    
    Atlasspace_r_nativeinflationscale_imports = ['import os', 'import subprocess']
    Atlasspace_r_nativeinflationscale = pe.Node(util.Function(input_names=['InflateExtraScale','in_val'],
                                              output_names=['out_file'],
                                              function=get_nativeinflationscale,
                                              imports=Atlasspace_r_nativeinflationscale_imports),
                                name='Atlasspace_r_nativeinflationscale')
    
    
    Atlasspace_r_nativeinflationscale.inputs.InflateExtraScale = 1
    wf.connect(Atlasspace_r_nativeverts,'out_file',Atlasspace_r_nativeinflationscale,'in_val')
    
    
    

	 #wb_command -surface-generate-inflated "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii       #"$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii -iterations-scale $NativeInflationScale



    Atlasspace_r_sur_inflated_imports = ['import os', 'import subprocess']
    Atlasspace_r_sur_inflated = pe.Node(util.Function(input_names=['in_file_1', 'in_file_2'],
                                            output_names=['out_file_1','out_file_2'],
                                            function=run_generate_sur_inflated,
                                            imports=Atlasspace_r_sur_inflated_imports),
                               name='Atlasspace_r_sur_inflated')


    wf.connect(Atlasspace_r_set_structure,'out_file',Atlasspace_r_sur_inflated,'in_file_1')
    wf.connect(Atlasspace_r_nativeinflationscale,'out_file',Atlasspace_r_sur_inflated,'in_file_2')

    
    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii
    Atlasspace_r_spec_file_imports_2 = ['import os', 'import subprocess']
    Atlasspace_r_spec_file_2 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=Atlasspace_r_spec_file_imports_2),
                                name='Atlasspace_r_spec_file_2')

    Atlasspace_r_spec_file_2.inputs.structure = "CORTEX_LEFT"

    wf.connect(Atlasspace_r_sur_inflated,'out_file_1',Atlasspace_r_spec_file_2,'in_file')
    wf.connect(Atlasspace_r_spec_file_1,'out_file',Atlasspace_r_spec_file_2,'spec_file')

    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii

    #wb_command -add-to-spec-file "$Folder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$Folder"/"$NativeFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii
    Atlasspace_r_spec_file_imports_3 = ['import os', 'import subprocess']
    Atlasspace_r_spec_file_3 = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                             output_names=['out_file'],
                                             function=add_spec_file,
                                             imports=Atlasspace_r_spec_file_imports_3),
                                name='Atlasspace_r_spec_file_3')

    Atlasspace_r_spec_file_3.inputs.structure = "CORTEX_LEFT"

    wf.connect(Atlasspace_r_sur_inflated,'out_file_2',Atlasspace_r_spec_file_3,'in_file')
    wf.connect(Atlasspace_r_spec_file_2,'out_file',Atlasspace_r_spec_file_3,'spec_file')
    

     ### Left hemisphere ####
    # For sphere.reg
    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
     
    l_spherereg_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_spherereg_mris_convert')
    l_spherereg_mris_convert.inputs.in_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.sphere.reg'
    l_spherereg_mris_convert.inputs.out_file = 'L.spherereg.native.surf.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type SPHERICAL
    l_spherereg_set_structure_imports = ['import os', 'import subprocess']
    l_spherereg_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_2,
                                              imports=l_spherereg_set_structure_imports),
                                name='l_spherereg_set_structure')
                                
    l_spherereg_set_structure.inputs.structure = "CORTEX_LEFT"
    wf.connect(l_spherereg_mris_convert,'converted',l_spherereg_set_structure,'in_file')
    
    # For sphere
     
    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
     
    l_sphere_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_sphere_mris_convert')
    l_sphere_mris_convert.inputs.in_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.sphere'
    l_sphere_mris_convert.inputs.out_file = 'L.sphere.native.surf.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type SPHERICAL
    l_sphere_set_structure_imports = ['import os', 'import subprocess']
    l_sphere_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_2,
                                              imports=l_sphere_set_structure_imports),
                                name='l_sphere_set_structure')
                                
    l_sphere_set_structure.inputs.structure = "CORTEX_LEFT"
    wf.connect(l_sphere_mris_convert,'converted',l_sphere_set_structure,'in_file')
    
    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    l_sphere_add_spec_file_import = ['import os', 'import subprocess']
    l_sphere_add_spec_file = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=l_sphere_add_spec_file_import),
                                name='l_sphere_add_spec_file')

    
    l_sphere_add_spec_file.inputs.structure = 'CORTEX_LEFT'
    wf.connect(l_sphere_set_structure,'out_file',l_sphere_add_spec_file,'in_file')
    wf.connect(Atlasspace_r_spec_file_3, 'out_file',l_sphere_add_spec_file, 'spec_file')
    
    # For sulc 
    #mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    l_sulc_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_sulc_mris_convert')
    l_sulc_mris_convert.inputs.scalarcurv_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.sulc'
    l_sulc_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.white'
    l_sulc_mris_convert.inputs.out_file = 'L.sulc.native.shape.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}
    
    l_sulc_set_structure_imports = ['import os', 'import subprocess']
    l_sulc_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=l_sulc_set_structure_imports),
                                name='l_sulc_set_structure')
                                
    l_sulc_set_structure.inputs.structure = "CORTEX_LEFT"
    wf.connect(l_sulc_mris_convert,'converted',l_sulc_set_structure,'in_file')   
    
    
    #wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    l_sulc_metric_maths_imports = ['import os', 'import subprocess']
    l_sulc_metric_maths = pe.Node(util.Function(input_names=['in_file', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_3,
                                              imports=l_sulc_metric_maths_imports),
                                name='l_sulc_metric_maths')
                                
    l_sulc_metric_maths.inputs.outfile_name = 'sulc.native.shape.gii'
    wf.connect(l_sulc_set_structure,'out_file',l_sulc_metric_maths,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
    l_sulc_set_maps_imports = ['import os', 'import subprocess']
    l_sulc_set_maps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_maps,
                                              imports=l_sulc_set_maps_imports),
                                name='l_sulc_set_maps')
                                
    l_sulc_set_maps.inputs.map_name = 'subject_L_Sulc'
    wf.connect(l_sulc_metric_maths,'out_file',l_sulc_set_maps,'in_file')   
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
    l_sulc_metric_palette_imports = ['import os', 'import subprocess']
    l_sulc_metric_palette = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=l_sulc_metric_palette_imports),
                                name='l_sulc_metric_palette')
                                
    
    wf.connect(l_sulc_set_maps,'out_file',l_sulc_metric_palette,'in_file') 
    
    
    # For thickness
    #mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    l_thickness_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_thickness_mris_convert')
    l_thickness_mris_convert.inputs.scalarcurv_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.thickness'
    l_thickness_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.white'
    l_thickness_mris_convert.inputs.out_file = 'L.thickness.native.shape.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}
    l_thickness_set_structure_imports = ['import os', 'import subprocess']
    l_thickness_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=l_thickness_set_structure_imports),
                                name='l_thickness_set_structure')
                                
    l_thickness_set_structure.inputs.structure = "CORTEX_LEFT"
    wf.connect(l_thickness_mris_convert,'converted',l_thickness_set_structure,'in_file')   
    
    #wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    l_thickness_metric_maths_imports = ['import os', 'import subprocess']
    l_thickness_metric_maths = pe.Node(util.Function(input_names=['in_file', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_3,
                                              imports=l_thickness_metric_maths_imports),
                                name='l_thickness_metric_maths')
                                
    l_thickness_metric_maths.inputs.outfile_name = 'thickness.native.shape.gii'
    wf.connect(l_thickness_set_structure,'out_file',l_thickness_metric_maths,'in_file')   
    
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
    l_thickness_set_maps_imports = ['import os', 'import subprocess']
    l_thickness_set_maps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_maps,
                                              imports=l_thickness_set_maps_imports),
                                name='l_thickness_set_maps')
                                
    l_thickness_set_maps.inputs.map_name = 'subject_L_Thickness'
    wf.connect(l_thickness_metric_maths,'out_file',l_thickness_set_maps,'in_file')   
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
    l_thickness_metric_palette_imports_2 = ['import os', 'import subprocess']
    l_thickness_metric_palette_2 = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=l_thickness_metric_palette_imports_2),
                                name='l_thickness_metric_palette_2')
                                
    
    wf.connect(l_thickness_set_maps,'out_file',l_thickness_metric_palette_2,'in_file') 
    
    
    # For curvature
    #mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    l_curv_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_curv_mris_convert')
    l_curv_mris_convert.inputs.scalarcurv_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.curv'
    l_curv_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.white'
    l_curv_mris_convert.inputs.out_file = 'L.curvature.native.shape.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}

    l_curv_set_structure_imports = ['import os', 'import subprocess']
    l_curv_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=l_curv_set_structure_imports),
                                name='l_curv_set_structure')
                                
    l_curv_set_structure.inputs.structure = "CORTEX_LEFT"
    wf.connect(l_curv_mris_convert,'converted',l_curv_set_structure,'in_file')   
    
    #wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    l_curv_metric_maths_imports = ['import os', 'import subprocess']
    l_curv_metric_maths = pe.Node(util.Function(input_names=['in_file', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_3,
                                              imports=l_curv_metric_maths_imports),
                                name='l_curv_metric_maths')
                                
    l_curv_metric_maths.inputs.outfile_name = 'curv.native.shape.gii'
    wf.connect(l_curv_set_structure,'out_file',l_curv_metric_maths,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
    l_curv_set_maps_imports = ['import os', 'import subprocess']
    l_curv_set_maps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_maps,
                                              imports=l_curv_set_maps_imports),
                                name='l_curv_set_maps')
                                
    l_curv_set_maps.inputs.map_name = 'subject_L_Curvature'
    wf.connect(l_curv_metric_maths,'out_file',l_curv_set_maps,'in_file')   
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
    
    l_curv_metric_palette_imports = ['import os', 'import subprocess']
    l_curv_metric_palette = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=l_curv_metric_palette_imports),
                                name='l_curv_metric_palette')
                                
    
    wf.connect(l_curv_set_maps,'out_file',l_curv_metric_palette,'in_file') 
    

    
    
    
    ### Right hemisphere ####
    #For sphere.reg
    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
     
    r_spherereg_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_spherereg_mris_convert')
    r_spherereg_mris_convert.inputs.in_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.sphere.reg'
    r_spherereg_mris_convert.inputs.out_file = 'R.spherereg.native.surf.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type SPHERICAL
    r_spherereg_set_structure_imports = ['import os', 'import subprocess']
    r_spherereg_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_2,
                                              imports=r_spherereg_set_structure_imports),
                                name='r_spherereg_set_structure')
                                
    r_spherereg_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(r_spherereg_mris_convert,'converted',r_spherereg_set_structure,'in_file')
    
    # For sphere
     
    #mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
     
    r_sphere_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_sphere_mris_convert')
    r_sphere_mris_convert.inputs.in_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.sphere'
    r_sphere_mris_convert.inputs.out_file = 'R.sphere.native.surf.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type SPHERICAL
    r_sphere_set_structure_imports = ['import os', 'import subprocess']
    r_sphere_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_2,
                                              imports=r_sphere_set_structure_imports),
                                name='r_sphere_set_structure')
                                
    r_sphere_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(r_sphere_mris_convert,'converted',r_sphere_set_structure,'in_file')
    
    #wb_command -add-to-spec-file "$T1wFolder"/"$NativeFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    r_sphere_add_spec_file_import = ['import os', 'import subprocess']
    r_sphere_add_spec_file = pe.Node(util.Function(input_names=['spec_file','in_file','structure'],
                                              output_names=['out_file'],
                                              function=add_spec_file,
                                              imports=r_sphere_add_spec_file_import),
                                name='r_sphere_add_spec_file')

    
    r_sphere_add_spec_file.inputs.structure = 'CORTEX_RIGHT'
    wf.connect(l_sphere_add_spec_file,'out_file',r_sphere_add_spec_file,'spec_file')
    wf.connect(r_sphere_set_structure,'out_file',r_sphere_add_spec_file,'in_file')

    
    # For sulc 
    #mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    r_sulc_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_sulc_mris_convert')
    r_sulc_mris_convert.inputs.scalarcurv_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.sulc'
    r_sulc_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.white'
    r_sulc_mris_convert.inputs.out_file = 'R.sulc.native.shape.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}
    
    r_sulc_set_structure_imports = ['import os', 'import subprocess']
    r_sulc_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=r_sulc_set_structure_imports),
                                name='r_sulc_set_structure')
                                
    r_sulc_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(r_sulc_mris_convert,'converted',r_sulc_set_structure,'in_file')   
    
    
    #wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    r_sulc_metric_maths_imports = ['import os', 'import subprocess']
    r_sulc_metric_maths = pe.Node(util.Function(input_names=['in_file', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_3,
                                              imports=r_sulc_metric_maths_imports),
                                name='r_sulc_metric_maths')
                                
    r_sulc_metric_maths.inputs.outfile_name = 'sulc.native.shape.gii'
    wf.connect(r_sulc_set_structure,'out_file',r_sulc_metric_maths,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
    r_sulc_set_maps_imports = ['import os', 'import subprocess']
    r_sulc_set_maps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],                                              function=run_set_maps,
                                              imports=r_sulc_set_maps_imports),
                                name='r_sulc_set_maps')
                                
    r_sulc_set_maps.inputs.map_name = 'subject_R_Sulc'
    wf.connect(r_sulc_metric_maths,'out_file',r_sulc_set_maps,'in_file')   
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
    r_sulc_metric_palette_imports = ['import os', 'import subprocess']
    r_sulc_metric_palette = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=r_sulc_metric_palette_imports),
                                name='r_sulc_metric_palette')
                                
    
    wf.connect(r_sulc_set_maps,'out_file',r_sulc_metric_palette,'in_file') 
    
    
    # For thickness
    #mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    r_thickness_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_thickness_mris_convert')
    r_thickness_mris_convert.inputs.scalarcurv_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.thickness'
    r_thickness_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.white'
    r_thickness_mris_convert.inputs.out_file = 'R.thickness.native.shape.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}
    r_thickness_set_structure_imports = ['import os', 'import subprocess']
    r_thickness_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=r_thickness_set_structure_imports),
                                name='r_thickness_set_structure')
                                
    r_thickness_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(r_thickness_mris_convert,'converted',r_thickness_set_structure,'in_file')   
    
    #wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    r_thickness_metric_maths_imports = ['import os', 'import subprocess']
    r_thickness_metric_maths = pe.Node(util.Function(input_names=['in_file', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_3,
                                              imports=r_thickness_metric_maths_imports),
                                name='r_thickness_metric_maths')
                                
    r_thickness_metric_maths.inputs.outfile_name = 'thickness.native.shape.gii'
    wf.connect(r_thickness_set_structure,'out_file',r_thickness_metric_maths,'in_file')   
    
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
    r_thickness_set_maps_imports = ['import os', 'import subprocess']
    r_thickness_set_maps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_maps,
                                              imports=r_thickness_set_maps_imports),
                                name='r_thickness_set_maps')
                                
    r_thickness_set_maps.inputs.map_name = 'subject_R_Thickness'
    wf.connect(r_thickness_metric_maths,'out_file',r_thickness_set_maps,'in_file')   
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
    r_thickness_metric_palette_imports = ['import os', 'import subprocess']
    r_thickness_metric_palette = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=r_thickness_metric_palette_imports),
                                name='r_thickness_metric_palette')
                                
    
    wf.connect(r_thickness_set_maps,'out_file',r_thickness_metric_palette,'in_file') 
    
    
    # For curvature
    #mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    r_curv_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_curv_mris_convert')
    r_curv_mris_convert.inputs.scalarcurv_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.curv'
    r_curv_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.white'
    r_curv_mris_convert.inputs.out_file = 'R.curvature.native.shape.gii'
    
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}

    r_curv_set_structure_imports = ['import os', 'import subprocess']
    r_curv_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=r_curv_set_structure_imports),
                                name='r_curv_set_structure')
                                
    r_curv_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(r_curv_mris_convert,'converted',r_curv_set_structure,'in_file')   
    
    #wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    r_curv_metric_maths_imports = ['import os', 'import subprocess']
    r_curv_metric_maths = pe.Node(util.Function(input_names=['in_file', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_3,
                                              imports=r_curv_metric_maths_imports),
                                name='r_curv_metric_maths')
                                
    r_curv_metric_maths.inputs.outfile_name = 'curv.native.shape.gii'
    wf.connect(r_curv_set_structure,'out_file',r_curv_metric_maths,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
    r_curv_set_maps_imports = ['import os', 'import subprocess']
    r_curv_set_maps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_maps,
                                              imports=r_curv_set_maps_imports),
                                name='r_curv_set_maps')
                                
    r_curv_set_maps.inputs.map_name = 'subject_R_Curvature'
    wf.connect(r_curv_metric_maths,'out_file',r_curv_set_maps,'in_file')   
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
    
    r_curv_metric_palette_imports = ['import os', 'import subprocess']
    r_curv_metric_palette = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=r_curv_metric_palette_imports),
                                name='r_curv_metric_palette')
                                
    
    wf.connect(r_curv_set_maps,'out_file',r_curv_metric_palette,'in_file') 
    
       ####Right hemisphere ###
   
    #wb_command -metric-math "abs(thickness)" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -var thickness  "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
    r_thickness_metric_maths_imports_1 = ['import os', 'import subprocess']
    r_thickness_metric_maths_1 = pe.Node(util.Function(input_names=['in_file','file_name'],
                                             output_names=['out_file'],
                                             function=run_metric_maths_2,
                                             imports=r_thickness_metric_maths_imports_1),
                               name='r_thickness_metric_maths_1')
                               
    r_thickness_metric_maths_1.inputs.in_file = '/data3/cnl/tgeorge/pfreesurfer-bash-runs/dcan_postfreesurfer_wf/dcan_postfreesurfer_block8/r_thickness_mris_convert/R.thickness.native.shape.gii'  
    r_thickness_metric_maths_1.inputs.file_name = 'R.thickness.native.shape.gii'
    
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
    r_thickness_metric_palette_imports_2 = ['import os', 'import subprocess']
    r_thickness_metric_palette_2 = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=r_thickness_metric_palette_imports_2),
                                name='r_thickness_metric_palette_2')
                                
    wf.connect(r_thickness_metric_maths_1, 'out_file', r_thickness_metric_palette_2 , 'in_file')

    
    #wb_command -metric-math "thickness > 0" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -var thickness "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
    r_thickness_metric_maths_imports_2 = ['import os', 'import subprocess']
    r_thickness_metric_maths_2 = pe.Node(util.Function(input_names=['in_file', 'file_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_1,
                                              imports=r_thickness_metric_maths_imports_2),
                                name='r_thickness_metric_maths_2')
    r_thickness_metric_maths_2.inputs.file_name = "R.roi.native.shape.gii"                        
    wf.connect(r_thickness_metric_palette, 'out_file', r_thickness_metric_maths_2, 'in_file')
    
    
    #wb_command -metric-fill-holes "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    r_thickness_metric_fillholes_imports = ['import os', 'import subprocess']
    r_thickness_metric_fillholes = pe.Node(util.Function(input_names=['surface_file','metric_infile','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_fillholes,
                                              imports=r_thickness_metric_fillholes_imports),
                                name='r_thickness_metric_fillholes')
                                
    r_thickness_metric_fillholes.inputs.outfile_name = "R.roi.native.shape.gii"   
    
    wf.connect(r_thickness_metric_maths_2,'out_file',r_thickness_metric_fillholes, 'metric_infile')
    wf.connect(Atlasspace_r_surface_avg, 'out_file',r_thickness_metric_fillholes, 'surface_file')
    
    #wb_command -metric-remove-islands "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
   
   
    r_thickness_metric_removeislands_imports = ['import os', 'import subprocess']
    r_thickness_metric_removeislands = pe.Node(util.Function(input_names=['surface_file','metric_infile', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_removeislands,
                                              imports=r_thickness_metric_removeislands_imports),
                                name='r_thickness_metric_removeislands')
                                
    r_thickness_metric_removeislands.inputs.outfile_name = 'R.roi.native.shape.gii'
    wf.connect(r_thickness_metric_fillholes,'out_file',r_thickness_metric_removeislands, 'metric_infile')
    wf.connect(Atlasspace_r_surface_avg,'out_file',r_thickness_metric_removeislands, 'surface_file')
  
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -map 1 "$Subject"_"$Hemisphere"_ROI
    r_thickness_metric_setmaps_imports = ['import os', 'import subprocess']
    r_thickness_metric_setmaps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_setmaps,
                                              imports=r_thickness_metric_setmaps_imports),
                                name='r_thickness_metric_setmaps')
                                
    r_thickness_metric_setmaps.inputs.map_name = 'R_ROI'
    wf.connect(r_thickness_metric_removeislands,'out_file',r_thickness_metric_setmaps, 'in_file')
   
   
    #wb_command -metric-dilate "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -nearest
    r_thickness_metric_dilate_imports = ['import os', 'import subprocess']
    r_thickness_metric_dilate = pe.Node(util.Function(input_names=['in_file','surface_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_dilate,
                                              imports=r_thickness_metric_dilate_imports),
                                name='r_thickness_metric_dilate')
                                
    r_thickness_metric_dilate.inputs.outfile_name = 'R.thickness.native.shape.gii'
    wf.connect(r_thickness_metric_palette,'out_file',r_thickness_metric_dilate, 'in_file')
    wf.connect(Atlasspace_r_surface_avg,'out_file',r_thickness_metric_dilate,'surface_file')
   
   
    #wb_command -metric-dilate "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii -nearest
   
    #wb_command -metric-dilate "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -nearest
    r_curvature_metric_dilate_imports = ['import os', 'import subprocess']
    r_curvature_metric_dilate = pe.Node(util.Function(input_names=['in_file','surface_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_dilate,
                                              imports=r_curvature_metric_dilate_imports),
                                name='r_curvature_metric_dilate')
                                
    r_curvature_metric_dilate.inputs.outfile_name = 'curvature.native.shape.gii'
    wf.connect(r_curv_mris_convert,'converted',r_curvature_metric_dilate, 'in_file')
    wf.connect(Atlasspace_r_surface_avg, 'out_file', r_curvature_metric_dilate, 'surface_file')
    
####left hemesphere 

  #wb_command -metric-math "abs(thickness)" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -var thickness  "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
    l_thickness_metric_maths_imports_1 = ['import os', 'import subprocess']
    l_thickness_metric_maths_1 = pe.Node(util.Function(input_names=['in_file','file_name'],
                                             output_names=['out_file'],
                                             function=run_metric_maths_2,
                                             imports=l_thickness_metric_maths_imports_1),
                               name='l_thickness_metric_maths_1')
                               
    
    l_thickness_metric_maths_1.inputs.file_name = 'L.thickness.native.shape.gii'
    wf.connect(l_thickness_mris_convert, 'converted',l_thickness_metric_maths_1, 'in_file')
    
    #wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
    l_thickness_metric_palette_imports = ['import os', 'import subprocess']
    l_thickness_metric_palette = pe.Node(util.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=run_metric_palette,
                                              imports=l_thickness_metric_palette_imports),
                                name='l_thickness_metric_palette')
                                
    wf.connect(l_thickness_metric_maths_1, 'out_file', l_thickness_metric_palette , 'in_file')

    
    #wb_command -metric-math "thickness > 0" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -var thickness "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
    l_thickness_metric_maths_imports_2 = ['import os', 'import subprocess']
    l_thickness_metric_maths_2 = pe.Node(util.Function(input_names=['in_file', 'file_name'],
                                              output_names=['out_file'],
                                              function=run_metric_maths_1,
                                              imports=l_thickness_metric_maths_imports_2),
                                name='l_thickness_metric_maths_2')
                                
    l_thickness_metric_maths_2.inputs.file_name = 'L.roi.native.shape.gii'
    wf.connect(l_thickness_metric_palette, 'out_file', l_thickness_metric_maths_2, 'in_file')
    
    
    #wb_command -metric-fill-holes "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    l_thickness_metric_fillholes_imports = ['import os', 'import subprocess']
    l_thickness_metric_fillholes = pe.Node(util.Function(input_names=['surface_file','metric_infile','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_fillholes,
                                              imports=l_thickness_metric_fillholes_imports),
                                name='l_thickness_metric_fillholes')
                                
    l_thickness_metric_fillholes.inputs.outfile_name = 'L.roi.native.shape.gii'
    wf.connect(l_thickness_metric_maths_2,'out_file',l_thickness_metric_fillholes, 'metric_infile')
    wf.connect(Atlasspace_l_surface_avg,'out_file',l_thickness_metric_fillholes,'surface_file')
    
    #wb_command -metric-remove-islands "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
   
   
    l_thickness_metric_removeislands_imports = ['import os', 'import subprocess']
    l_thickness_metric_removeislands = pe.Node(util.Function(input_names=['surface_file','metric_infile', 'outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_removeislands,
                                              imports=l_thickness_metric_removeislands_imports),
                                name='l_thickness_metric_removeislands')
                                
    l_thickness_metric_removeislands.inputs.outfile_name = 'L.roi.native.shape.gii'
   
    wf.connect(l_thickness_metric_fillholes,'out_file',l_thickness_metric_removeislands, 'metric_infile')
    wf.connect(Atlasspace_l_surface_avg, 'out_file',l_thickness_metric_removeislands, 'surface_file')
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -map 1 "$Subject"_"$Hemisphere"_ROI
    l_thickness_metric_setmaps_imports = ['import os', 'import subprocess']
    l_thickness_metric_setmaps = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_setmaps,
                                              imports=l_thickness_metric_setmaps_imports),
                                name='l_thickness_metric_setmaps')
                                
    l_thickness_metric_setmaps.inputs.map_name = 'R_ROI'
    wf.connect(l_thickness_metric_removeislands,'out_file',l_thickness_metric_setmaps, 'in_file')
   
   
    #wb_command -metric-dilate "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -nearest
    l_thickness_metric_dilate_imports = ['import os', 'import subprocess']
    l_thickness_metric_dilate = pe.Node(util.Function(input_names=['in_file','surface_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_dilate,
                                              imports=l_thickness_metric_dilate_imports),
                                name='l_thickness_metric_dilate')
                                
    l_thickness_metric_dilate.inputs.outfile_name = 'L.thickness.native.shape.gii'
    wf.connect(l_thickness_metric_palette,'out_file',l_thickness_metric_dilate, 'in_file')
    wf.connect(Atlasspace_l_surface_avg, 'out_file',l_thickness_metric_dilate,'surface_file')
   
   
    #wb_command -metric-dilate "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii -nearest
   
    l_curvature_metric_dilate_imports = ['import os', 'import subprocess']
    l_curvature_metric_dilate = pe.Node(util.Function(input_names=['in_file','surface_file','outfile_name'],
                                              output_names=['out_file'],
                                              function=run_metric_dilate,
                                              imports=l_curvature_metric_dilate_imports),
                                name='l_curvature_metric_dilate')
                                
    l_curvature_metric_dilate.inputs.outfile_name = 'curvature.native.shape.gii'
   
    wf.connect(Atlasspace_l_surface_avg, 'out_file',l_curvature_metric_dilate, 'surface_file')
    wf.connect(l_curv_mris_convert, 'converted', l_curvature_metric_dilate, 'in_file')


    # Right aparc

    #mris_convert --annot "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii

 
    r_aparc_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_aparc_mris_convert')
    r_aparc_mris_convert.inputs.annot_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/label/rh.aparc.annot'
    r_aparc_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.white'
    r_aparc_mris_convert.inputs.out_file = 'R.aparc.native.label.gii'
    
	
	
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii $Structure
     
    r_aparc_set_structure_imports = ['import os', 'import subprocess']
    r_aparc_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=r_aparc_set_structure_imports),
                                name='r_aparc_set_structure')
                                
    r_aparc_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(r_aparc_mris_convert,'converted',r_aparc_set_structure,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii -map 1 "$Subject"_"$Hemisphere"_"$Map"
    
    r_aparc_set_mapnames_imports = ['import os', 'import subprocess']
    r_aparc_set_mapnames = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_mapnames,
                                              imports=r_aparc_set_mapnames_imports),
                                name='r_aparc_set_mapnames')
                                
    r_aparc_set_mapnames.inputs.map_name = 'R_aparc'
    wf.connect(r_aparc_set_structure,'out_file',r_aparc_set_mapnames,'in_file')  
    
    
  
    
    #wb_command -gifti-label-add-prefix "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii "${Hemisphere}_" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
    r_aparc_giftilabels_imports = ['import os', 'import subprocess']
    r_aparc_giftilabels = pe.Node(util.Function(input_names=['in_file','prefix','file_name'],
                                              output_names=['out_file'],
                                              function=run_set_giftilabels,
                                              imports=r_aparc_giftilabels_imports),
                                name='r_aparc_giftilabels')
                                
    r_aparc_giftilabels.inputs.prefix = 'R_'
    r_aparc_giftilabels.inputs.file_name = 'R.aparc.native.label.gii'
    wf.connect(r_aparc_set_mapnames,'out_file',r_aparc_giftilabels,'in_file')  
    
    
       # Right aparc.a2009s

    #mris_convert --annot "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii

 
    r_aparc2009s_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='r_aparc2009s_mris_convert')
    r_aparc2009s_mris_convert.inputs.annot_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/label/rh.aparc.a2009s.annot'
    r_aparc2009s_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/rh.white'
    r_aparc2009s_mris_convert.inputs.out_file = 'R.aparc.a2009s.native.label.gii'
    
	
	
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii $Structure
     
    r_aparc2009s_set_structure_imports = ['import os', 'import subprocess']
    r_aparc2009s_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=r_aparc2009s_set_structure_imports),
                                name='r_aparc2009s_set_structure')
                                
    r_aparc2009s_set_structure.inputs.structure = "CORTEX_RIGHT"
    wf.connect(r_aparc2009s_mris_convert,'converted',r_aparc2009s_set_structure,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii -map 1 "$Subject"_"$Hemisphere"_"$Map"
    
    r_aparc2009s_set_mapnames_imports = ['import os', 'import subprocess']
    r_aparc2009s_set_mapnames = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_mapnames,
                                              imports=r_aparc2009s_set_mapnames_imports),
                                name='r_aparc2009s_set_mapnames')
                                
    r_aparc2009s_set_mapnames.inputs.map_name = 'R_aparc2009s'
    wf.connect(r_aparc2009s_set_structure,'out_file',r_aparc2009s_set_mapnames,'in_file')  
    
    
  
    
    #wb_command -gifti-label-add-prefix "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii "${Hemisphere}_" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
    r_aparc2009s_giftilabels_imports = ['import os', 'import subprocess']
    r_aparc2009s_giftilabels = pe.Node(util.Function(input_names=['in_file','prefix', 'file_name'],
                                              output_names=['out_file'],
                                              function=run_set_giftilabels,
                                              imports=r_aparc2009s_giftilabels_imports),
                                name='r_aparc2009s_giftilabels')
                                
    r_aparc2009s_giftilabels.inputs.prefix = 'R_'
    r_aparc2009s_giftilabels.inputs.file_name = 'R.aparc.a2009s.native.label.gii'
    wf.connect(r_aparc2009s_set_mapnames,'out_file',r_aparc2009s_giftilabels,'in_file')  
    
    
    
     # Left aparc

    #mris_convert --annot "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii

 
    l_aparc_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_aparc_mris_convert')
    l_aparc_mris_convert.inputs.annot_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/label/lh.aparc.annot'
    l_aparc_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.white'
    l_aparc_mris_convert.inputs.out_file = 'L.aparc.native.label.gii'
    
	
	
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii $Structure
     
    l_aparc_set_structure_imports = ['import os', 'import subprocess']
    l_aparc_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_1,
                                              imports=l_aparc_set_structure_imports),
                                name='l_aparc_set_structure')
                                
    l_aparc_set_structure.inputs.structure = "CORTEX_LEFT"
    wf.connect(l_aparc_mris_convert,'converted',l_aparc_set_structure,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii -map 1 "$Subject"_"$Hemisphere"_"$Map"
    
    l_aparc_set_mapnames_imports = ['import os', 'import subprocess']
    l_aparc_set_mapnames = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_mapnames,
                                              imports=l_aparc_set_mapnames_imports),
                                name='l_aparc_set_mapnames')
                                
    l_aparc_set_mapnames.inputs.map_name = 'L_aparc'
    wf.connect(l_aparc_set_structure,'out_file',l_aparc_set_mapnames,'in_file')  
    
    
  
    
    #wb_command -gifti-label-add-prefix "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii "${Hemisphere}_" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
    l_aparc_giftilabels_imports = ['import os', 'import subprocess']
    l_aparc_giftilabels = pe.Node(util.Function(input_names=['in_file','prefix','file_name'],
                                              output_names=['out_file'],
                                              function=run_set_giftilabels,
                                              imports=l_aparc_giftilabels_imports),
                                name='l_aparc_giftilabels')
    
    l_aparc_giftilabels.inputs.file_name = 'L.aparc.native.label.gii'                            
    l_aparc_giftilabels.inputs.prefix = 'L_'
    wf.connect(l_aparc_set_mapnames,'out_file',l_aparc_giftilabels,'in_file')  
    
    
       # Left aparc.a2009s

    #mris_convert --annot "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii

 
    l_aparc2009s_mris_convert = pe.Node(interface=freesurfer.MRIsConvert(),name='l_aparc2009s_mris_convert')
    l_aparc2009s_mris_convert.inputs.annot_file = '/data3/cnl/tgeorge/postfreesurfer/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/label/lh.aparc.a2009s.annot'
    l_aparc2009s_mris_convert.inputs.in_file = '/home/tgeorge/postfreesurfer/postfs-output/working/cpac_sub-0050952_ses-1/anat_preproc_freesurfer_52/anat_freesurfer/recon_all/surf/lh.white'
    l_aparc2009s_mris_convert.inputs.out_file = 'L.aparc.a2009s.native.label.gii'
    
	
	
    #wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii $Structure
     
    l_aparc2009s_set_structure_imports = ['import os', 'import subprocess']
    l_aparc2009s_set_structure = pe.Node(util.Function(input_names=['in_file','structure'],
                                              output_names=['out_file'],
                                              function=run_set_structure_2,
                                              imports=l_aparc2009s_set_structure_imports),
                                name='l_aparc2009s_set_structure')
                                
    l_aparc2009s_set_structure.inputs.structure = "CORTEX_LEFT"
    wf.connect(l_aparc2009s_mris_convert,'converted',l_aparc2009s_set_structure,'in_file')   
    
    #wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii -map 1 "$Subject"_"$Hemisphere"_"$Map"
    
    l_aparc2009s_set_mapnames_imports = ['import os', 'import subprocess']
    l_aparc2009s_set_mapnames = pe.Node(util.Function(input_names=['in_file','map_name'],
                                              output_names=['out_file'],
                                              function=run_set_mapnames,
                                              imports=l_aparc2009s_set_mapnames_imports),
                                name='l_aparc2009s_set_mapnames')
                                
    l_aparc2009s_set_mapnames.inputs.map_name = 'L_aparc2009s'
    wf.connect(l_aparc2009s_set_structure,'out_file',l_aparc2009s_set_mapnames,'in_file')  
    
    
  
    
    #wb_command -gifti-label-add-prefix "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii "${Hemisphere}_" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
    l_aparc2009s_giftilabels_imports = ['import os', 'import subprocess']
    l_aparc2009s_giftilabels = pe.Node(util.Function(input_names=['in_file','prefix','file_name'],
                                              output_names=['out_file'],
                                              function=run_set_giftilabels,
                                              imports=l_aparc2009s_giftilabels_imports),
                                name='l_aparc2009s_giftilabels')
                                
    l_aparc2009s_giftilabels.inputs.prefix = 'L_'
    l_aparc2009s_giftilabels.inputs.file_name = 'L.aparc.a2009s.native.label.gii'
    wf.connect(l_aparc2009s_set_mapnames,'out_file',l_aparc2009s_giftilabels,'in_file')  

    
    ds = pe.Node(DataSink(), name='ds')
    ds.inputs.parameterization= True
    ds.inputs.base_directory= '/data3/cnl/tgeorge/pfreesurfer-bash-runs'
    ds.inputs.container = 'output_block_combined'
    ds.inputs.container = 'output_specfile'
    wf.connect(l_aparc2009s_giftilabels,'out_file', ds,'output_block_combined')
    wf.connect(r_sphere_add_spec_file,'out_file',ds,'output_specfile')
    wf.run()
    
    
    
   

run_workflow()     

