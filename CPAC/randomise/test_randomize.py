import os
import glob
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
from nipype.interfaces.fsl import ImageStats


def run_randomize(inputs,output_dir=None,run=True):
   import pipeline
   randomise_workflow = pe.Workflow(name = 'preproc')
   if output_dir == None:
     output_dir = '/home/nrajamani'
        
   workflow_dir = os.path.join(output_dir,"randomise_results")
   randomise_workflow.base_dir = workflow_dir
   # taken from QAP files 
   #resource_pool = {}
   
   num_of_cores = 1
   #resource_pool({'epireg': (warp_nipype2.warp_nipype, 'outputspec.epireg')})
   t_node = pipeline.create_randomise()
   t_node.inputs.inputspec.subjects=  '/home/nrajamani/randomise_files/single_grp_age_reho_to_standard_zstd_smooth_merged.nii.gz'
   t_node.inputs.inputspec.design_matrix_file= '/home/nrajamani/randomise_files/single_grp_age.mat'
   t_node.inputs.inputspec.constrast_file= '/home/nrajamani/randomise_files/single_grp_age.con'
   #t_node.inputs.inputspec.f_constrast_file= '/home/nrajamani/'
   t_node.inputs.inputspec.permutations=5000
   #t_node.inputs.inputspec.mask = '/home/nrajamani/group_analysis_results_benchmark-FNIRT/group_model_single_grp_age/alff_to_standard_zstd/rest_run-1/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic0.gm0.compcor1.csf1/_hp_0.01/_lp_0.1/_target_angle_deg_90/alff_warp_maths/model_files/single_grp_age_alff_to_standard_zstd_merged_mask.nii.gz'# 


   dataSink = pe.Node(nio.DataSink(), name='dataSink_file')
   dataSink.inputs.base_directory = workflow_dir
   #node, out_file = resource_pool["epireg"]
   #warp_workflow.connect(t_node,'outputspec.roi_file',dataSink,'roi_file')
   randomise_workflow.connect(t_node,'outputspec.index_file',dataSink, 'index_file')
   #randomise_workflow.connect(t_node,'outputspec.thresh_out',dataSink,'threshold_file')
   randomise_workflow.connect(t_node,'outputspec.localmax_txt_file',dataSink,'localmax_txt_file')
   randomise_workflow.connect(t_node,'outputspec.localmax_vol_file',dataSink,'localmax_vol_file')
   randomise_workflow.connect(t_node,'outputspec.max_file',dataSink,'max_file')
   randomise_workflow.connect(t_node,'outputspec.mean_file',dataSink,'mean_file')
   randomise_workflow.connect(t_node,'outputspec.pval_file',dataSink,'pval_file')
   randomise_workflow.connect(t_node,'outputspec.size_file',dataSink,'size_file')
   randomise_workflow.connect(t_node,'outputspec.tstat_files',dataSink,'tstat_files')
   randomise_workflow.connect(t_node,'outputspec.t_corrected_p_files',dataSink,'t_corrected_p_files')
   if run == True:
       randomise_workflow.run(plugin='MultiProc', plugin_args ={'n_procs': num_of_cores})
       #outpath = glob.glob(os.path.join(workflow_dir, "EPI_DistCorr","*"))[0]
       #return outpath
   else:
       return randomise_workflow, randomise_workflow.base_dir

run_randomize(['subjects','design_matrix_file','constrast_file','permutations'],output_dir=None,run=True)
