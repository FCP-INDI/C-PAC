import os
import glob
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
from nipype.interfaces.fsl import ImageStats


def run_randomize(inputs,output_dir=None,run=True):
   import pipeline
   randomise_workflow = pe.Workflow(name = 'preproc')
   if output_dir == None:
     output_dir = ''
        
   workflow_dir = os.path.join(output_dir,"randomise_results")
   randomise_workflow.base_dir = workflow_dir
  
   #resource_pool = {}
   
   num_of_cores = 1
   
   t_node = pipeline.create_randomise()
   t_node.inputs.inputspec.subjects=  ''
   t_node.inputs.inputspec.design_matrix_file= ''
   t_node.inputs.inputspec.constrast_file= ''
   #t_node.inputs.inputspec.f_constrast_file= ''
   t_node.inputs.inputspec.permutations=5000
   #t_node.inputs.inputspec.mask = ''# 


   dataSink = pe.Node(nio.DataSink(), name='dataSink_file')
   dataSink.inputs.base_directory = workflow_dir
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
   else:
       return randomise_workflow, randomise_workflow.base_dir

run_randomize(['subjects','design_matrix_file','constrast_file','permutations'],output_dir=None,run=True)
