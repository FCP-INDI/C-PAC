# CPAC/pipeline/cluter_templates.py
#
# Author: Daniel Clark, 2016

# Start task ID string
start_taskid_str = 'echo "Start - TASKID " %(env_arr_idx)s " : " $(date)'

# End task ID string
end_taskid_str = 'echo "End - TASKID " %(env_arr_idx)s " : " $(date)'

# Run CPAC via python -c print command
python_cpac_str = 'python -c "from CPAC.pipeline.cpac_pipeline import run; '\
                  'run(\'%(config_file)s\', \'%(subject_list_file)s\', '\
                      '%(env_arr_idx)s, \'%(strategies_file)s\', '\
                      '\'%(pipeline_name)s\', plugin=\'ResourceMultiProc\', '\
                      'plugin_args=%(plugin_args)s)"'

# SGE template string
pbs_template = \
'''#! %(shell)s
## C-PAC PBS batch file - %(timestamp)s
#PBS -S %(shell)s
#PBS -N cpac_pipeline_%(pipeline_name)s
#PBS -t 1-%(num_subs)d
#PBS -q %(queue)s
#PBS -l nodes=1:ppn=%(cores_per_sub)d
#PBS -A %(user)s
#PBS -V
#PBS -wd %(work_dir)s
'''
# Add in start, python CPAC, end
pbs_template = '\n'.join([pbs_template,
                          start_taskid_str, python_cpac_str, end_taskid_str])

# SGE template string
sge_template = \
'''#! %(shell)s
## C-PAC SGE batch file - %(timestamp)s
#$ -S %(shell)s
#$ -N cpac_pipeline_%(pipeline_name)s
#$ -t 1-%(num_subs)d
#$ -q %(queue)s
#$ -pe %(par_env)s %(cores_per_sub)d
#$ -A %(user)s
#$ -V
#$ -wd %(work_dir)s
##$ -e %(err_log)s
##$ -o %(out_log)s
'''
# Add in start, python CPAC, end
sge_template = '\n'.join([sge_template,
                          start_taskid_str, python_cpac_str, end_taskid_str])

# SLURM template string
slurm_template = \
'''#! %(shell)s
## C-PAC SLURM batch file - %(timestamp)s
#SBATCH --job-name=cpac_pipeline_%(pipeline_name)s
#SBATCH --array=1-%(num_subs)d
#SBATCH --cpus-per-task=%(cores_per_sub)d
#SBATCH --uid=%(user)s
#SBATCH --get-user-env
#SBATCH --workdir=%(work_dir)s
##SBATCH --error=%(err_log)s
##SBATCH --output=%(out_log)s
'''
# Add in start, python CPAC, end
slurm_template = '\n'.join([slurm_template,
                            start_taskid_str, python_cpac_str, end_taskid_str])
