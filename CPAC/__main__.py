#!/usr/bin/env python

import CPAC
import click

# CLI tree
#
# cpac
#     cpac gui
# cpac individual
#     cpac --pipe <pipeline config> --data <data config>
#     cpac --data <data config>  (runs default pipeline)
#     cpac --benchmark
# cpac group
#     cpac group feat
#         cpac group feat <pipeline config>
#         cpac group feat load_preset <args...>
#     cpac group basc
#         cpac group basc <pipeline config>
#     cpac group mdmr
#         cpac group mdmr <pipeline config>
# cpac utils
#     cpac utils data_config
#         cpac utils data_config new_template
#         cpac utils data_config build_data_config <data settings file>
#     cpac utils pipe_config
#         cpac utils pipe_config new_template
#     cpac utils feat_load_preset


@click.group()
def main():
    pass


@main.command()
def gui():
    import CPAC.GUI
    CPAC.GUI.run()


@main.group()
def individual():
    pass


@individual.command()
@click.argument('pipe_config')
@click.argument('data_config')
def run(pipe_config, data_config):
    import CPAC.pipeline.cpac_runner as cpac_runner
    cpac_runner.run(pipe_config, data_config)


@main.group()
@click.pass_context
def participant(ctx):
    ctx.forward(individual)


# Group analysis
@main.group()
def group():
    pass


# Group analysis - FSL FEAT
@group.group()
def feat():
    pass


@feat.command()
@click.argument('pipe_config')
@click.argument('pipe_output_dir')
def run(pipe_config, pipe_output_dir):
    import CPAC.pipeline.cpac_group_runner as cpac_group_runner
    cpac_group_runner.run_feat(pipe_config, pipe_output_dir)


@feat.group()
def load_preset():
    pass

@load_preset.command()
@click.argument('group_participants')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('model_name')
@click.option('--output_dir', default=None)
def single_grp_avg(group_participants, z_thresh, p_thresh, model_name,
                   output_dir=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(group_participants, 'all', z_thresh, p_thresh,
                                'single_grp', output_dir=output_dir,
                                model_name=model_name)

@load_preset.command()
@click.argument('group_participants')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('pheno_file')
@click.argument('pheno_sub')
@click.argument('covariate')
@click.argument('model_name')
@click.option('--output_dir', default=None)
def single_grp_cov(group_participants, z_thresh, p_thresh, pheno_file,
                   pheno_sub, covariate, model_name, output_dir=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(group_participants, 'all', z_thresh, p_thresh,
                                'single_grp_cov', pheno_file=pheno_file,
                                pheno_sub_label=pheno_sub,
                                covariate=covariate, output_dir=output_dir,
                                model_name=model_name)

@load_preset.command()
@click.argument('group_participants')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('pheno_file')
@click.argument('pheno_sub')
@click.argument('covariate')
@click.argument('model_name')
@click.option('--output_dir', default=None)
def unpaired_two(group_participants, z_thresh, p_thresh, pheno_file,
                 pheno_sub, covariate, model_name, output_dir=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(group_participants, 'all', z_thresh, p_thresh,
                                'unpaired_two', pheno_file=pheno_file,
                                pheno_sub_label=pheno_sub,
                                covariate=covariate, output_dir=output_dir,
                                model_name=model_name)

@load_preset.command()
@click.argument('group_participants')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('conditions')
@click.argument('condition_type')
@click.argument('model_name')
@click.option('--output_dir', default=None)
def paired_two(group_participants, z_thresh, p_thresh, conditions,
               condition_type, model_name, output_dir=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(group_participants, 'all', z_thresh, p_thresh,
                                'paired_two', covariate=conditions,
                                condition_type=condition_type,
                                output_dir=output_dir, model_name=model_name)


@load_preset.command()
@click.argument('group_participants')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('conditions')
@click.argument('condition_type')
@click.argument('model_name')
@click.option('--output_dir', default=None)
def tripled_two(group_participants, z_thresh, p_thresh, conditions,
               condition_type, model_name, output_dir=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(group_participants, 'all', z_thresh, p_thresh,
                                'tripled_two', covariate=conditions,
                                condition_type=condition_type,
                                output_dir=output_dir, model_name=model_name)


# Group analysis - PyBASC
@group.group()
def basc():
    pass


@basc.command()
@click.argument('pipe_config')
def run(pipe_config):
    import CPAC.pipeline.cpac_group_runner as cpac_group_runner
    cpac_group_runner.run_basc(pipe_config)


@basc.command()
@click.argument('pipeline_dir')
@click.argument('roi_file')
@click.option('--roi_file_two', default=None)
@click.option('--ref_file', default=None)
@click.option('--output_dir', default=None)
@click.option('--output_size', default=800)
@click.option('--basc_proc', default=2)
@click.option('--basc_memory', default=4)
@click.option('--scan', default=None)
def quickrun(pipeline_dir, roi_file, roi_file_two=None, ref_file=None,
             output_size=800, output_dir=None, basc_proc=2, basc_memory=4,
             scan=None):
    import CPAC.pipeline.cpac_group_runner as cpac_group_runner
    cpac_group_runner.run_basc_quickrun(pipeline_dir, roi_file, roi_file_two,
                                        ref_file, output_size, output_dir,
                                        basc_proc, basc_memory, scan=scan)

    
@group.command(name="mdmr")
@click.argument("pipeline_config", type=click.Path(exists=True))
def group_mdmr(pipeline_config):
    from CPAC.pipeline.cpac_group_runner import run_cwas
    
    run_cwas(pipeline_config)
    

# Utilities
@main.group()
def utils():
    pass


@utils.command(name="data_config")
def utils_data_config():
    pass


@utils.group()
def data_config():
    pass


@data_config.command()
def new_template():
    from CPAC.utils.build_data_config import util_copy_data_settings_template
    util_copy_data_settings_template()


@data_config.command()
@click.argument('data_settings_file')
def build_data_config(data_settings_file):
    from CPAC.utils.build_data_config import run
    run(data_settings_file)


@utils.group()
def pipe_config():
    pass


@pipe_config.command()
def new_template():
    # TODO: this
    pass


if __name__ == "__main__":
    main()
