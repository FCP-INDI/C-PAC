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
#         cpac utils data_config build <data settings file>
#     cpac utils pipe_config
#         cpac utils pipe_config new_template


@click.group()
def main():
    pass


@main.command()
def version():
    import CPAC
    print('C-PAC version: {0}'.format(CPAC.__version__))


@main.command()
def gui():
    import CPAC.GUI
    CPAC.GUI.run()


@main.command()
@click.argument('data_config')
@click.option('--pipe_config')
def run(data_config, pipe_config=None):
    if not pipe_config:
        import os
        import pkg_resources as p
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_template.yml"))
    if data_config == 'ADHD200':
        import os
        import pkg_resources as p
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-ADHD200.yml"))
    if data_config == 'ADHD200_2':
        import os
        import pkg_resources as p
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-ADHD200_only2.yml"))
    if data_config == 'ABIDE':
        import os
        import pkg_resources as p
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-ABIDE.yml"))
    if data_config == 'NKI-RS':
        import os
        import pkg_resources as p
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-NKI-RocklandSample.yml"))

    import CPAC.pipeline.cpac_runner as cpac_runner
    cpac_runner.run(pipe_config, data_config)


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
@click.option('--pipe_output_dir', default=None)
def run(pipe_config, pipe_output_dir=None):
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
@group.command()
@click.argument('pipe_config')
def basc(pipe_config):
    import CPAC.pipeline.cpac_group_runner as cpac_group_runner
    cpac_group_runner.run_basc(pipe_config)


@group.command(name="mdmr")
@click.argument("pipeline_config", type=click.Path(exists=True))
def group_mdmr(pipeline_config):
    from CPAC.pipeline.cpac_group_runner import run_cwas
    
    run_cwas(pipeline_config)
    

# Utilities
@main.group()
def utils():
    pass


@utils.group()
def data_config():
    pass


@data_config.command()
def new_template():
    from CPAC.utils.build_data_config import util_copy_template
    util_copy_template()


@data_config.command()
@click.argument('data_settings_file')
def build(data_settings_file):
    from CPAC.utils.build_data_config import run
    run(data_settings_file)


@utils.group()
def pipe_config():
    pass


@pipe_config.command()
def new_template():
    from CPAC.utils.build_data_config import util_copy_template
    util_copy_template('pipeline_config')


@utils.group()
def test():
    pass


@test.command()
def run_suite():
    import os
    import pkg_resources as p
    import CPAC.pipeline.cpac_runner as cpac_runner

    test_config_dir = \
        p.resource_filename("CPAC",
                            os.path.join("resources",
                                         "configs", "test_configs"))

    data_test = \
        p.resource_filename("CPAC",
                            os.path.join("resources",
                                         "configs", "test_configs",
                                         "data-test_S3-ADHD200_1.yml"))

    data_test_no_scan_param = \
        p.resource_filename("CPAC",
                            os.path.join("resources",
                                         "configs", "test_configs",
                                         "data-test_S3-ADHD200_no-params.yml"))

    data_test_fmap = \
        p.resource_filename("CPAC",
                            os.path.join("resources",
                                         "configs", "test_configs",
                                         "data-test_S3-NKI-RS_fmap.yml"))

    no_params = False
    for config_file in os.listdir(test_config_dir):
        if 'pipe-test' in config_file:
            pipe = os.path.join(test_config_dir, config_file)
            if 'DistCorr' in pipe:
                data = data_test_fmap
            elif not no_params:
                data = data_test_no_scan_param
                no_params = True
            else:
                data = data_test
            # run
            cpac_runner.run(pipe, data)


if __name__ == "__main__":
    main()
