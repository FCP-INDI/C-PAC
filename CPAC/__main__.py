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
#     cpac group isc
#         cpac group isc <pipeline config>
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
@click.option('--num_cores')
def run(data_config, pipe_config=None, num_cores=None):
    if not pipe_config:
        import os
        import pkg_resources as p
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_template.yml"))

    if pipe_config == 'benchmark-ants':
        import os
        import pkg_resources as p
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_benchmark-ANTS.yml"))

    if pipe_config == 'benchmark-fnirt':
        import os
        import pkg_resources as p
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_benchmark-FNIRT.yml"))

    if data_config == 'benchmark-data':
        import os
        import pkg_resources as p
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_cpac_benchmark.yml"))

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
    cpac_runner.run(pipe_config, data_config, num_subs_at_once=num_cores)


# Group analysis
@main.group()
def group():
    pass


# Group analysis - FSL FEAT
@group.group()
def feat():
    pass


@feat.command()
@click.argument('group_config')
def build_models(group_config):
    import CPAC.pipeline.cpac_group_runner as cgr
    cgr.build_feat_models(group_config)


@feat.command()
@click.argument('group_config')
def run(group_config):
    import CPAC.pipeline.cpac_group_runner as cgr
    cgr.run_feat(group_config)


@feat.command()
@click.argument('group_config')
def randomise(group_config):
    import CPAC.pipeline.cpac_group_runner as cgr
    cgr.run_feat(group_config, feat=False)


@feat.group()
def load_preset():
    pass

@load_preset.command()
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('model_name')
@click.option('--output_dir', default=None)
@click.option('--group_participants', default=None)
def single_grp_avg(pipeline_dir, z_thresh, p_thresh, model_name,
                   output_dir=None, group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'single_grp', group_participants,
                                output_dir=output_dir,
                                model_name=model_name)

@load_preset.command()
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('pheno_file')
@click.argument('pheno_sub')
@click.argument('covariate')
@click.argument('model_name')
@click.option('--output_dir', default=None)
@click.option('--group_participants', default=None)
def single_grp_cov(pipeline_dir, z_thresh, p_thresh, pheno_file,
                   pheno_sub, covariate, model_name, output_dir=None,
                   group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'single_grp_cov', group_participants,
                                pheno_file=pheno_file,
                                pheno_sub_label=pheno_sub,
                                covariate=covariate, output_dir=output_dir,
                                model_name=model_name)

@load_preset.command()
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('pheno_file')
@click.argument('pheno_sub')
@click.argument('covariate')
@click.argument('model_name')
@click.option('--output_dir', default=None)
@click.option('--group_participants', default=None)
def unpaired_two(pipeline_dir, z_thresh, p_thresh, pheno_file,
                 pheno_sub, covariate, model_name, output_dir=None,
                 group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'unpaired_two', group_participants,
                                pheno_file=pheno_file,
                                pheno_sub_label=pheno_sub,
                                covariate=covariate, output_dir=output_dir,
                                model_name=model_name)

@load_preset.command()
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('conditions')
@click.argument('condition_type')
@click.argument('model_name')
@click.option('--output_dir', default=None)
@click.option('--group_participants', default=None)
def paired_two(pipeline_dir, z_thresh, p_thresh, conditions,
               condition_type, model_name, output_dir=None,
               group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'paired_two', group_participants,
                                covariate=conditions,
                                condition_type=condition_type,
                                output_dir=output_dir, model_name=model_name)


@load_preset.command()
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('conditions')
@click.argument('condition_type')
@click.argument('model_name')
@click.option('--output_dir', default=None)
@click.option('--group_participants', default=None)
def tripled_two(pipeline_dir, z_thresh, p_thresh, conditions,
                condition_type, model_name, output_dir=None,
                group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        import os
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'tripled_two', group_participants,
                                covariate=conditions,
                                condition_type=condition_type,
                                output_dir=output_dir, model_name=model_name)


# Group analysis - PyBASC
@group.command()
@click.argument('group_config')
def basc(group_config):
    import CPAC.pipeline.cpac_group_runner as cpac_group_runner
    cpac_group_runner.run_basc(group_config)


@group.command(name="mdmr")
@click.argument("group_config", type=click.Path(exists=True))
def group_mdmr(group_config):
    from CPAC.pipeline.cpac_group_runner import run_cwas
    run_cwas(group_config)


@group.command(name="isc")
@click.argument("group_config", type=click.Path(exists=True))
def group_isc(group_config):
    from CPAC.pipeline.cpac_group_runner import run_isc
    run_isc(group_config)

# Utilities
@main.group()
def utils():
    pass


@utils.group()
def data_config():
    pass


@data_config.command()
def new_settings_template():
    from CPAC.utils.build_data_config import util_copy_template
    util_copy_template('data_settings')


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
def group_config():
    pass


@group_config.command()
def new_template():
    from CPAC.utils.build_data_config import util_copy_template
    util_copy_template('group_config')


@utils.group()
def tools():
    pass


@tools.command()
@click.argument('moving_image')
@click.argument('reference')
@click.option('--initial', default=None)
@click.option('--rigid', default=None)
@click.option('--affine', default=None)
@click.option('--nonlinear', default=None)
@click.option('--func_to_anat', default=None)
@click.option('--dim', default=3)
@click.option('--interp', default='Linear')
@click.option('--inverse', default=False)
def ants_apply_warp(moving_image, reference, initial=None, rigid=None,
                    affine=None, nonlinear=None, func_to_anat=None, dim=3,
                    interp='Linear', inverse=False):
    from CPAC.registration.utils import run_ants_apply_warp
    run_ants_apply_warp(moving_image, reference, initial, rigid, affine,
                        nonlinear, func_to_anat, dim, interp, inverse)


@utils.group()
def workflows():
    pass


@workflows.command()
@click.argument('func_ts')
@click.argument('func_brain_mask')
@click.option('--hp', default=0.01)
@click.option('--lp', default=0.1)
def alff(func_ts, func_brain_mask, hp=0.01, lp=0.1):
    from CPAC.alff.alff import run_alff
    paths = run_alff(func_ts, func_brain_mask, hp, lp)
    print(paths)


@utils.group()
def test():
    pass


@test.command()
@click.option('--list', '-l', 'show_list', is_flag=True)
@click.option('--filter', '-f', 'pipeline_filter', default='')
def run_suite(show_list=False, pipeline_filter=''):
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

    if show_list:
        print("")
        print("Availables pipelines:")

    no_params = False
    for config_file in os.listdir(test_config_dir):
        if config_file.startswith('pipe-test_'):
            if pipeline_filter not in config_file:
                continue

            if show_list:
                print("- " + config_file[len('pipe-test_'):])
                continue

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

    if show_list:
        print("")


@test.group()
def functions():
    pass


@functions.command()
@click.argument('pipe_config')
def gather_outputs_func(pipe_config):
    #from CPAC.pipeline.
    #run_gather_outputs_func(pipe_config)
    pass


if __name__ == "__main__":
    main()
