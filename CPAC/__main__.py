#!/usr/bin/env python
# Copyright (C) 2018-2022  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
import os
import pkg_resources as p
import click
from click_aliases import ClickAliasedGroup
from CPAC.utils.docs import version_report

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
    """Display environment version information"""
    import CPAC
    print('\n'.join(['Environment', '===========', version_report(),
                     f'C-PAC version: {CPAC.__version__}']))


@main.command()
@click.argument('data_config')
@click.option('--pipe-config', '--pipe_config')
@click.option('--num-cores', '--num_cores')
@click.option('--ndmg-mode', '--ndmg_mode', is_flag=True)
@click.option('--debug', is_flag=True)
def run(data_config, pipe_config=None, num_cores=None, ndmg_mode=False,
        debug=False):
    if not pipe_config:
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_template.yml"))

    if pipe_config == 'benchmark-ants':
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_benchmark-ANTS.yml"))

    if pipe_config == 'benchmark-fnirt':
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_benchmark-FNIRT.yml"))

    if pipe_config == 'anat-only':
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_anat-only.yml"))

    if data_config == 'benchmark-data':
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_cpac_benchmark.yml"))

    if data_config == 'ADHD200':
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-ADHD200.yml"))
    if data_config == 'ADHD200_2':
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-ADHD200_only2.yml"))
    if data_config == 'ABIDE':
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-ABIDE.yml"))
    if data_config == 'NKI-RS':
        data_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_config_S3-BIDS-NKI-RocklandSample.yml"))

    if ndmg_mode:
        pipe_config = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_ndmg.yml"))

    import CPAC.pipeline.cpac_runner as cpac_runner
    cpac_runner.run(data_config, pipe_config, num_subs_at_once=num_cores,
                    debug=debug)


# Group analysis
@main.group()
def group():
    pass


# Group analysis - FSL FEAT
@group.group(cls=ClickAliasedGroup)
def feat():
    pass


@feat.command()
@click.argument('group_config')
def build_models(group_config):
    import CPAC.pipeline.cpac_group_runner as cgr
    cgr.build_feat_models(group_config)


@feat.command(name='run')
@click.argument('group_config')
def run_feat(group_config):
    import CPAC.pipeline.cpac_group_runner as cgr
    cgr.run_feat(group_config)


@feat.command()
@click.argument('group_config')
def randomise(group_config):
    import CPAC.pipeline.cpac_group_runner as cgr
    cgr.run_feat(group_config, feat=False)


@feat.group(aliases=['load_preset'], cls=ClickAliasedGroup)
def load_preset():
    pass

@load_preset.command(aliases=['single_grp_avg'])
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('model_name')
@click.option('--output-dir', '--output_dir', default=None)
@click.option('--group-participants', '--group_participants', default=None)
def single_grp_avg(pipeline_dir, z_thresh, p_thresh, model_name,
                   output_dir=None, group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'single_grp', group_participants,
                                output_dir=output_dir,
                                model_name=model_name)

@load_preset.command(aliases=['single_grp_cov'])
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('pheno_file')
@click.argument('pheno_sub')
@click.argument('covariate')
@click.argument('model_name')
@click.option('--output-dir', '--output_dir', default=None)
@click.option('--group-participants', '--group_participants', default=None)
def single_grp_cov(pipeline_dir, z_thresh, p_thresh, pheno_file,
                   pheno_sub, covariate, model_name, output_dir=None,
                   group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'single_grp_cov', group_participants,
                                pheno_file=pheno_file,
                                pheno_sub_label=pheno_sub,
                                covariate=covariate, output_dir=output_dir,
                                model_name=model_name)

@load_preset.command(aliases=['unpaired_two'])
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('pheno_file')
@click.argument('pheno_sub')
@click.argument('covariate')
@click.argument('model_name')
@click.option('--output-dir', '--output_dir', default=None)
@click.option('--group-participants', '--group_participants', default=None)
def unpaired_two(pipeline_dir, z_thresh, p_thresh, pheno_file,
                 pheno_sub, covariate, model_name, output_dir=None,
                 group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'unpaired_two', group_participants,
                                pheno_file=pheno_file,
                                pheno_sub_label=pheno_sub,
                                covariate=covariate, output_dir=output_dir,
                                model_name=model_name)

@load_preset.command(aliases=['paired_two'])
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('conditions')
@click.argument('condition_type')
@click.argument('model_name')
@click.option('--output-dir', '--output_dir', default=None)
@click.option('--group-participants', '--group_participants', default=None)
def paired_two(pipeline_dir, z_thresh, p_thresh, conditions,
               condition_type, model_name, output_dir=None,
               group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), 'cpac_group_analysis')
    create_fsl_flame_preset.run(pipeline_dir, 'all', z_thresh, p_thresh,
                                'paired_two', group_participants,
                                covariate=conditions,
                                condition_type=condition_type,
                                output_dir=output_dir, model_name=model_name)


@load_preset.command(aliases=['tripled_two'])
@click.argument('pipeline_dir')
@click.argument('z_thresh')
@click.argument('p_thresh')
@click.argument('conditions')
@click.argument('condition_type')
@click.argument('model_name')
@click.option('--output-dir', '--output_dir', default=None)
@click.option('--group-participants', '--group_participants', default=None)
def tripled_two(pipeline_dir, z_thresh, p_thresh, conditions,
                condition_type, model_name, output_dir=None,
                group_participants=None):
    from CPAC.utils import create_fsl_flame_preset
    if not output_dir:
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

# Group analysis - QPP
@group.command()
@click.argument('group_config')
def qpp(group_config):
    from CPAC.pipeline.cpac_group_runner import run_qpp
    run_qpp(group_config)


# Utilities
@main.group(cls=ClickAliasedGroup)
def utils():
    pass


@utils.command()
@click.argument('crash_file')
def crash(crash_file):

    import mock

    def accept_all(object, name, value):
        return value

    with mock.patch('nipype.interfaces.base.traits_extension.File.validate', side_effect=accept_all):
        from nipype.scripts.crash_files import display_crash_file
        display_crash_file(crash_file, False, False, None)


@utils.group(aliases=['data_config'], cls=ClickAliasedGroup)
@click.option('--tracking-opt-out', '--tracking_opt-out', is_flag=True,
              help='Disable usage tracking.')
def data_config(tracking_opt_out):
    if not tracking_opt_out:
        # pylint: disable=import-outside-toplevel
        from CPAC.utils.ga import track_config
        track_config('cli')


@data_config.command(aliases=['new_settings_template'])
def new_settings_template():
    from CPAC.utils.build_data_config import util_copy_template
    util_copy_template('data_settings')


@data_config.command()
@click.argument('data_settings_file')
def build(data_settings_file):
    from CPAC.utils.build_data_config import run
    run(data_settings_file)


@utils.group(aliases=['pipe_config'], cls=ClickAliasedGroup)
@click.option('--tracking-opt-out', '--tracking_opt-out', is_flag=True,
              help='Disable usage tracking.')
def pipe_config(tracking_opt_out):
    if not tracking_opt_out:
        # pylint: disable=import-outside-toplevel
        from CPAC.utils.ga import track_config
        track_config('cli')


@pipe_config.command(name='new-template', aliases=['new_template'])
def new_pipeline_template():
    from CPAC.utils.build_data_config import util_copy_template
    util_copy_template('pipeline_config')


@utils.group(aliases=['group_config'], cls=ClickAliasedGroup)
@click.option('--tracking-opt-out', '--tracking_opt-out', is_flag=True,
              help='Disable usage tracking.')
def group_config(tracking_opt_out):
    if not tracking_opt_out:
        # pylint: disable=import-outside-toplevel
        from CPAC.utils.ga import track_config
        track_config('cli')


@group_config.command(name='new-template', aliases=['new_template'])
def new_group_template():
    from CPAC.utils.build_data_config import util_copy_template
    util_copy_template('group_config')


@utils.group(cls=ClickAliasedGroup)
def tools():
    pass


@tools.command(aliases=['ants_apply_warp'])
@click.argument('moving_image')
@click.argument('reference')
@click.option('--initial', default=None)
@click.option('--rigid', default=None)
@click.option('--affine', default=None)
@click.option('--nonlinear', default=None)
@click.option('--func-to-anat', '--func_to_anat', default=None)
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


@utils.command()
@click.argument('directory')
def repickle(directory):
    from CPAC.utils import repickle as repickle_util
    if os.path.exists(directory):
        repickle_util(directory)


@utils.group(cls=ClickAliasedGroup)
def test():
    pass


@test.command(aliases=['run_suite'])
@click.option('--list', '-l', 'show_list', is_flag=True)
@click.option('--filter', '-f', 'pipeline_filter', default='')
def run_suite(show_list=False, pipeline_filter=''):
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
            cpac_runner.run(data, pipe)

    if show_list:
        print("")


@test.group(cls=ClickAliasedGroup)
def functions():
    pass


@functions.command(aliases=['gather_outputs_func'])
@click.argument('pipe_config')
def gather_outputs_func(pipe_config):
    #from CPAC.pipeline.
    #run_gather_outputs_func(pipe_config)
    pass


if __name__ == "__main__":
    main()
