#!/usr/bin/env python

import click

# CLI tree
#
# cpac
#     cpac gui
# cpac individual
#     cpac individual run <pipeline config> <data config>
#     cpac individual quickrun
#         cpac individual quickrun benchmark
#         cpac individual quickrun default <data config>
#         cpac individual quickrun preproc <data config>
# cpac group
#     cpac group feat
#         cpac group feat run <pipeline config>
#         cpac group feat load_preset <args...>
#     cpac group basc
#         cpac group basc run <pipeline config>
#         cpac group basc quickrun <pipeline dir> <working dir> ....
# cpac utils
#     cpac utils data_config
#         cpac utils data_config new_template
#         cpac utils data_config build_data_config <data settings file>
#     cpac utils pipe_config
#         cpac utils pipe_config new_template


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


@feat.command()
def load_preset():
    # TODO: this
    pass


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
