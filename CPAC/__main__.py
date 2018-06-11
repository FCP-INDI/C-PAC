#!/usr/bin/env python

import click

# cpac
#     cpac gui
# cpac individual
# cpac group
#     cpac group feat run <pipeline config>
#     cpac group basc run <pipeline config>
# cpac utils
#     cpac utils data_config


@click.group()
def main():
    pass


#@main.command()
#def gui():
#    import CPAC.GUI
#    CPAC.GUI.run()


@main.group()
def utils():
    pass


@utils.command(name="data_config")
def utils_data_config():
    pass
    

@main.group()
def individual():
    pass


@main.group()
@click.pass_context
def participant(ctx):
    ctx.forward(individual)


@main.group()
def group():
    pass


@group.group()
def feat():
    pass


@feat.command()
@click.argument('pipe_config')
@click.argument('pipe_output_dir')
def run(pipe_config, pipe_output_dir):
    import CPAC.pipeline.cpac_group_runner as cpac_group_runner
    cpac_group_runner.run_feat(pipe_config, pipe_output_dir)


@group.group()
def basc():
    pass


@basc.command()
@click.argument('pipe_config')
def run(pipe_config):
    import CPAC.pipeline.cpac_group_runner as cpac_group_runner
    cpac_group_runner.run_basc(pipe_config)


if __name__ == "__main__":
    main()
