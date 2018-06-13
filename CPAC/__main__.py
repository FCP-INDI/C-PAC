#!/usr/bin/env python

import CPAC
import click


@click.group()
def main():
    pass


@main.command()
def gui():
    CPAC.GUI.run()


@main.group()
def utils():
    pass


@utils.command(name="data_config_setup")
def utils_data_config_setup():
    CPAC.GUI.run()


@main.group()
def individual():
    pass


@main.group()
def participant():
    ctx.forward(individual)


@main.group()
def group():
    pass


@group.command(name="cwas")
@click.argument("pipeline_config", type=click.Path(exists=True))
def group_cwas(pipeline_config):
    from CPAC.pipeline.cpac_group_runner import run_cwas
    
    run_cwas(pipeline_config)


if __name__ == "__main__":
    main()