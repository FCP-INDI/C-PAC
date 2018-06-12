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


def validate_mdmr(ctx, param, value):
    if param.name == "parallel":
        if value < 1:
            raise click.BadParameter("Need to be greater than zero.")
        return value

    if param.name == "permutations":
        if value < 1:
            raise click.BadParameter("Need to be greater than zero.")
        return value


@group.command(name="mdmr")
@click.argument("pipeline_config", type=click.Path(exists=True))
def group_mdmr(pipeline_config):
    from CPAC.pipeline.cpac_group_runner import run_mdmr
    
    run_mdmr(pipeline_config)


if __name__ == "__main__":
    main()