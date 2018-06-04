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
            raise click.BadParameter("need to be greater than zero.")


@group.command(name="mdmr")
@click.argument("data_config", type=click.Path(exists=True))
@click.argument("regressors", type=click.Path(exists=True))
@click.option("--permutations", "-p", type=int, default=2000)
@click.option("--column", "-c", type=str, multiple=True)
@click.option("--parallel", callback=validate_mdmr, type=int, default=1)
def group_mdmr(data_config, regressors, permutations, column, parallel):
    print( data_config, permutations)


if __name__ == "__main__":
    main()