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


@utils.command(name="data_config")
def utils_data_config():
    CPAC.GUI.run()


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


if __name__ == "__main__":
    main()