#!/usr/bin/env python

import sys
import os
import subprocess
import click

def run_snakemake(configfile, verbose=False, extra_args=[]):
    """Run Snakemake with the specified options and configuration."""

    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')

    # Basic Snakemake command
    cmd = ["snakemake", "-s", snakefile, "--use-conda"]

    # Add additional Snakemake arguments
    cmd += list(extra_args)

    if configfile:
        # Only add the specified config file without defaults and system confs
        cmd += ["--configfile", configfile]

    # Print the final command if verbose with cmd list as a string
    if verbose:
        print('Command executed:', ' '.join(cmd))

    # Run Snakemake
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in Snakemake invocation: {e}', file=sys.stderr)
        return e.returncode
    except FileNotFoundError as e:
        print(f'Snakemake not found: {e}', file=sys.stderr)
        return 1 


def run_snakemake_plan():
    """Preview the Snakemake plan before running the pipeline."""
    os.system("snakemake -s /mnt/d/omar_wsl/rnaseq_pipeline_snakemake/workflow/Snakefile \
            --use-conda --dag --configfile workflow/config_new_tuxedo.yml --quiet \
            | dot -Tpng > results/dag.png")
    os.system("snakemake -s /mnt/d/omar_wsl/rnaseq_pipeline_snakemake/workflow/Snakefile \
            --use-conda --rulegraph --configfile workflow/config_new_tuxedo.yml --quiet \
            | dot -Tpng > results/rulegraph.png")


@click.group()
def cli():
    """Define the CLI group."""
    pass

@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('--verbose', is_flag=True, help="Enable verbose output.")
@click.argument('snakemake_args', nargs=-1)
def run(configfile, snakemake_args, verbose):
    """Execute workflow (using Snakemake underneath)."""
    run_snakemake_plan()
    run_snakemake(configfile, verbose=verbose,
                  extra_args=snakemake_args)



# @click.command()
# @click.argument('configfile')
# def check(configfile):
#     """Check configuration."""
#     run_snakemake(configfile, extra_args=['check'])

# @click.command()
# @click.argument('configfile')
# def showconf(configfile):
#     """Show full configuration across project config files."""
#     run_snakemake(configfile, extra_args=['showconf'])

# @click.command()
# def info():
#     """Provide basic install/config file info."""
#     try:
#         from version import version 
#     except ImportError:
#         version = "unknown"
#     print(f"""
# This is your package version v{version}

# Package install path: {os.path.dirname(__file__)}
# Snakemake Snakefile: 'workflow/Snakefile'
# """)


# Register commands with the Click CLI
cli.add_command(run)
# cli.add_command(check)
# cli.add_command(showconf)
# cli.add_command(info)


# ANSI color codes
GRE = '\033[92m'  # Green color
NC = '\033[0m'    # No color


def main():
    """Main entry point."""
    print(f"""
 _____________________________________________________________
|                                                             |
|         ██████  █████  ██████  ██████  ██  ██████           |
|        ██      ██   ██ ██   ██ ██   ██ ██ ██    ██          |
|        ██      ███████ ██████  ██   ██ ██ ██    ██          |
|        ██      ██   ██ ██   ██ ██   ██ ██ ██    ██          |
|         ██████ ██   ██ ██   ██ ██████  ██  ██████           |
|                                                             |
|     ██ ███    ██ ██████  ██████  ██████  ███      ███       |
|     ██ ████   ██ ██     ██    ██ ██   ██ ████    ████       |
|     ██ ██ ██  ██ ██████ ██    ██ ██████  ██ ██  ██ ██       |
|     ██ ██  ██ ██ ██     ██    ██ ██   ██ ██  ████  ██       |
|     ██ ██   ████ ██      ██████  ██   ██ ██   ██   ██       |
|                                                             |
|   ██████  ███    ██  █████       ███████ ███████  ██████    |
|   ██   ██ ████   ██ ██   ██      ██      ██      ██    ██   |
|   ██████  ██ ██  ██ ███████ ████ ███████ ███████ ██    ██   |
|   ██   ██ ██  ██ ██ ██   ██           ██ ██      ██    ██   |
|   ██   ██ ██   ████ ██   ██      ███████ ███████  ██████▄   |
|                                                             |
| {GRE}      RNAseq Analysis Toolkit for Cardiology Research{NC}       |
|_____________________________________________________________|

""")
    cli()

if __name__ == '__main__':
    main()
