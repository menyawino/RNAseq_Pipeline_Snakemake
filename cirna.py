#!/usr/bin/env python

import sys
import os
import subprocess
import click


# run snakemake with the specified options and configuration
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


# run snakemake plan to preview the pipeline
def run_snakemake_plan():
    """Preview the Snakemake plan before running the pipeline."""
    
    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')
    
    os.system("snakemake -s " + snakefile + " --use-conda --dag \
        --configfile workflow/config_new_tuxedo.yml --quiet \
        | dot -Tpng > results/dag.png")
    os.system("snakemake -s " + snakefile + " --use-conda --rulegraph \
        --configfile workflow/config_new_tuxedo.yml --quiet \
        | dot -Tpng > results/rulegraph.png")


# build folders for the pipeline: analysis, benchmarks, results, logs if they don't exist
def build_folders():
    """Build folders for the pipeline if they don't exist."""
    folders = ['analysis', 'benchmarks', 'results', 'logs']
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
    

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
    build_folders()
    run_snakemake_plan()
    run_snakemake(configfile, verbose=verbose,
                  extra_args=snakemake_args)


cli.add_command(run)


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
