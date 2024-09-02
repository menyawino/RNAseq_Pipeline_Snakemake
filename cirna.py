#!/usr/bin/env python

import sys
import os
import subprocess
import click

def get_snakefile_path(name):
    """Get the path of the Snakefile."""
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, name)
    return snakefile

def run_snakemake(configfile, no_use_conda=False, verbose=False,
                  snakefile_name='workflow/Snakefile', extra_args=[]):
    """Run Snakemake with the specified options and configuration."""
    # Find the Snakefile relative to the package path
    snakefile = get_snakefile_path(snakefile_name)

    # Basic Snakemake command
    cmd = ["snakemake", "-s", snakefile]

    # Add --use-conda if not disabled
    if not no_use_conda:
        cmd += ["--use-conda"]

    # Set default -j to 1; can be overridden later
    cmd += ["-j", "1"]

    # Add additional Snakemake arguments
    cmd += list(extra_args)

    if configfile:
        # Only add the specified config file without defaults and system confs
        cmd += ["--configfile", configfile]

    if verbose:
        print('Final command:', cmd)

    # Run Snakemake
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in Snakemake invocation: {e}', file=sys.stderr)
        return e.returncode

@click.group()
def cli():
    """Define the CLI group."""
    pass

@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('--no-use-conda', is_flag=True, default=False)
@click.option('--verbose', is_flag=True)
@click.argument('snakemake_args', nargs=-1)
def run(configfile, snakemake_args, no_use_conda, verbose):
    """Execute workflow (using Snakemake underneath)."""
    run_snakemake(configfile, snakefile_name='workflow/Snakefile',
                  no_use_conda=no_use_conda, verbose=verbose,
                  extra_args=snakemake_args)

@click.command()
@click.argument('configfile')
def check(configfile):
    """Check configuration."""
    run_snakemake(configfile, extra_args=['check'])

@click.command()
@click.argument('configfile')
def showconf(configfile):
    """Show full configuration across project config files."""
    run_snakemake(configfile, extra_args=['showconf'])

@click.command()
def info():
    """Provide basic install/config file info."""
    from .version import version
    print(f"""
This is your package version v{version}

Package install path: {os.path.dirname(__file__)}
Snakemake Snakefile: {get_snakefile_path('Snakefile')}
""")

cli.add_command(run)
cli.add_command(check)
cli.add_command(showconf)
cli.add_command(info)

def main():
    """Main entry point."""
    cli()

if __name__ == '__main__':
    main()
