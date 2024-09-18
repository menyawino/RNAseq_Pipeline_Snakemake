#!/usr/bin/env python

import sys
import os
import subprocess
import click
import time
import psutil
import shutil
from datetime import timedelta

# build folders for the pipeline: analysis, benchmarks, results, logs if they don't exist
def build_folders():
    """Build folders for the pipeline if they don't exist."""
    folders = ['analysis', 'benchmarks', 'results', 'logs']
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
       
# Track the resources used by the pipeline     
def track_resources(start_time, net_start, verbose=False):
    """Track and display resource usage."""
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Get CPU and memory usage
    cpu_usage = psutil.cpu_percent(interval=1)
    memory_info = psutil.virtual_memory()
    memory_usage = memory_info.used / (1024 ** 3)  # Convert to GB
    total_memory = memory_info.total / (1024 ** 3)  # Convert to GB

    # Get network usage
    net_end = psutil.net_io_counters()
    bytes_sent = (net_end.bytes_sent - net_start.bytes_sent) / (1024 ** 2)  # Convert to MB
    bytes_recv = (net_end.bytes_recv - net_start.bytes_recv) / (1024 ** 2)  # Convert to MB

    # Get the size of the files in Gb
    output_folder_size = get_folder_size('results') / (1024 ** 3) + get_folder_size('analysis') / (1024 ** 3)


    # Print resource usage stats
    print(f"""
     _____________________________________________________________
    |                                                             |
    |                 Pipeline Resource Usage Report              |
    |_____________________________________________________________|
    |   Runtime:                    {str(timedelta(seconds=elapsed_time))}                |
    |   CPU Usage:                  {cpu_usage}%                          |
    |   Memory Usage:               {memory_usage:.2f} GB / {total_memory:.2f} GB           |
    |   Network Sent:               {bytes_sent:.2f} MB                       |
    |   Network Received:           {bytes_recv:.2f} MB                       |
    |   Output Files Size:         {output_folder_size:.2f} GB                       |
    |_____________________________________________________________|
    """)
    

#  get the size of a folder in bytes 
def get_folder_size(folder):
    """Return the size of a folder in bytes."""
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(folder):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size


# run snakemake with the specified options and configuration
def run_snakemake(configfile, verbose=False, extra_args=[]):
    """Run Snakemake with the specified options and configuration."""

    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')

    # Basic Snakemake command
    cmd = ["snakemake", "-s", snakefile, "--use-conda", "-k", "--benchmark-extended"]

    # Add additional Snakemake arguments
    cmd += list(extra_args)

    if configfile:
        # Only add the specified config file without defaults and system confs
        cmd += ["--configfile", configfile]

    # Print the final command if verbose with cmd list as a string
    if verbose:
        print('Command executed:', ' '.join(cmd))

    start_time = time.time()
    net_start = psutil.net_io_counters()  # Capture network usage at start

    # run snakemake plan to preview the pipeline
    run_snakemake_plan(configfile)
    
    # generate snakemake report for the pipeline
    get_snakemake_report(configfile)
    
    # run Snakemake
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in Snakemake invocation: {e}', file=sys.stderr)
        return e.returncode
    except FileNotFoundError as e:
        print(f'Snakemake not found: {e}', file=sys.stderr)
        return 1 
    
    # display resource usage
    


# run snakemake plan to preview the pipeline
def run_snakemake_plan(configfile):
    """Preview the Snakemake plan before running the pipeline."""
    
    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')
    
    os.system("snakemake -s " + snakefile + " --use-conda --dag \
        --configfile " + configfile + " --quiet \
        | dot -Tpng > results/dag.png")
    os.system("snakemake -s " + snakefile + " --use-conda --rulegraph \
        --configfile " + configfile + " --quiet \
        | dot -Tpng > results/rulegraph.png")


# generate snakemake report for the pipeline
def get_snakemake_report(configfile):
    """Generate a Snakemake report for the pipeline."""
    
    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')
    
    os.system("snakemake -s " + snakefile + " --use-conda --report " + 
              "--configfile " + configfile + " --quiet")
    

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
    
    start_time = time.time()
    net_start = psutil.net_io_counters()
    
    run_snakemake(configfile, verbose=verbose,
                  extra_args=snakemake_args)
    
    track_resources(start_time, net_start, verbose=verbose)
    


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
