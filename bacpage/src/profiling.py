import argparse
import sys
from pathlib import Path

import snakemake

from bacpage.src import common_funcs


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Reconstructs maximum likelihood phylogeny from consensus sequences."

    parser.add_argument(
        "directory", type=str, nargs="?", default=".", help="Path to valid project directory [current directory]."
    )
    parser.add_argument(
        "--configfile", type=str, default=".", help="Path to assembly configuration file ['config.yaml']."
    )
    parser.add_argument(
        "--database", type=str, default=None, help="Database to use for antimicrobial resistance profiling."
    )
    parser.add_argument( "--threads", type=int, default=-1, help="Number of threads available for command [all]." )
    parser.add_argument( "--verbose", action="store_true", help="Print lots of stuff to screen." )

    parser.set_defaults( command=profiling_entrypoint )


def find_sequences( directory: Path ) -> dict[str, Path]:
    assert directory.exists(), f"{directory} does not exist."

    # search for assemblies
    assembly_loc = directory / "results/assembly"
    consensus_loc = directory / "results/consensus"
    extensions = [".fa", ".fasta"]
    for search_location, exts in [(assembly_loc, [".gff"]), (consensus_loc, extensions), (directory, extensions)]:
        if search_location.exists():
            files = common_funcs.find_files( search_location, exts )
            if len( files ) > 0:
                return files
    sys.stderr.write(
        f"Unable to find assemblies, consensus sequences, or loose fastas in {directory}. Please indicate a valid project directory."
    )
    sys.exit( -8 )


def profile_sequences( project_directory, configfile, database, threads, verbose ):
    project_path = Path( project_directory ).absolute()
    assert project_path.exists() and project_path.is_dir(), f"Specified project directory {project_path} does not exist. Please specify a valid directory."

    print( "Searching for sequences...", end="" )
    try:
        input_sequences = find_sequences( project_path )
    except:
        print( "Error" )
        raise
    print( f"Done. {len( input_sequences )} found." )

    # Check config file
    print( "Loading and validating configuration file...", end="" )
    try:
        config = common_funcs.load_configfile( configfile, project_path )
    except Exception:
        print( "Error" )
        raise
    print( "Done" )

    config["SAMPLES"] = input_sequences

    if database:
        config["antibiotic_resistance"]["database"] = database

    # Calculate number of threads
    useable_threads = common_funcs.calculate_threads( threads )

    snakefile = common_funcs.PACKAGE_DIR / "rules/profiling.smk"
    assert snakefile.exists(), f"Snakefile {snakefile} does not exist."
    status = snakemake.snakemake(
        snakefile, printshellcmds=True, forceall=True, force_incomplete=True, workdir=project_directory,
        restart_times=common_funcs.RESTART_TIMES, config=config, cores=useable_threads, lock=False, quiet=not verbose
    )
    if not status:
        sys.stderr.write( "Snakemake pipeline did not complete successfully. Check for error messages and rerun." )
        sys.exit( -2 )


def profiling_entrypoint( args: argparse.Namespace ):
    profile_sequences(
        project_directory=args.directory,
        configfile=args.configfile,
        database=args.database,
        threads=args.threads,
        verbose=args.verbose
    )

# Identify reference genes
# print( "Identifying gene sequences for typing...", end="" )
# try:
#    GENES = get_genes( config["reference_genes"] )
# except Exception:
#    print( "Error" )
#    raise
# print( "Done" )
# print( f"The following genes will be used: [{', '.join( GENES.keys() )}]\n" )
