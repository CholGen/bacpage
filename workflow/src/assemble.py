import argparse
import sys
from pathlib import Path

import pandas as pd
import snakemake
import yaml
from snakemake.utils import validate

from workflow.src import common


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Assembles consensus sequence from raw sequencing reads."

    parser.add_argument(
        "directory", type=str, default=".", help="Path to valid project directory [current directory]."
    )
    parser.add_argument(
        "--configfile", type=str, default=".", help="Path to assembly configuration file ['config.yaml']."
    )
    parser.add_argument(
        "--samples", type=str, default=".",
        help="Path to file detailing raw sequencing reads for all samples ['sample_data.csv']."
    )
    parser.add_argument( "--threads", type=int, default=-1, help="Number of threads available for assembly [all]." )
    parser.add_argument( "--verbose", action="store_true", help="Print lots of stuff to screen." )

    parser.set_defaults( command=assemble_entrypoint )


def assemble_entrypoint( args: argparse.Namespace ):
    run_assemble(
        project_directory=args.directory,
        configfile=args.configfile,
        sample_data=args.samples,
        threads=args.threads,
        verbose=args.verbose,
    )


def run_assemble( project_directory: str, configfile: str, sample_data: str, threads: int, verbose: bool = False ):
    # Check project directory
    project_directory = Path( project_directory ).absolute()
    assert project_directory.exists() and project_directory.is_dir(), f"Specified project directory {project_directory} does not exist. Please specify a valid directory."

    # Check config file
    print( "Loading and validating configuration file...", end="" )
    try:
        config = load_configfile( configfile, project_directory )
    except Exception:
        print( "Error" )
        raise
    print( "Done" )

    # Check sample data
    print( "Loading and validating samples metadata...", end="" )
    try:
        metadata, skipped_samples = load_sampledata(
            sample_data, project_directory, check_size=config["preprocessing"]["check_size"],
            minimum_size=config["preprocessing"]["minimum_size"]
        )
    except Exception:
        print( "Error" )
        raise
    print( "Done" )

    if len( skipped_samples ) > 0:
        print( f"Skipping samples [{', '.join( skipped_samples )}] because they have no reads." )

    config["SAMPLES"] = metadata.set_index( "sample" )[["read1", "read2"]].to_dict( orient="index" )

    # Identify reference genes
    # print( "Identifying gene sequences for typing...", end="" )
    # try:
    #    GENES = get_genes( config["reference_genes"] )
    # except Exception:
    #    print( "Error" )
    #    raise
    # print( "Done" )
    # print( f"The following genes will be used: [{', '.join( GENES.keys() )}]\n" )

    # Calculate number of threads if not specified
    useable_threads = common.calculate_threads( threads )

    # Run snakemake command
    snakefile = common.PACKAGE_DIR / "workflow/rules/assemble.smk"
    assert snakefile.exists(), f"Snakefile does not exist. Checking {snakefile}"
    if verbose:
        status = snakemake.snakemake(
            snakefile, printshellcmds=True, forceall=True, force_incomplete=True, workdir=project_directory,
            restart_times=common.RESTART_TIMES,
            config=config, cores=useable_threads, lock=False,
        )
    else:
        status = snakemake.snakemake(
            snakefile, printshellcmds=True, forceall=True, force_incomplete=True, workdir=project_directory,
            restart_times=common.RESTART_TIMES,
            config=config, cores=useable_threads, lock=False, quiet=True
        )
    if not status:
        sys.stderr.write( "Snakemake pipeline did not complete successfully. Check for error messages and rerun." )
        sys.exit( -2 )


def load_configfile( specified_loc: str, project_directory: Path ) -> dict:
    configfile_loc = Path( specified_loc ).absolute()
    if specified_loc == ".":
        configfile_loc = project_directory / common.DEFAULT_CONFIG
        assert configfile_loc.exists(), "Unable to automatically find config in project directory (searching for 'config.yaml'). Please specify a valid configuration file."
    assert configfile_loc.exists(), f"{configfile_loc} does not exist. Please specify a valid file."

    with open( configfile_loc, "r" ) as cf:
        configfile = yaml.safe_load( cf )

    schema_location = common.PACKAGE_DIR / "workflow/schemas/Illumina_config.schema.yaml"
    validate( configfile, schema_location )

    for key in common.CONFIG_PATHS:
        configfile[key] = str( common.normalize_path( configfile[key], common.PACKAGE_DIR / "resources" ) )

    return configfile


def load_sampledata( specified_loc: str, project_directory: Path, check_size: bool = False,
                     minimum_size: int = 100 ) -> (pd.DataFrame, list):
    sampledata_loc = Path( specified_loc )
    if specified_loc == ".":
        sampledata_loc = project_directory / common.DEFAULT_SAMPLEDATA
        if not sampledata_loc.exists():
            sys.stderr.write(
                f"Unable to automatically find sample data file in project directory (searching for '{common.DEFAULT_SAMPLEDATA}'). Please specify a valid sample data file."
            )
            sys.exit( -3 )
    elif not sampledata_loc.is_absolute():
        sampledata_loc = project_directory / sampledata_loc
    if not sampledata_loc.exists():
        sys.stderr.write( f"{sampledata_loc} does not exist. Please specify a valid file." )
        sys.exit( -4 )

    md = pd.read_csv( sampledata_loc )

    schema_location = common.PACKAGE_DIR / "workflow/schemas/Illumina_metadata.schema.yaml"
    validate( md, schema_location )

    md["read1"] = md["read1"].apply( lambda x: str( common.normalize_path( x, project_directory ) ) )
    md["read2"] = md["read2"].apply( lambda x: str( common.normalize_path( x, project_directory ) ) )

    duplicate_samples = md.loc[md["sample"].duplicated(), "sample"]
    if len( duplicate_samples ) > 0:
        sys.stderr.write(
            f"Sample data contains duplicate samples ({', '.join( duplicate_samples.to_list() )}) at rows {duplicate_samples.index.to_list()}"
        )
        sys.exit( -1 )

    skipped_samples = list()
    if check_size:
        md["size"] = md["read1"].apply( lambda x: Path( x ).stat().st_size )
        skipped_samples = md.loc[md["size"] < minimum_size, "sample"].to_list()
        md = md.loc[md["size"] >= minimum_size]

    return md, skipped_samples
