import argparse
import sys
from pathlib import Path

import pandas as pd
import snakemake
from snakemake.utils import validate

from bacpage.src import common_funcs, identify


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Assembles consensus sequence from raw sequencing reads."

    parser.add_argument(
        "directory", type=str, nargs="?", default=".", help="Path to valid project directory [current directory]."
    )
    parser.add_argument(
        "--denovo", action="store_true", help="Perform de novo assembly rather than reference-based assembly."
    )
    parser.add_argument(
        "--configfile", type=str, default=".", help="Path to assembly configuration file ['config.yaml']."
    )
    parser.add_argument(
        "--samples", type=str, default=".",
        help="Path to file detailing raw sequencing reads for all samples ['sample_data.csv']."
    )
    parser.add_argument(
        "--delim", default="_", type=str, help="deliminator to extract sample name from file name [_]"
    )
    parser.add_argument(
        "--index", default=0, type=int, help="index of sample name after splitting file name by delim [0]"
    )
    parser.add_argument(
        "--no-qa", action="store_true", help="Whether to skip quality assessment of assemblies [False]"
    )
    parser.add_argument( "--threads", type=int, default=-1, help="Number of threads available for assembly [all]." )
    parser.add_argument( "--verbose", action="store_true", help="Print lots of stuff to screen." )

    parser.set_defaults( command=assemble_entrypoint )


def postamble( denovo: bool, directory: Path ):
    if denovo:
        print()
        print( "Successfully performed de novo assembly of your samples" )
        print( f"Annotated assemblies are available at {directory / 'results/assembly/'}." )
        print(
            f"Quality metrics of your input data and assemblies are available at {directory / 'reports/qc_report.html'}. Open this file in a web browser to view." )
        print()
        print( "To perform antimicrobial resistance profiling, run `bacpage profile`." )
    else:
        print()
        print( "Successfully performed reference-based assembly of your samples." )
        print( f"Consensus sequences are available at {directory / 'results/consensus'}." )
        print(
            f"Quality metrics of your input data and consensus sequences are available at {directory / 'reports/qc_report.html'}. Open this file in a web browser to view." )
        print()
        print( "Generate a phylogenetic tree incorporating these samples using `bacpage phylogeny`. " )
        print( "Or, determine the presense of antimicrobial resistance genes using `bacpage profile`." )


def assemble_entrypoint( args: argparse.Namespace ):
    run_assemble(
        project_directory=args.directory,
        configfile=args.configfile,
        sample_data=args.samples,
        denovo=args.denovo,
        qc=not args.no_qa,
        delim=args.delim,
        index=args.index,
        threads=args.threads,
        verbose=args.verbose,
    )


def run_assemble( project_directory: str, configfile: str, sample_data: str, denovo: bool, qc: bool, threads: int,
                  verbose: bool = False, delim: str = "_", index: int = 0 ):
    # Check project directory
    project_directory = Path( project_directory ).absolute()
    assert project_directory.exists() and project_directory.is_dir(), f"Specified project directory {project_directory} does not exist. Please specify a valid directory."

    # Check config file
    print( "Loading and validating configuration file...", end="" )
    try:
        config = common_funcs.load_configfile( configfile, project_directory )
    except Exception:
        print( "Error" )
        raise
    print( "Done" )

    # Perform QC if enabled
    config["QC"] = qc

    # Check sample data
    print( "Loading and validating samples data...", end="" )
    try:
        metadata, skipped_samples = load_sampledata(
            specified_loc=sample_data,
            project_directory=project_directory,
            check_size=config["preprocessing"]["check_size"],
            minimum_size=config["preprocessing"]["minimum_size"],
            delim=delim,
            index=index
        )
    except Exception:
        print( "Error" )
        raise
    print( f"Found {len( metadata )} samples. " )

    if len( skipped_samples ) > 0:
        print( f"Skipping samples [{', '.join( skipped_samples )}] because they have no reads." )

    # Add samples dictionary to config, so snakemake can utilize it.
    config["SAMPLES"] = metadata.set_index( "sample" )[["read1", "read2"]].to_dict( orient="index" )

    # Calculate number of threads if not specified
    useable_threads = common_funcs.calculate_threads( threads )

    # Run the appropriate snakemake command
    config["DENOVO"] = denovo
    if denovo:
        snakefile = common_funcs.PACKAGE_DIR / "rules/denovo.smk"
    else:
        snakefile = common_funcs.PACKAGE_DIR / "rules/assemble.smk"
    assert snakefile.exists(), f"Snakefile does not exist. Checking {snakefile}"
    if verbose:
        status = snakemake.snakemake(
            snakefile, force_incomplete=True, workdir=project_directory,
            restart_times=common_funcs.RESTART_TIMES,
            config=config, cores=useable_threads, lock=False,
            quiet=False, printshellcmds=True, printreason=False
        )
    else:
        status = snakemake.snakemake(
            snakefile, force_incomplete=True, workdir=project_directory,
            restart_times=common_funcs.RESTART_TIMES,
            config=config, cores=useable_threads, lock=False,
            quiet=True, printshellcmds=False, printreason=False
        )
    if not status:
        sys.stderr.write( "Snakemake pipeline did not complete successfully. Check for error messages and rerun.\n" )
        sys.exit( -2 )

    postamble( denovo, project_directory )


def find_sampledata( specified_loc: str, project_directory: Path, delim: str = "_", index: int = 0 ) -> Path:
    if specified_loc != ".":
        specified_path = Path( specified_loc )
        if not specified_path.exists():
            sys.stderr.write(
                f"Cannot parse sample data because {specified_loc} does not exists. Please specify a valid path." )
            sys.exit( -1 )
        return specified_path

    detected_path = project_directory / common_funcs.DEFAULT_SAMPLEDATA
    if detected_path.exists():
        return detected_path

    search_directory = project_directory / "input"
    if not search_directory.exists():
        sys.stderr.write(
            f"Cannot automatically find samples because {common_funcs.DEFAULT_SAMPLEDATA} and {project_directory / 'input'} do not exist. Please specify a sample data file or add input files to `input/`." )
        sys.exit( -1 )

    samples = identify.generate_sample_data( directory=search_directory, delim=delim, index=index )
    samples_cache = detected_path
    identify.write_samples_to_file( sample_data=samples, output=samples_cache )
    return samples_cache


def load_sampledata( specified_loc: str, project_directory: Path, check_size: bool = False,
                     minimum_size: int = 100, delim: str = "_", index: int = 0 ) -> (pd.DataFrame, list):
    """ Attempts for find sample data using user supplied information. If sample data file is directly specified, use it, else
    search for the file in the project directory.

    Parameters
    ----------
    specified_loc : str
        Location of sample data file. Must be in CSV format. Use "." to attempt to automatically find file in project_direcctory.
    project_directory : pathlib.Path
        Location of project directory.
    check_size : bool
        Indicates whether to validate input file sizes.
    minimum_size : int
        Minimum file size for an input file to be considered valid. Not considered if check_size is False.
    delim : str
        When automatically searching for samples from sequencing files, use this delim to split the file name. Default is "_".
    index : int
        When automatically searching for samples from seuqencing files, take this field from the split file name as the sample name. Default is 0.

    Returns
    -------
    pd.DataFrame
        Dataframe matching sample names to input raw sequencing data.
    list
        Sample names which were removed from sample data because their input files are smaller than minimum_size.
    """
    sampledata_loc = find_sampledata(
        specified_loc=specified_loc,
        project_directory=project_directory,
        delim=delim,
        index=index
    )

    md = pd.read_csv( sampledata_loc )

    schema_location = common_funcs.PACKAGE_DIR / "schemas/Illumina_metadata.schema.yaml"
    validate( md, schema_location )

    md["read1"] = md["read1"].apply( lambda x: str( common_funcs.normalize_path( x, project_directory ) ) )
    md["read2"] = md["read2"].apply( lambda x: str( common_funcs.normalize_path( x, project_directory ) ) )

    duplicate_samples = md.loc[md["sample"].duplicated(), "sample"]
    if len( duplicate_samples ) > 0:
        sys.stderr.write(
            f"Sample data contains duplicate samples ({', '.join( duplicate_samples.to_list() )}) at rows {duplicate_samples.index.to_list()}\n"
        )
        sys.exit( -1 )

    skipped_samples = list()
    if check_size:
        md["size"] = md["read1"].apply( lambda x: Path( x ).stat().st_size )
        skipped_samples = md.loc[md["size"] < minimum_size, "sample"].to_list()
        md = md.loc[md["size"] >= minimum_size]

    return md, skipped_samples
