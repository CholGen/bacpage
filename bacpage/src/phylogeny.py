import argparse
import sys
from pathlib import Path

import snakemake
from Bio import SeqIO

from bacpage.src import common_funcs

OTHER_IUPAC = {'r', 'y', 's', 'w', 'k', 'm', 'd', 'h', 'b', 'v'}
VALID_CHARACTERS = [{'a'}, {'c'}, {'g'}, {'t'}, {'n'}, OTHER_IUPAC, {'-'}, {'?'}]


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Reconstructs maximum likelihood phylogeny from consensus sequences."

    parser.add_argument(
        "directory", type=str, default=".", help="Path to valid project directory [current directory]."
    )
    parser.add_argument(
        "--configfile", type=str, default=".", help="Path to assembly configuration file ['config.yaml']."
    )
    parser.add_argument(
        "--minimum-completeness", type=float, default=0.9,
        help="minimum coverage required to be included in analysis [0.9]"
    )
    parser.add_argument( "--terra", type=bool, action="store_true",
                         help="Generate input for Terra and print instructions. [False; runs pipeline locally]" )
    parser.add_argument( "--no-detect", type=bool, action="store_true",
                         help="Skip performing recombinant region detection [False; perform detection]" )
    parser.add_argument( "--mask", type=str, default=None,
                         help="gff file used to mask to all sequences prior to tree building. If not specified, sequences will not be masked [False]" )
    parser.add_argument( "--threads", type=int, default=-1, help="Number of threads available for command [all]." )
    parser.add_argument( "--verbose", action="store_true", help="Print lots of stuff to screen." )

    parser.set_defaults( command=phylogeny_entrypoint )


def calculate_completeness( sequence_loc: Path ) -> float:
    record = SeqIO.read( sequence_loc, "fasta" )
    seq = record.seq.lower()
    l = len( seq )
    counts = []

    for v in VALID_CHARACTERS:
        counts.append( sum( map( lambda x: seq.count( x ), v ) ) )
    invalid_nucleotides = l - sum( counts )

    if invalid_nucleotides > 0:
        print( "Invalid characters in sequence. Might not be a valid nucleotide sequence." )

    return 1.0 - (counts[4] / l)


def load_input( directory: str, minimum_completeness: float ) -> dict[str, Path]:
    search_directory = Path( directory ).absolute()
    if common_funcs.is_project_directory( search_directory ):
        search_directory = search_directory / "results/consensus"

    fastas = common_funcs.find_files( directory=search_directory, extensions=[".fa", ".fasta"] )

    if minimum_completeness > 0:
        for name in list( fastas.keys() ):
            completeness = calculate_completeness( fastas[name] )
            if completeness < minimum_completeness:
                del fastas[name]
    return fastas


def validate_sequences( sequence_paths: dict[str, Path], background_dataset: Path = None ):
    seen_names = list()

    if background_dataset and (background_dataset != ""):
        for record in SeqIO.parse( background_dataset, "fasta" ):
            seen_names.append( record.name )

    for path in sequence_paths.values():
        for record in SeqIO.parse( path, "fasta" ):
            seen_names.append( record.name )

    if len( seen_names ) != len( set( seen_names ) ):
        sys.stderr.write( "Sequence names are not unique. Remove duplicate sequences.\n" )
        sys.exit( -6 )


def postamble( directory: Path ):
    print()
    print( f"Reconstructed phylogeny is available at {directory / 'results/phylogeny/phylogeny.tree'}" )
    print(
        "Open the file in any tree viewer to visualize. We recommend figtree (http://tree.bio.ed.ac.uk/software/figtree/)." )


def reconstruct_phylogeny( project_directory: str, configfile: str, minimum_completeness: float = 0,
                           terra: bool = False, mask: str = None, skip_detect: bool = False, threads: int = 1,
                           verbose: bool = True,
                           dryrun: bool = False ):
    project_path = Path( project_directory ).absolute()

    if not project_path.exists() or not project_path.is_dir():
        sys.stderr.write(
            f"Specified project directory {project_path} does not exist. Please specify a valid directory.\n" )
        sys.exit( -1 )

    input_sequences = load_input( project_directory, minimum_completeness=minimum_completeness )

    # Check config file
    print( "Loading and validating configuration file...", end="" )
    try:
        config = common_funcs.load_configfile( configfile, project_path )
    except Exception:
        print( "Error" )
        raise
    print( "Done" )

    if config["background_dataset"] not in ["", "<background-dataset-path>"]:
        background_dataset = Path( config["background_dataset"] ).absolute()
        config["BACKGROUND"] = background_dataset
    else:
        config["BACKGROUND"] = ""

    validate_sequences( input_sequences, config["BACKGROUND"] )

    config["SAMPLES"] = input_sequences

    # Add command line arguments
    config["TERRA"] = terra
    config["DETECT"] = not skip_detect
    config["MASK"] = False
    config["MASK_LOCATION"] = ""
    if mask:
        mask_loc = Path( mask ).absolute()
        if not mask_loc.exists():
            sys.stderr.write( f"{mask_loc} does not exist. Please specify a valid GFF file." )
        config["MASK"] = True
        config["MASK_LOCATION"] = str( mask_loc )

    # Calculate number of threads
    useable_threads = common_funcs.calculate_threads( threads )

    # load snakefile
    snakefile = common_funcs.PACKAGE_DIR / "rules/phylogeny.smk"
    assert snakefile.exists(), f"Snakefile does not exist. Searching for {snakefile}"

    if dryrun:
        return config, snakefile

    status = snakemake.snakemake(
        snakefile, printshellcmds=True, forceall=True, force_incomplete=True, workdir=project_directory,
        restart_times=common_funcs.RESTART_TIMES, config=config, cores=useable_threads, lock=False, quiet=not verbose
    )
    if not status:
        sys.stderr.write( "Snakemake pipeline did not complete successfully. Check for error messages and rerun.\n" )
        sys.exit( -2 )

    postamble( project_path )


def phylogeny_entrypoint( args: argparse.Namespace ):
    reconstruct_phylogeny(
        project_directory=args.directory,
        configfile=args.configfile,
        minimum_completeness=args.mininum_completeness,
        terra=args.terra,
        mask=args.mask,
        skip_detect=args.no_detect,
        threads=args.threads,
        verbose=args.verbose,
    )
