import sys
from importlib.resources import files
from pathlib import Path

import yaml
from Bio import SeqIO
from snakemake import WorkflowError
from snakemake.utils import validate

DEFAULT_CONFIG = "config.yaml"
DEFAULT_SAMPLEDATA = "sample_data.csv"
PACKAGE_DIR = files( "bacpage" )
CONFIG_PATHS = {"Illumina": ["reference"], "phylogeny": ["reference", "recombinant_mask"], "assemble": ["reference"]}
RESTART_TIMES = 0
OTHER_IUPAC = {'r', 'y', 's', 'w', 'k', 'm', 'd', 'h', 'b', 'v'}
VALID_CHARACTERS = [{'a'}, {'c'}, {'g'}, {'t'}, {'n'}, OTHER_IUPAC, {'-'}, {'?'}]


def find_files( directory: Path, extensions: list[str] ) -> dict[str, Path]:
    search_path = directory.absolute()
    print( search_path )
    found = dict()
    for file in search_path.iterdir():
        if file.suffix in extensions:
            name = file.stem
            if name in found:
                sys.stderr.write( f"{name} was found in two or more files. Duplicate sequences cannot be processed.\n" )
                sys.exit( -5 )
            path = search_path / file
            found[name] = path
    return found


def is_project_directory( directory ):
    return (directory / "sample_data.csv").exists() and (directory / "config.yaml").exists()


def normalize_path( value: str, working_directory: Path ) -> Path:
    path = Path( value )
    if not path.is_absolute():
        return working_directory / path
    else:
        return path


def load_configfile( specified_loc: str, project_directory: Path, schema: str = "Illumina" ) -> dict:
    """ Attempts for find config file using user supplied information. If config file is directly specified, use it, else
    search for the config file in the project directory.

    Parameters
    ----------
    specified_loc: str
        Path to config file. Pass "." to automatically search for file in project directory.
    project_directory: pathlib.Path
        Path to project directory. Used if config file path is not specified and to normalize relative paths in the config file.
    schema: str
        Specify a different schema to use to validate configuration file.

    Returns
    -------
    dict
        Config file loaded as a python object.
    """
    configfile_loc = Path( specified_loc ).absolute()
    if specified_loc == ".":
        configfile_loc = project_directory / DEFAULT_CONFIG
        assert configfile_loc.exists(), "Unable to automatically find config in project directory (searching for 'config.yaml'). Please specify a valid configuration file."
    assert configfile_loc.exists(), f"{configfile_loc} does not exist. Please specify a valid file."

    with open( configfile_loc, "r" ) as cf:
        configfile = yaml.safe_load( cf )

    schema_location = PACKAGE_DIR / f"schemas/{schema}_config.schema.yaml"

    try:
        validate( configfile, schema_location )
    except WorkflowError as err:
        # message = err.args[0].split( "\n" )
        print( "Error" )
        # sys.stderr.write( f"{message[0]} {message[1].split( ': ' )[1]}.\n" )
        sys.stderr.write( err.args[0] )
        sys.exit( -9 )

    config_paths = CONFIG_PATHS.get( schema, [] )
    for key in config_paths:
        configfile[key] = str( normalize_path( configfile[key], PACKAGE_DIR / "resources" ) )

    return configfile


def calculate_threads( threads ):
    if threads == 0:
        sys.stderr.write(
            "Pipeline cannot function without threads. Please specify a non-zero number of threads with the '--threads' "
            "option or run comman without '--threads' option to automatically detected the number of available threads.\n"
        )
        sys.exit( -4 )
    elif threads < 0:
        import multiprocessing
        threads = multiprocessing.cpu_count()
    print( f"Using {threads} threads." )

    return threads


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
    """ Searches the indicated directory for fasta files. Specifically, it will search in 'results/consensus' if the
    directory is identified as a project directory, otherwise the root directory. Additionally, the function can filter
    out fasta files which contain more than 'minimum_completeness' ambiguous characters.

    Parameters
    ----------
    directory: str
        Directory to search within. Can be a project directory.
    minimum_completeness: float
        Proportion of non-ambiguous characters a fasta file must have to be returned. Set to 0 to return all fasta
        files.

    Returns
    -------
    dict[str, Path]
        Dictionary containing sample names mapped to absolute path for each fasta file found in 'directory.'
    """
    search_directory = Path( directory ).absolute()
    if is_project_directory( search_directory ):
        search_directory = search_directory / "results/consensus"

    fastas = find_files( directory=search_directory, extensions=[".fa", ".fasta"] )

    if minimum_completeness > 0:
        for name in list( fastas.keys() ):
            completeness = calculate_completeness( fastas[name] )
            if completeness < minimum_completeness:
                del fastas[name]
    return fastas
