import sys
from importlib.resources import files
from pathlib import Path

import yaml
from snakemake.utils import validate

DEFAULT_CONFIG = "config.yaml"
DEFAULT_SAMPLEDATA = "sample_data.csv"
PACKAGE_DIR = files( "bacpage" )
CONFIG_PATHS = ["reference", "reference_genes", "recombinant_mask"]
RESTART_TIMES = 0


def find_files( directory: Path, extensions: list[str] ) -> dict[str, Path]:
    search_path = directory.absolute()
    found = dict()
    for file in search_path.iterdir():
        if file.suffix in extensions:
            name = file.stem
            if name in found:
                sys.stderr.write( f"{name} was found in two or more files. Duplicate sequences cannot be processed." )
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


def load_configfile( specified_loc: str, project_directory: Path ) -> dict:
    configfile_loc = Path( specified_loc ).absolute()
    if specified_loc == ".":
        configfile_loc = project_directory / DEFAULT_CONFIG
        assert configfile_loc.exists(), "Unable to automatically find config in project directory (searching for 'config.yaml'). Please specify a valid configuration file."
    assert configfile_loc.exists(), f"{configfile_loc} does not exist. Please specify a valid file."

    with open( configfile_loc, "r" ) as cf:
        configfile = yaml.safe_load( cf )

    schema_location = PACKAGE_DIR / "schemas/Illumina_config.schema.yaml"
    validate( configfile, schema_location )

    for key in CONFIG_PATHS:
        configfile[key] = str( normalize_path( configfile[key], PACKAGE_DIR / "resources" ) )

    return configfile


def calculate_threads( threads ):
    if threads == 0:
        sys.stderr.write(
            "Pipeline cannot function without threads. Please specify a non-zero number of threads with the '--threads' "
            "option or run comman without '--threads' option to automatically detected the number of available threads."
        )
        sys.exit( -4 )
    elif threads < 0:
        import multiprocessing
        threads = multiprocessing.cpu_count()
    print( f"Using {threads} threads." )
    return threads
