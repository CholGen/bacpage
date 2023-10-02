import os.path

import pandas as pd
import pytest
import yaml
from snakemake.utils import validate

from workflow.src import example


@pytest.fixture( scope="session" )
def project_directory( tmp_path_factory ) -> str:
    directory = tmp_path_factory.mktemp( "testing_directory" )
    project_directory = os.path.join( directory, "project" )
    example.create_project_directory( directory=project_directory )
    return project_directory


def test_creation_of_project_directory( project_directory: str ):
    assert os.path.exists( project_directory ), f"{project_directory} does not exist."
    assert os.path.isdir( project_directory ), f"{project_directory} exists but is not a directory."


def test_creation_of_input_directory( project_directory: str ):
    input_directory = os.path.join( project_directory, "input" )
    assert os.path.exists( input_directory ), f"{input_directory} does not exist."
    assert os.path.isdir( input_directory ), f"{input_directory} exists but is not a directory."


def test_config_file_exists( project_directory: str ):
    config_file = os.path.join( project_directory, "config.yaml" )
    assert os.path.exists( config_file ), f"{config_file} does not exist."
    assert os.path.isfile( config_file ), f"{config_file} exists but is not a file."
    assert os.path.getsize( config_file ) > 100, f"{config_file} should contain text but is empty."


def test_config_file_is_valid_yaml( project_directory: str ):
    config_file = os.path.join( project_directory, "config.yaml" )
    with open( config_file, "r" ) as parameters:
        project_yaml = yaml.safe_load( parameters )
        validate( project_yaml, "workflow/schemas/Illumina_config.schema.yaml" )


def test_sample_file_exists( project_directory: str ):
    samples_file = os.path.join( project_directory, "sample_data.csv" )
    assert os.path.exists( samples_file ), f"{samples_file} does not exist."
    assert os.path.isfile( samples_file ), f"{samples_file} exists but is not a file."
    assert os.path.getsize( samples_file ) > 100, f"{samples_file} should contain text but is empty."


def test_sample_file_is_valid_csv( project_directory: str ):
    sf = pd.read_csv( os.path.join( project_directory, "sample_data.csv" ) )
    for col in ["sample", "read1", "read2"]:
        assert col in sf.columns, f"{col} not in columns."
    assert sf.shape == (3, 3)
