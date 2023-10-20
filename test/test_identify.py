import os
from pathlib import Path

import pandas as pd
import pytest

from workflow.src import identify

SAMPLES = [
    "ERR976393",
    "MOZ-PMB0645689",
    "DRC-HL423",
    "CMR-CEN020JA",
    "ERR976399",
    "ERR976398"
]


@pytest.fixture( scope="session" )
def sample_data( tmp_path_factory ) -> str:
    directory = tmp_path_factory.mktemp( "testing_directory" )
    sample_data = os.path.join( directory, "sample_data.csv" )
    identify.generate_sample_data( "test/test_fastqs", sample_data )
    return sample_data


def test_if_exists( sample_data ):
    assert os.path.exists( sample_data ), f"{sample_data} does not exist."
    assert os.path.isfile( sample_data ), f"{sample_data} exists but is not a file."
    assert os.path.getsize( sample_data ) > 100, f"{sample_data} should contain text but is empty."


def test_if_valid_csv( sample_data ):
    sample_df = pd.read_csv( sample_data )
    for col in ["sample", "read1", "read2"]:
        assert col in sample_df.columns, f"{col} not in columns."
    assert sample_df.shape[1] > 0, "sample_data.csv contains no entries."


def test_if_all_samples_parsed( sample_data ):
    sample_df = pd.read_csv( sample_data )
    results = dict()
    for sample in SAMPLES:
        results[sample] = any( sample_df["sample"] == sample )
    not_found = [sample for sample in results if not results[sample]]
    assert len( not_found ) == 0, f"Not all samples found in samples_data.csv ({', '.join( not_found )})"


def test_if_sample_path_valid( sample_data ):
    sample_df = pd.read_csv( sample_data )
    for idx, entry in sample_df.iterrows():
        read1 = Path( entry["read1"] )
        read2 = Path( entry["read2"] )
        assert read1.exists(), f"{read1} does not exist."
        assert read1.is_absolute(), f"{read1} should be absolute but isn't."
        assert read2.exists(), f"{read2} does not exist."
        assert read2.is_absolute(), f"{read2} should be absolute but isn't."
