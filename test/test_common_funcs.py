from pathlib import Path

import pytest

from bacpage.src import common_funcs


def test_recognize_folder_of_fastas():
    search_directory = "test/test_tree_fasta_directory"
    found = common_funcs.load_input( directory=search_directory, minimum_completeness=0 )
    search_directory = Path( search_directory ).absolute()
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 21 )]
    assert sorted( found.values() ) == sorted( expected )


def test_recognize_project_directory():
    search_directory = "test/test_tree_project_directory"
    found = common_funcs.load_input( directory=search_directory, minimum_completeness=0 )
    search_directory = Path( search_directory ).absolute() / "results/consensus"
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 21 )]
    assert sorted( found.values() ) == sorted( expected )


def test_remove_low_coverage_sequences():
    search_directory = "test/test_tree_fasta_directory"
    found = common_funcs.load_input( directory=search_directory, minimum_completeness=0.9 )
    search_directory = Path( search_directory ).absolute()
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 15 )]
    assert sorted( found.values() ) == sorted( expected )


def test_raise_error_with_duplicate_samples_in_directory():
    search_directory = "test/test_tree_duplicate_samples"
    with pytest.raises( SystemExit ) as excinfo:
        found = common_funcs.load_input( directory=search_directory, minimum_completeness=0.9 )
    assert excinfo.value.code == -5


def test_config_finds_local_parameters():
    config = common_funcs.load_configfile(
        "test/configs/local_parameters.yaml", Path( "test/test_pipeline" ).absolute(), schema="phylogeny"
    )
    reference = Path( config["reference"] )
    assert reference.is_absolute() and reference.exists(), f"Local reference in resources directory was not correctly found ({reference})."
    recombinant_mask = Path( config["recombinant_mask"] )
    assert recombinant_mask.is_absolute() and recombinant_mask.exists(), f"Local reference gene directory was not correctly found ({recombinant_mask})."
