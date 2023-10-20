import shutil
from pathlib import Path

import pytest

from workflow.src import phylogeny


def test_recognize_folder_of_fastas():
    search_directory = "test/test_tree_fasta_directory"
    found = phylogeny.load_input( directory=search_directory, minimum_completeness=0 )
    search_directory = Path( search_directory ).absolute()
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 21 )]
    assert sorted( found.values() ) == sorted( expected )


def test_recognize_project_directory():
    search_directory = "test/test_tree_project_directory"
    found = phylogeny.load_input( directory=search_directory, minimum_completeness=0 )
    search_directory = Path( search_directory ).absolute() / "results/consensus"
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 21 )]
    assert sorted( found.values() ) == sorted( expected )


def test_remove_low_coverage_sequences():
    search_directory = "test/test_tree_fasta_directory"
    found = phylogeny.load_input( directory=search_directory, minimum_completeness=0.9 )
    search_directory = Path( search_directory ).absolute()
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 15 )]
    assert sorted( found.values() ) == sorted( expected )


def test_raise_error_with_duplicate_samples_in_directory():
    search_directory = "test/test_tree_duplicate_samples"
    with pytest.raises( SystemExit ) as excinfo:
        found = phylogeny.load_input( directory=search_directory, minimum_completeness=0.9 )
    assert excinfo.value.code == -5


def test_raise_error_with_duplicate_sample_names():
    search_directory = "test/test_tree_duplicate_sequences"
    with pytest.raises( SystemExit ) as excinfo:
        phylogeny.reconstruct_phylogeny( search_directory, "test/test_tree_fasta_directory/config.yaml", 0, 1, True )
    assert excinfo.value.code == -6


# @pytest.fixture()
# def build_tree( scope="session" ):
#    pass
@pytest.mark.slow
def test_build_tree_successfully():
    project_directory = Path( "test/test_tree_fasta_directory" )
    phylogeny.reconstruct_phylogeny( str( project_directory ), ".", minimum_completeness=0.9, threads=-1, verbose=True )

    # Cleanup
    if (project_directory / "results").exists():
        shutil.rmtree( project_directory / "results" )
    if (project_directory / "intermediates").exists():
        shutil.rmtree( project_directory / "intermediates" )
