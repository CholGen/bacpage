import shutil
from pathlib import Path

import pytest
from Bio import Phylo

from bacpage.src import phylogeny


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


@pytest.fixture()
def phylogeny_run( scope="session" ):
    project_directory = Path( "test/test_tree_fasta_directory" )
    phylogeny.reconstruct_phylogeny( str( project_directory ), ".", minimum_completeness=0.9, threads=-1,
                                     verbose=False )
    yield project_directory

    if (project_directory / "results").exists():
        shutil.rmtree( project_directory / "results" )
    if (project_directory / "intermediates").exists():
        shutil.rmtree( project_directory / "intermediates" )


@pytest.mark.slow
def test_tree_reconstruction_successfully( phylogeny_run ):
    tree = phylogeny_run / "results/phylogeny/phylogeny.tree"
    assert tree.exists() and tree.is_file(), "Phylogeny was either not created or is not a file."


@pytest.mark.slow
def test_tree_reconstruction_all_taxa_present( phylogeny_run ):
    tree = phylogeny_run / "results/phylogeny/phylogeny.tree"
    tree = Phylo.read( tree, "newick" )

    found = [clade.name for clade in tree.get_terminals()]
    expected = [f"t{i}" for i in range( 1, 15 )]

    assert sorted( found ) == sorted( expected ), "Not all taxa where found in tree."
    assert 1 < tree.total_branch_length() < 100, f"Branch length ({tree.total_branch_length()}) is not a reasonable magnitude."
