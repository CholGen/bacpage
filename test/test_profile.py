from pathlib import Path

import pytest

from bacpage.src import profiling


def test_find_folder_of_fastas():
    search_directory = Path( "test/test_tree_fasta_directory" ).absolute()
    found = profiling.find_sequences( directory=search_directory )
    search_directory = Path( search_directory ).absolute()
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 21 )]
    assert sorted( found.values() ) == sorted( expected )


def test_find_consensus_sequences():
    search_directory = Path( "test/test_tree_project_directory" ).absolute()
    found = profiling.find_sequences( directory=search_directory )
    search_directory = Path( search_directory ).absolute() / "results/consensus"
    expected = [search_directory / f"t{i}.fasta" for i in range( 1, 21 )]
    assert sorted( found.values() ) == sorted( expected )


def test_find_assemblies():
    search_directory = Path( "test/test_find_assemblies" ).absolute()
    found = profiling.find_sequences( directory=search_directory )
    search_directory = Path( search_directory ).absolute() / "results/assembly"
    expected = [search_directory / f"t{i}.annotated.gff" for i in range( 1, 10 )]
    assert sorted( found.values() ) == sorted( expected )


def test_raise_error_with_duplicate_samples_in_directory():
    search_directory = Path( "test/test_tree_duplicate_samples" ).absolute()
    with pytest.raises( SystemExit ) as excinfo:
        found = profiling.find_sequences( directory=search_directory )
    assert excinfo.value.code == -5


def test_profiling_snakemake_works():
    project_directory = Path( "test/test_profiling/" )
    expected_output = [project_directory / "results/reports/antibiotic_resistance.tsv"]
    profiling.profile_sequences(
        project_directory,
        ".",
        "card", -1,
        False
    )

    results = dict()
    for file in expected_output:
        results[file] = file.exists()
    assert all( results.values() ), f"Not all expected files generated. Expected {results}"

    expected_output[0].unlink( missing_ok=True )
