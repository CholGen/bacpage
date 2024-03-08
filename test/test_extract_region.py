import shutil
from pathlib import Path

import pytest
from Bio import AlignIO

from bacpage.src.utils import extract_region


@pytest.fixture
def extraction_run():
    regions = "test/test_extract_script/regions.bed"
    directory = "test/test_extract_script/"
    directory_path = Path( directory )

    extract_region.extract_regions( regions=regions, directory=directory )
    yield directory_path / extract_region.OUTPUT_PATH

    if (directory_path / extract_region.OUTPUT_PATH).exists():
        shutil.rmtree( directory_path / extract_region.OUTPUT_PATH )


def tests_output_directory_created( extraction_run: Path ):
    assert extraction_run.exists(), f"Output directory {extraction_run} was not created."


def test_correct_number_of_regions_extracted( extraction_run: Path ):
    files = [i for i in extraction_run.iterdir() if i.suffix == ".fasta"]
    got = len( files )
    want = 2
    assert got == want, f"Incorrect number of regions extracted. Want {want}, got {got}."


def test_results_are_valid_fastas( extraction_run: Path ):
    for file in extraction_run.iterdir():
        if file.suffix == ".fasta":
            records = AlignIO.read( file, "fasta" )


def test_correct_number_of_sequences_extracted( extraction_run: Path ):
    number_of_sequences = []
    expected_number_of_sequencs = 4
    for file in extraction_run.iterdir():
        if file.suffix == ".fasta":
            records = AlignIO.read( file, "fasta" )
            number_of_sequences.append( (file, len( records )) )
    assert all( map( lambda x: x[1] == expected_number_of_sequencs,
                     number_of_sequences ) ), f"Incorrect number of sequences in a result: {number_of_sequences}"


def test_correct_extraction_length_ctxA( extraction_run: Path ):
    ctxA = AlignIO.read( extraction_run / "ctxA.extracts.fasta", "fasta" )
    got = ctxA.get_alignment_length()
    want = 776
    assert got == want


def test_correct_extraction_length_ctxB( extraction_run: Path ):
    ctxb = AlignIO.read( extraction_run / "ctxB.extracts.fasta", "fasta" )
    got = ctxb.get_alignment_length()
    want = 374
    assert got == want
