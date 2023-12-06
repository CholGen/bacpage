import shutil
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import snakemake
from Bio import Phylo

from bacpage.src import phylogeny, common_funcs


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


def test_config_finds_local_parameters():
    config = common_funcs.load_configfile(
        "test/configs/local_parameters.yaml", Path( "test/test_pipeline" ).absolute(), schema="phylogeny"
    )
    reference = Path( config["reference"] )
    assert reference.is_absolute() and reference.exists(), f"Local reference in resources directory was not correctly found ({reference})."
    recombinant_mask = Path( config["recombinant_mask"] )
    assert recombinant_mask.is_absolute() and recombinant_mask.exists(), f"Local reference gene directory was not correctly found ({recombinant_mask})."


def test_error_if_vcf_background_doesnt_match_reference():
    search_directory = "test/test_tree_vcf"
    with pytest.raises( SystemExit ) as excinfo:
        phylogeny.reconstruct_phylogeny(
            project_directory=search_directory,
            configfile=".",
            minimum_completeness=0,
            threads=1,
            verbose=True
        )
    assert excinfo.value.code < 0


def test_raise_error_if_reference_present_already():
    search_directory = "test/test_tree_project_directory"
    with pytest.raises( SystemExit ) as excinfo:
        phylogeny.reconstruct_phylogeny(
            project_directory=search_directory,
            configfile="test/test_tree_project_directory/duplicate_reference.yaml",
            minimum_completeness=0,
            threads=1,
            verbose=True,
            dryrun=True
        )
    assert excinfo.value.code < 0


def test_raise_error_with_duplicate_sample_names():
    search_directory = "test/test_tree_duplicate_sequences"
    with pytest.raises( SystemExit ) as excinfo:
        phylogeny.reconstruct_phylogeny(
            project_directory=search_directory,
            configfile="test/test_tree_fasta_directory/config.yaml",
            minimum_completeness=0,
            threads=1,
            verbose=True
        )
    assert excinfo.value.code == -6


def test_parse_vcf_names():
    vcf = Path( "test/test_vcf.vcf" )
    found = phylogeny.parse_names_vcf( vcf )
    expected = [
        "Africa|TZA|SAMN19110433|T10|2017",
        "Africa|TZA|SAMN19110437|T13|2016",
        "Africa|TZA|SAMN19110438|T13|2016",
        "Africa|TZA|SAMN19110441|T13|2017",
        "Africa|TZA|SAMN19110439|T13|2016",
        "Africa|TZA|SAMN19110428|T13|2017",
        "Africa|TZA|SAMN19110430|T13|2017",
        "Africa|TZA|SAMN19110431|T13|2017",
        "Africa|TZA|SAMN19110436|T13|2017",
        "Africa|TZA|SAMN19110443|T13|2017"
    ]
    assert sorted( found ) == sorted( expected )


def test_vcf_to_fasta():
    from bacpage.scripts.vcf_to_fasta import convert_vcf
    import filecmp
    vcf_file = "test/test_vcf_to_fasta/input.vcf.gz"
    reference = "test/test_vcf_to_fasta/reference.fasta"
    output = "test/test_vcf_to_fasta/output.fasta"

    convert_vcf( vcf=vcf_file, reference=reference, output=output )
    assert Path( output ).exists() and Path( output ).is_file(), f"{output} does not exist or is not a file."
    assert filecmp.cmp( "test/test_vcf_to_fasta/input.fasta", output, shallow=True )


def get_rules_dryrun( snakefile: Path, config: dict[str, Any], workdir: str ):
    with redirect_stdout( StringIO() ) as f:
        value = snakemake.snakemake( snakefile, forceall=True, workdir=workdir, config=config, summary=True )
    df = pd.read_csv( StringIO( f.getvalue() ), sep="\t" )
    return df["rule"].unique()


def test_correct_rules_run_if_fasta_background():
    search_directory = "test/test_tree_project_directory"
    config, snakefile = phylogeny.reconstruct_phylogeny(
        project_directory=search_directory,
        configfile="test/test_tree_project_directory/background_fasta.yaml",
        minimum_completeness=0,
        threads=1,
        verbose=True,
        dryrun=True,
    )
    estimated_rules = get_rules_dryrun( snakefile, config, search_directory )
    expected_rules = ["concatenate_sequences", "concatenate_reference", "convert_to_vcf", "generate_alignment_from_vcf",
                      "sparsify_alignment", "run_gubbins", "generate_tree", "move_tree_and_rename",
                      "move_recombinant_mask"]

    assert sorted( estimated_rules ) == sorted( expected_rules )


def test_correct_rules_run_if_vcf_background():
    search_directory = "test/test_tree_project_directory"
    config, snakefile = phylogeny.reconstruct_phylogeny(
        project_directory=search_directory,
        configfile="test/test_tree_project_directory/vcf_background.yaml",
        minimum_completeness=0,
        threads=1,
        verbose=True,
        dryrun=True,
    )
    estimated_rules = get_rules_dryrun( snakefile, config, search_directory )
    expected_rules = ["concatenate_sequences", "concatenate_reference", "convert_to_vcf",
                      "combine_sequences_and_background_vcf", "generate_alignment_from_vcf", "run_gubbins",
                      "sparsify_alignment", "generate_tree", "move_tree_and_rename", "move_recombinant_mask"]

    assert sorted( estimated_rules ) == sorted( expected_rules )


def test_correct_rules_run_if_masking_specified():
    search_directory = "test/test_tree_project_directory"
    config, snakefile = phylogeny.reconstruct_phylogeny(
        project_directory=search_directory,
        configfile=".",
        minimum_completeness=0,
        threads=1,
        mask="bacpage/resources/cholera_mask.gff",
        verbose=True,
        dryrun=True,
    )
    estimated_rules = get_rules_dryrun( snakefile, config, search_directory )
    expected_rules = ["concatenate_sequences", "concatenate_reference", "convert_to_vcf", "mask_vcf",
                      "generate_alignment_from_vcf", "run_gubbins", "sparsify_alignment", "generate_tree",
                      "move_recombinant_mask", "move_tree_and_rename"]

    assert sorted( estimated_rules ) == sorted( expected_rules )


def test_correct_rules_run_if_recombinant_masking_skipped():
    search_directory = "test/test_tree_project_directory"
    config, snakefile = phylogeny.reconstruct_phylogeny(
        project_directory=search_directory,
        configfile=".",
        minimum_completeness=0,
        threads=1,
        skip_detect=True,
        verbose=True,
        dryrun=True,
    )
    estimated_rules = get_rules_dryrun( snakefile, config, search_directory )
    expected_rules = ["concatenate_sequences", "concatenate_reference", "convert_to_vcf",
                      "generate_alignment_from_vcf", "sparsify_alignment", "generate_tree",
                      "move_tree_and_rename"]

    assert sorted( estimated_rules ) == sorted( expected_rules )


def test_correct_rules_run_if_terra_specified():
    search_directory = "test/test_tree_project_directory"
    config, snakefile = phylogeny.reconstruct_phylogeny(
        project_directory=search_directory,
        configfile=".",
        minimum_completeness=0,
        threads=1,
        terra=True,
        verbose=True,
        dryrun=True,
    )
    estimated_rules = get_rules_dryrun( snakefile, config, search_directory )
    expected_rules = ["concatenate_sequences", "concatenate_reference", "convert_to_vcf"]

    assert sorted( estimated_rules ) == sorted( expected_rules )


def test_correct_rules_run_if_gubbins_mask_supplied():
    search_directory = "test/test_tree_project_directory"
    mask = Path( "bacpage/resources/cholera_gubbins_mask.gff" ).absolute()
    config, snakefile = phylogeny.reconstruct_phylogeny(
        project_directory=search_directory,
        configfile=".",
        minimum_completeness=0,
        threads=1,
        detect_file=str( mask ),
        verbose=True,
        dryrun=True,
    )
    estimated_rules = get_rules_dryrun( snakefile, config, search_directory )
    expected_rules = ["concatenate_sequences", "concatenate_reference", "convert_to_vcf",
                      "bypass_gubbins", "generate_alignment_from_bypassed_gubbins", "sparsify_alignment",
                      "generate_tree", "move_tree_and_rename"]
    assert sorted( estimated_rules ) == sorted( expected_rules )


def test_abridged_config_accepted():
    search_directory = "test/test_tree_project_directory"
    config = "test/configs/phylogeny_only.yaml"
    config, snakefile = phylogeny.reconstruct_phylogeny(
        project_directory=search_directory,
        configfile=config,
        dryrun=True
    )
    assert config is not None


@pytest.fixture
def phylogeny_run( scope="session" ):
    project_directory = Path( "test/test_tree_fasta_directory" )
    phylogeny.reconstruct_phylogeny(
        str( project_directory ), ".", minimum_completeness=0.9, threads=-1,
        verbose=False
    )
    yield project_directory

    if (project_directory / "results").exists():
        shutil.rmtree( project_directory / "results" )
    if (project_directory / "intermediates").exists():
        shutil.rmtree( project_directory / "intermediates" )


def test_phylogeny_postamble():
    project_directory = Path( "test/test_tree_fasta_directory" ).absolute()
    phylogeny.postamble( project_directory )


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
