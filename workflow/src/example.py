import argparse
import os

import yaml

DEFAULT_CONFIG = {
    'run_type': 'Illumina',
    'samples': '<project-path>/sample_data.csv',
    'output_directory': '<project-path>',
    'reference': '<sequencing-path>/resources/vc_reference.fasta',
    'reference_genes': '<sequencing-path>/resources/cholera_ref_genes/',
    'recombinant_mask': '<sequencing-path>/resources/cholera_gubbins_mask.gff',
    'background_dataset': '<background-dataset-path>',
    'generate': {
        'consensus_sequences': True,
        'typing': False,
        'quality_control_report': True,
        'phylogeny': False
    },
    'preprocessing': {
        'check_size': True,
        'minimum_size': 100
    },
    'alignment_bwa': {
        'bwa_params': '-M'
    },
    'trimming': {
        'minimum_length': 30,
        'minimum_quality': 20,
        'window_length': 4
    },
    'coverage_mask': {
        'required_depth': 10
    },
    'plot_coverage': {
        'bin_size': 10000
    },
    'call_variants': {
        'maximum_depth': 2000,
        'minimum_mapping_quality': 30,
        'minimum_base_quality': 20,
        'mpileup_parameters': '-B -a INFO/AD,INFO/ADF,INFO/ADR -Ou',
        'call_parameters': '-mv -Ov --ploidy 1'
    },
    'filter_variants': {
        'minimum_depth': 10,
        'minimum_support': 0.5,
        'minimum_strand_depth': 5
    },
    'call_consensus': { 'consensus_parameters': '--mark-del N' },
    'mlst_profiling': {
        'scheme': 'vcholerae',
        'mlst_params': '--quiet --csv --legacy'
    },
    'antibiotic_resistance': { 'database': 'card' },
    'tree_building': {
        'minimum_completeness': 0.9,
        'outgroup': 'Asia|IDN|ERR025382|UNK|1957',
        'iqtree_parameters': '-nt AUTO -m TEST -bb 1000'
    }
}


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Set up project directory for analysis."
    parser.add_argument( "--destination", required=True, help="Location to create project directory" )
    parser.set_defaults( command=example_entrypoint )


def example_entrypoint( args: argparse.Namespace ):
    create_project_directory( args.destination )


def create_project_directory( directory: str ):
    assert not os.path.exists( directory ), f"{directory} already exists."

    # create project directory
    os.mkdir( directory )

    # create input directory
    os.mkdir( os.path.join( directory, "input" ) )

    # create samples file
    with open( os.path.join( directory, "sample_data.csv" ), "w" ) as samples_file:
        samples_file.write( "sample,read1,read2\n" )
        samples_file.write( "a,path-to-a-read1,path-to-a-read2\n" )
        samples_file.write( "b,path-to-b-read1,path-to-b-read2\n" )
        samples_file.write( "c,path-to-c-read1,path-to-c-read2\n" )

    # create config file
    project_config = DEFAULT_CONFIG
    project_config["samples"] = os.path.join( directory, "sample_data.csv" )
    project_config['output_directory'] = directory
    with open( os.path.join( directory, "config.yaml" ), "w" ) as config_file:
        yaml.dump( project_config, config_file )
