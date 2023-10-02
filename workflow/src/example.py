import argparse
import os

from jinja2 import Environment, FileSystemLoader


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
    environment = Environment( loader=FileSystemLoader( "workflow/schemas/" ) )
    template = environment.get_template( "illumina_config.template.yaml" )

    project_config = dict()
    project_config["sample_data"] = os.path.abspath( os.path.join( directory, "sample_data.csv" ) )
    project_config['project_path'] = os.path.abspath( directory )
    content = template.render( project_config )
    with open( os.path.join( directory, "config.yaml" ), "w" ) as config_file:
        config_file.write( content )
