import argparse
from pathlib import Path

from jinja2 import Environment, FileSystemLoader


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Set up project directory for analysis."
    parser.add_argument( "directory", help="Location to create project directory" )
    parser.add_argument( "--quiet", action="store_true", help="Do not display helpful messages during creation." )
    parser.set_defaults( command=example_entrypoint )


def example_entrypoint( args: argparse.Namespace ):
    create_project_directory( args.directory, quiet=args.quiet )


def create_project_directory( directory: str, quiet: bool = False ):
    project_directory = Path( directory ).absolute()
    assert not project_directory.exists(), f"{project_directory} already exists."

    # create project directory
    project_directory.mkdir()

    # create input directory
    (project_directory / "input").mkdir()

    # create samples file
    sample_data = project_directory / "sample_data.csv"
    with sample_data.open( "w" ) as samples_file:
        samples_file.write( "sample,read1,read2\n" )
        samples_file.write( "a,path-to-a-read1,path-to-a-read2\n" )
        samples_file.write( "b,path-to-b-read1,path-to-b-read2\n" )
        samples_file.write( "c,path-to-c-read1,path-to-c-read2\n" )

    # create config file
    environment = Environment( loader=FileSystemLoader( "workflow/schemas/" ) )
    template = environment.get_template( "illumina_config.template.yaml" )

    project_config = dict()
    project_config["sample_data"] = sample_data
    project_config['project_path'] = project_directory
    content = template.render( project_config )
    with (project_directory / "config.yaml").open( "w" ) as config_file:
        config_file.write( content )

    if not quiet:
        print( f"Creating project directory at {directory}." )
        print( f"The absolute path of the project directory is {project_directory}" )
        print()
        print( "The structure of the project directory is as follows:" )
        print( project_directory.name + "/" )
        print( "├── input/" )
        print( "├── sample_data.csv" )
        print( "└── config.yaml" )
        print( "" )
        print( "To assemble genomes, place raw sequencing FASTQs in the input directory." )
        print( "Next, use `bacpage identify` to fill in sample_data.csv with information about your samples." )
        print( "Finally, use `bacpage assemble` to generate consensus sequences." )
