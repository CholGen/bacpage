import argparse


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Set up project directory for analysis."
    parser.add_argument( "--destination", required=True, help="Location to create project directory" )
    parser.set_defaults( command=example_entrypoint )


def example_entrypoint( args: argparse.Namespace ):
    pass
