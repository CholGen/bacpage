import argparse

from bacpage import __version__


def print_version( args: argparse.Namespace ):
    print( __version__ )


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.set_defaults( command=print_version )
