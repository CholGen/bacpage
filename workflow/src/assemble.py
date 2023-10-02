import argparse


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Assembles consensus sequence from raw sequencing reads."
    parser.add_argument( "--read1", help="location of read1" )
    parser.add_argument( "--read2", help="location of read2" )
    parser.add_argument( "--sample-name", help="name of sample" )

    parser.set_defaults( command=assemble_entrypoint )


def assemble_entrypoint( args: argparse.Namespace ):
    pass
