import argparse
import sys

from .src import assemble, example, identify, phylogeny

COMMANDS = {
    "assemble"      : [assemble, "Assembles consensus sequence from raw sequencing reads."],
    "example"       : [example, "Set up project directory for analysis."],
    "identify_files": [identify, "Generate a valid sample_data.csv from a directory of FASTQs."],
    "phylogeny"     : [phylogeny, "Align sequences and construct a maximum likelihood tree."],
    # "submit" : [submit, "Prepare files for submission to online repositories."],
    # "typing" : [typing, "Classify consensus sequences based on the presence or absense of various genes."],
}


def main( sysargs=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        prog="bacpage",
        description="""██████╗  █████╗  ██████╗██████╗  █████╗  ██████╗ ███████╗
██╔══██╗██╔══██╗██╔════╝██╔══██╗██╔══██╗██╔════╝ ██╔════╝
██████╔╝███████║██║     ██████╔╝███████║██║  ███╗█████╗
██╔══██╗██╔══██║██║     ██╔═══╝ ██╔══██║██║   ██║██╔══╝
██████╔╝██║  ██║╚██████╗██║     ██║  ██║╚██████╔╝███████╗
╚═════╝ ╚═╝  ╚═╝ ╚═════╝╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚══════╝

    A bioinformatics toolkit to assemble and analyze BACterial PAthogen GEnomes""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(
        title="Available commands",
        description="One of the following commands must be specified:",
        required=True
    )

    for command, values in COMMANDS.items():
        parser_subcommand = subparsers.add_parser( command, help=values[1] )
        values[0].add_command_arguments( parser_subcommand )

    if len( sysargs ) < 1:
        parser.print_help()
        sys.exit( -1 )
    else:
        args = parser.parse_args( sysargs )
    args.command( args )


if __name__ == "__main__":
    main()
