import argparse
import sys

from .src import assemble, example, identify, phylogeny, profiling, version
from .src.utils import extract_region

COMMANDS = {
    "assemble"      : [assemble, "Assembles consensus sequence from raw sequencing reads."],
    "setup"         : [example, "Set up project directory for analysis."],
    "identify_files": [identify, "Generate a valid sample_data.csv from a directory of FASTQs."],
    "phylogeny"     : [phylogeny, "Align sequences and construct a maximum likelihood tree."],
    "profile"       : [profiling, "Classify consensus sequences based on the presence or absense of various genes."],
    "version"       : [version, "Prints the version of bacpage and exits."]
    # "submit" : [submit, "Prepare files for submission to online repositories."],
}

UTILITIES = {
    "extract_region": [extract_region, "Extract substring(s) from all consensus sequences."]
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

    util_command = subparsers.add_parser( "utilities", help="Miscellaneous tools." )
    util_subparser = util_command.add_subparsers(
        dest="command",
        title="Available utilities",
        description="One of the following utilities must be specified:"
    )
    for command, values in UTILITIES.items():
        util_subcommand = util_subparser.add_parser( command, help=values[1] )
        values[0].add_command_arguments( util_subcommand )

    if len( sysargs ) < 1:
        parser.print_help()
        sys.exit( -1 )
    else:
        args = parser.parse_args( sysargs )

    if not args.command:
        util_command.print_help()
        sys.exit( -1 )

    args.command( args )


if __name__ == "__main__":
    main()
