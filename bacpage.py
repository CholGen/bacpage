import argparse

from workflow.src import assemble, example, identify

COMMANDS = {
    "assemble": [assemble, "Assembles consensus sequence from raw sequencing reads."],
    # "phylogeny" : [phylogeny, "Align sequences and construct a maximum likelihood tree."],
    # "typing" : [typing, "Classify consensus sequences based on the presence or absense of various genes."],
    # "submit" : [submit, "Prepare files for submission to online repositories."],
    "identify_files": [identify, "Generate a valid sample_data.csv from a directory of FASTQs."],
    "example": [example, "Set up project directory for analysis."]
}

if __name__ == "__main__":
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
        description="One of the following commands must be specified:"
    )

    for command, values in COMMANDS.items():
        parser_subcommand = subparsers.add_parser( command, help=values[1] )
        values[0].add_command_arguments( parser_subcommand )

    args = parser.parse_args()
    args.command( args )
