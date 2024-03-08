import argparse
import sys
from pathlib import Path

from Bio import SeqIO

from bacpage.src import common_funcs

OUTPUT_PATH = "results/extractions/"


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Extract regions specified in a bed file from consensus sequences in a project directory."
    parser.add_argument(
        "directory", default=".", help="location of FASTQ files [current directory]"
    )
    parser.add_argument( "--region", help="BED file containing region(s) to extract." )
    parser.add_argument( "--output", help="Directory to save extraction alignments" )

    parser.set_defaults( command=extract_entrypoint )


def parse_bed_file( region ):
    region_path = Path( region ).absolute()

    if not region_path.exists():
        sys.stderr.write( f"{region_path} does not exist. Please specify a valid BED file with the '--region' option." )
        sys.exit( -1 )

    parsed_regions = dict()
    with open( region_path, "r" ) as region_file:
        for line in region_file:
            fields = line.strip().split( "\t" )
            try:
                name = fields[0]
                start = int( fields[1] )
                end = int( fields[2] )
            except IndexError:
                sys.stderr.write(
                    f"Unable to parse line: '{line}' in bed file due to not enough fields. Make sure {region_file} contains at least three fields: name, start, and end"
                )
                sys.exit( -2 )

            parsed_regions[name] = [start, end]
    return parsed_regions


def extract_regions( regions: str, directory: str ):
    # Load bed file
    bed = parse_bed_file( regions )
    region_names = list( bed.keys() )

    # Load fasta files
    input_sequences = common_funcs.load_input( directory=directory, minimum_completeness=0 )

    # Iterate through regions
    extractions = {name: [] for name in region_names}
    for sample, sequence_file in input_sequences.items():
        sequence = SeqIO.read( sequence_file, "fasta" )
        for name, region in bed.items():
            extraction = sequence[region[0]:region[1]]
            extractions[name].append( extraction )

    for name, extracts in extractions.items():
        outdir_path = Path( directory ).absolute() / OUTPUT_PATH
        outdir_path.mkdir( exist_ok=True )
        SeqIO.write( extracts, outdir_path / f"{name}.extracts.fasta", "fasta" )

    postamble( regions=region_names, directory=Path( directory ) )


def extract_entrypoint( args: argparse.Namespace ):
    extract_regions(
        regions=args.region,
        directory=args.directory
    )


def postamble( regions: list[str], directory: Path ):
    print()
    print( f"Extracted regions {', '.join( regions )} from files in {directory}." )
    print( f"Multi-sequence alignments for each regions have been saved to {directory / OUTPUT_PATH}" )
