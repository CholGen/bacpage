import argparse
from pathlib import Path


def add_command_arguments( parser: argparse.ArgumentParser ):
    parser.description = "Generate a valid sample_data.csv from a directory of FASTQs."
    parser.add_argument(
        "directory", default=".", help="location of FASTQ files [current directory]"
    )
    parser.add_argument(
        "--delim", default="_", type=str, help="deliminator to extract sample name from file name [_]"
    )
    parser.add_argument(
        "--index", default=0, type=int, help="index of sample name after splitting file name by delim [0]"
    )
    parser.add_argument(
        "--output", type=str, default="sample_data.csv", help="location to save sample data ['sample_data.csv']"
    )

    parser.set_defaults( command=identify_entrypoint )


def identify_entrypoint( args ):
    generate_sample_data( directory=args.directory, output=args.output, delim=args.delim, index=args.index )


def generate_sample_data( directory, output, delim="_", index=0 ):
    directory = Path( directory )
    samples = dict()
    for file in directory.iterdir():
        if file.name.endswith( ("fastq.gz", "fq.gz") ):
            sample_name = file.name.split( delim )[index]
            if sample_name in samples:
                samples[sample_name].append( file.absolute() )
            else:
                samples[sample_name] = [file.absolute()]

    with open( output, "w" ) as output_file:
        output_file.write( "sample,read1,read2\n" )
        for sample in samples:
            files = sorted( samples[sample] )
            output_file.write( f"{sample},{files[0]},{files[1]}\n" )
