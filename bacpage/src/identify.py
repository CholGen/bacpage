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


def write_samples_to_file( sample_data: dict[str, list[str]], output: Path ):
    with open( output, "w" ) as output_file:
        output_file.write( "sample,read1,read2\n" )
        for sample in sample_data:
            files = sorted( sample_data[sample] )
            output_file.write( f"{sample},{files[0]},{files[1]}\n" )


def identify_entrypoint( args: argparse.Namespace ):
    sample_data = generate_sample_data( directory=args.directory, delim=args.delim, index=args.index )
    write_samples_to_file( sample_data=sample_data, output=args.output )
    postamble( files=sample_data, directory=Path( args.directory ), output=Path( args.output ) )


def generate_sample_data( directory: str, delim: str = "_", index: int = 0 ) -> dict[str, list[str]]:
    directory = Path( directory )
    samples = dict()
    for file in directory.iterdir():
        if file.name.endswith( ("fastq.gz", "fq.gz") ):
            sample_name = file.name.split( delim )[index]
            if sample_name in samples:
                samples[sample_name].append( file.absolute() )
            else:
                samples[sample_name] = [file.absolute()]
    for sample in samples:
        samples[sample] = sorted( samples[sample] )

    return samples


def postamble( files: dict[str, list[str]], directory: Path, output: Path ):
    print()
    print( f"Identified {len( files )} samples in {directory}." )
    print( f"Saving the location of these files to {output}." )
    print( "This project directory is now set up for assembly using `bacpage assemble`." )
