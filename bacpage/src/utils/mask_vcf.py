import argparse
import gzip
import io
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


def add_command_arguments( pars: argparse.ArgumentParser ):
    pars.description = "Mask VCF file based on recombinant regions detected by Gubbins."
    pars.add_argument(
        "--regions", help="GFF file containing recombinant regions detected by Gubbins", required=True, type=Path
    )
    pars.add_argument( "--vcf", help="VCF file to be masked", required=True, type=Path )
    pars.add_argument( "--output", help="Output masked VCF file", required=True, type=Path )


def convert_to_vcf( input_vcf: Path ) -> Path:
    # 1. Validation: Check if it's a file
    if not input_vcf.is_file():
        sys.stderr.write( f"Input file {input_vcf} does not exist or is not a file.\n." )
        sys.exit( -100 )

    # 2. Validation: Check extensions
    # Note: .vcf.gz and .bcf.gz require checking multiple suffixes
    valid_extensions = { ".vcf", ".vcf.gz", ".bcf", ".bcf.gz" }

    # Extract full extension (e.g., .vcf.gz instead of just .gz)
    ext = "".join( input_vcf.suffixes[-2:] ) if len( input_vcf.suffixes ) >= 2 else input_vcf.suffix

    if ext not in valid_extensions:
        sys.stderr.write(
            f"File extension {ext} is not supported. Only {', '.join( valid_extensions )} are supported at the moment.\n"
        )
        sys.exit( -100 )

    # 3. Return immediately if already vcf.gz
    if ext == ".vcf.gz":
        return input_vcf

    # 4. Conversion logic
    # Create a temporary file path with the correct suffix
    fd, output_vcf_path = tempfile.mkstemp( suffix=".vcf.gz" )

    # We close the file descriptor immediately as bcftools will write to the path
    os.close( fd )

    output_vcf = Path( output_vcf_path )

    try:
        subprocess.run(
            f"bcftools view -Oz -o {output_vcf} {input_vcf}",
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        # Cleanup temp file if command fails
        if output_vcf.exists():
            output_vcf.unlink()
        sys.stderr.write( f"bcftools conversion failed: {e.stderr}" )
        sys.exit( -100 )

    return output_vcf


def read_vcf( path: Path ):
    with gzip.open( path, "rt" ) as f:
        lines = [l for l in f if not l.startswith( '##' )]
    return pd.read_csv(
        io.StringIO( ''.join( lines ) ),
        dtype={
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        },
        sep='\t',
        low_memory=False,
    )


def read_gff( regions: Path ):
    gff_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gff = pd.read_csv( regions, sep="\t", header=None, comment="#", names=gff_columns )
    gff["taxa"] = gff["attribute"].str.extract( r'taxa="([^"]*)"' )
    return gff


def mask_vcf( vcf_path: Path, regions: Path, outpath: Path ):
    # Convert alignment to gzipped VCF if not already
    input_vcf = convert_to_vcf( vcf_path )

    # extract header information
    fd, header_path = tempfile.mkstemp( suffix=".txt" )
    os.close( fd )
    subprocess.run(
        f"bcftools view -h {input_vcf} > {header_path}", shell=True, check=True, capture_output=True, text=True
    )

    gff = read_gff( regions=regions )

    vcf_df = read_vcf( input_vcf )
    for _, entry in gff.iterrows():
        vcf_df.loc[vcf_df["POS"].between( entry["start"], entry["end"] ), entry["taxa"].split()] = "."

    fd, temp_vcf = tempfile.mkstemp( suffix=".vcf" )
    os.close( fd )

    with open( temp_vcf, "wt" ) as f:
        f.write( "##fileformat=VCFv4.2\n" )
        vcf_df.to_csv( f, sep="\t", index=False )

    subprocess.run(
        f"bcftools reheader -h {header_path} {temp_vcf} | bcftools view -Oz -o {outpath}",
        shell=True,
        check=True,
        capture_output=True,
        text=True
    )


def mask_entrypoint( arguments: argparse.Namespace ):
    mask_vcf(
        vcf_path=arguments.vcf,
        regions=arguments.regions,
        outpath=arguments.output,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_command_arguments( parser )
    args = parser.parse_args()
    mask_entrypoint( arguments=args )
