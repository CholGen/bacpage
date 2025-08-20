import argparse
import tempfile
from io import StringIO
from subprocess import run

from Bio import SeqIO


def convert_vcf( vcf: str, reference: str, output: str ):
    with tempfile.NamedTemporaryFile( delete=False, suffix=".fasta" ) as temp_file:
        temp_file_path = temp_file.name

    # load reference and get id
    reference_seq = SeqIO.read( reference, "fasta" )
    reference_id = reference_seq.id
    run( f"union -filter {reference} > {temp_file_path}", shell=True )

    # Load samples in vcf and remove reference if present.
    query = run( f"bcftools query -l {vcf}", shell=True, capture_output=True, text=True )
    samples = query.stdout.strip().split( "\n" )
    samples = [sample for sample in samples if sample != reference_id]

    records = list()
    with open( output, "w" ) as output_file:
        for sample in samples:
            consensus = run(
                f"bcftools consensus --mark-del N -f {temp_file_path} -s '{sample}' {vcf}", shell=True,
                capture_output=True,
                text=True
            )
            record = SeqIO.read( StringIO( consensus.stdout ), "fasta" )
            record.id = sample
            record.name = ""
            record.description = ""
            record.seq = record.seq.replace( "*", "N" ).upper()
            SeqIO.write( record, output_file, "fasta" )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts a multi-sample VCF file to a multisequence alignment in FASTA format."
    )

    # Initialize optional arguments
    parser.add_argument( "--vcf", type=str, help="VCF file to convert" )
    parser.add_argument( "--reference", type=str, help="Reference sequence used when generating the VCF." )
    parser.add_argument( "--output", type=str, help="Location to save fasta file." )
    args = parser.parse_args()

    convert_vcf( args.vcf, args.reference, args.output )
