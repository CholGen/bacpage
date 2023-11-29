import argparse
from io import StringIO
from subprocess import run

from Bio import SeqIO


def convert_vcf( vcf: str, reference: str, output: str ):
    # load reference and get id
    reference_seq = SeqIO.read( reference, "fasta" )
    reference_id = reference_seq.id

    # Load samples in vcf and remove reference if present.
    query = run( f"bcftools query -l {vcf}", shell=True, capture_output=True, text=True )
    samples = query.stdout.strip().split( "\n" )
    samples = [sample for sample in samples if sample != reference_id]

    records = list()
    for sample in samples:
        consensus = run(
            f"bcftools consensus -f {reference} -s '{sample}' {vcf}", shell=True, capture_output=True, text=True
        )
        try:
            record = SeqIO.read( StringIO( consensus.stdout ), "fasta" )
        except ValueError:
            print( consensus )
            raise
        record.id = sample
        record.name = ""
        record.description = ""
        records.append( record )

    SeqIO.write( records, output, "fasta" )


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
