import argparse
import pandas as pd
from subprocess import run
from io import StringIO

IDX_COLUMNS = ["gene", "length", "reads_mapped", "reads_unmapped"]
COV_COLUMNS = ["gene", "position", "depth"]


def collect_reads( alignment: str ) -> pd.DataFrame:
    command = run( f"samtools idxstats {alignment}", shell=True, capture_output=True, text=True )
    reads = pd.read_csv( StringIO( command.stdout ), sep="\t", names=IDX_COLUMNS )
    reads["frac_mapped"] = reads["reads_mapped"] / (reads["reads_mapped"] + reads["reads_unmapped"]).sum()
    reads = reads.loc[reads["gene"] != "*"]
    return reads


def calculate_coverage( alignment: str, minimum_coverage: int ) -> pd.DataFrame:
    command = run( f"samtools depth -a {alignment}", shell=True, capture_output=True, text=True )
    cov = pd.read_csv( StringIO( command.stdout ), sep="\t", names=COV_COLUMNS )
    cov = cov.groupby( "gene" )["depth"].agg( ["count", "mean", lambda x: sum( x > minimum_coverage )] )
    cov.columns = ["pos_total", "mean_depth", "pos_covered"]
    cov["frac_covered"] = cov["pos_covered"] / cov["pos_total"]
    return cov


def calculate_stats( alignment: str, name: str, minimum_coverage: int, output: str ) -> None:
    reads_df = collect_reads( alignment )
    coverage = calculate_coverage( alignment, minimum_coverage=minimum_coverage )
    output_df = coverage.merge( reads_df, on="gene" )
    output_df = output_df.drop( columns=["reads_unmapped", "length"] )
    output_df["sample"] = name
    output_df = output_df[
        ["sample", "gene", "mean_depth", "pos_total", "pos_covered", "frac_covered", "reads_mapped", "frac_mapped"]]
    output_df.to_csv( output, index=False )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Calculates alignment stats for a given BAM file." )

    # Initialize optional arguments
    parser.add_argument( "-a", "--alignment", required=True, type=str, help="input alignment" )
    parser.add_argument( "-m", "--min-coverage", required=True, type=int,
                         help="minimum coverage required to call a base." )
    parser.add_argument( "-s", "--sample-name", required=False, default=None, help="" )
    parser.add_argument( "-o", "--output", required=True, type=str, help="location to save output" )

    args = parser.parse_args()
    sample_name = args.sample_name if args.sample_name else args.alignment
    calculate_stats( alignment=args.alignment, name=sample_name, minimum_coverage=args.min_coverage,
                     output=args.output )
