# Functions that are common to either Illumina or ONT pipelines.
import os


def normalize_path( path: str, working_directory: str ) -> str:
    if not os.path.isabs( path ):
        return os.path.join( working_directory, path )
    else:
        return path


def get_genes( search_directory: str ) -> dict[str, str]:
    genes = dict()
    for file in os.listdir( search_directory ):
        if file.endswith( (".fa", ".fasta") ):
            gene = os.path.splitext( file )[0]
            path = os.path.join( search_directory, file )
            genes[gene] = path
    return genes


def determine_outputs( samples: dict[str, str], config: dict[str, str] ) -> list[str]:
    outputs = list()
    if config["generate"]["consensus_sequences"]:
        for key in samples.keys():
            outputs.append( f"results/consensus/{key}.masked.fasta" )
            outputs.append( f"results/consensus/{key}.consensus.fasta" )
    if config["generate"]["typing"]:
        outputs.append( "results/reports/typing_information.csv" )
        #outputs.append( "results/reports/mlst_types.csv" )
        outputs.append( "results/reports/antibiotic_resistance.tsv" )
    if config["generate"]["quality_control_report"]:
        outputs.append( "results/reports/qc_report.html" )
    if config["generate"]["phylogeny"]:
        outputs.append( "results/phylogeny/sparse_alignment.fasta" )
    return outputs
