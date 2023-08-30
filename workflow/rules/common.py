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
