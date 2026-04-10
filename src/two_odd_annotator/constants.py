# Default config file path to use if not provided by user
DEFAULT_CONFIG_PATH = "configs/default_config.yml"

# Name of the output base directory for all results
RESULTS_DIR = "2ODD_annotation"

# output files for each species subdirectory
METADATA_YML = "metadata.yml"
CLEAN_FASTA_HEADERS_JSON = "clean_fasta_headers.json"

DIAMOND_RESULTS = "diamond_results.tsv"
HMMER_RESULTS = "hmmer_results.out"
BLASTP_RESULTS = "blastp_results.tsv"

FILTERED_DIAMOND_HITS = "filtered_diamond_hits.csv"
FILTERED_HMMER_HITS = "filtered_hmmer_hits.csv"
FILTERED_BLASTP_HITS = "filtered_blastp_hits.csv"

FILTERED_DIAMOND_FASTA = "filtered_diamond.fasta"
FILTERED_HMMER_FASTA = "filtered_hmmer.fasta"
FILTERED_BLASTP_FASTA = "filtered_blastp.fasta"

# combines all filtered FASTA files and ingroup sequences for annotation step
ANNOTATION_FASTA = "annotation.fasta" 
ANNOTATION_MSA = "annotation_msa.fasta"
ANNOTATION_MSA_TRIM = "annotation_msa_trim.fasta"
ANNOTATION_TREE = "annotation_tree.nwk"
ANNOTATION_CSV = "annotation_results.csv"
