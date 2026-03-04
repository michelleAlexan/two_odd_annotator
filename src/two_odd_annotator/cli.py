import argparse
from two_odd_annotator.pipeline import Runner
from two_odd_annotator.constants import DEFAULT_CONFIG_PATH

from importlib.metadata import version


def main():

    parser = argparse.ArgumentParser(
        prog="annodd", 
        description="Run 2ODD annotation pipeline for plant protein sequences."
        )
    
    parser.add_argument(
        "--version", 
        action="version", 
        version=f"annodd {version('two_odd_annotator')}"
        )
    
    parser.add_argument(
        "--input-path", 
        required=True,
        help="Path to the input FASTA file or directory containing FASTA files. " \
        "FASTA files should be named according to the Latin species name," \
        " e.g. 'Solanum_tuberosum.fasta' or 'Solanum_tuberosum.pep.fasta'."
    )
    
    parser.add_argument(
        "--output-dir", 
        required=True, 
        help="Path to the output base directory where results will be stored. " \
        "Each species will have its own subdirectory named according to the Latin species name, e.g. 'Solanum_tuberosum'."
    )

    parser.add_argument(
        "--config-path", 
        default=DEFAULT_CONFIG_PATH,
        help="Path to the config YAML file defining parameters for the pipeline run."
    )

    parser.add_argument(
        "--reuse-existing",
        choices=["true","false"],
        help="Override config: whether to reuse existing results or force rerun."
    )

    parser.add_argument(
        "--sp-name-mapping",
        help="Override config: path to json file mapping incorrect species names to correct ones." \
        "This is useful when mapping species names to the NCBI taxonomy id fails."
    )

    parser.add_argument(
    "--seq-sim-method",
    choices=["diamond","blastp","hmmer"],
    help="Override config: method for sequence similarity filtering."
    )

    parser.add_argument(
        "--compute-plots",
        choices=["true","false"],
        help="Override config: whether to compute summary plots after pipeline run."
    )

    args = parser.parse_args()


    # Initialize pipeline
    pipeline = Runner(
        input_path = args.input_path, 
        output_dir = args.output_dir, 
        config_path = args.config_path
        )

    
    pipeline.run(
        reuse_existing=args.reuse_existing.lower() == "true" if args.reuse_existing else None,
        sp_name_mapping=args.sp_name_mapping,
        seq_sim_method=args.seq_sim_method,
        compute_plots=args.compute_plots.lower() == "true" if args.compute_plots else None
    )
