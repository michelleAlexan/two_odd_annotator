import argparse
from two_odd_annotator.pipeline import Runner

from importlib.metadata import version


def main():
    """
    Command line interface for running the 2ODD annotation pipeline.
    """
    parser = argparse.ArgumentParser(
        prog="annodd", 
        description="Run 2ODD annotation pipeline for plant protein sequences."
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
        help="Path to custom config YAML file. If not provided, default config is used."
    )

    parser.add_argument(
        "--version", 
        action="version", 
        version=f"annodd {version('two_odd_annotator')}"
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
    help="Override config: method for sequence similarity filtering. Default is 'hmmer'.")


    parser.add_argument(
        "--compute-plots",
        choices=["true","false"],
        help="Override config: whether to compute summary plots after pipeline run."
    )

    parser.add_argument(
        "--seq-len-thresh",
        help="Override config: the length of the candidate sequence must be within +/- this value of the median 2ODD sequence length of 349 aa. " \
        "Default is 100." \
        "If set to -1, no length filtering will be applied."
    )

    parser.add_argument(
        "--delete-intermediate-files",
        choices=["true","false"],
        help="Override config (default: false): whether to delete intermediate files (e.g. alignment files) after annotation is complete."
    )

    args = parser.parse_args()



    # Initialize pipeline
    pipeline = Runner(
        input_path=args.input_path,
        output_base_dir=args.output_dir,
        config_path=args.config_path,
        reuse_existing=args.reuse_existing.lower() == "true" if args.reuse_existing else None,
        sp_name_mapping=args.sp_name_mapping,
        seq_sim_method=args.seq_sim_method,
        compute_plots=args.compute_plots.lower() == "true" if args.compute_plots else None,
        seq_len_thresh=args.seq_len_thresh,
        delete_intermediate_files=args.delete_intermediate_files.lower() == "true" if args.delete_intermediate_files else None,
    )

    pipeline.run()

