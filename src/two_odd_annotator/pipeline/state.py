from datetime import datetime
from pathlib import Path
import os

from bio_tools.files.fasta import write_clean_fasta_with_taxid
from bio_tools.taxa.taxonomy import map_scientific_notation_to_tax_id


from two_odd_annotator.constants import (
    RESULTS_DIR, 

    METADATA_YML, 
    CLEAN_FASTA_HEADERS_JSON, 

    DIAMOND_RESULTS, 
    HMMER_RESULTS, 
    BLASTP_RESULTS, 

    FILTERED_DIAMOND_HITS, 
    FILTERED_HMMER_HITS, 
    FILTERED_BLASTP_HITS, 

    FILTERED_DIAMOND_FASTA, 
    FILTERED_HMMER_FASTA, 
    FILTERED_BLASTP_FASTA,

    ANNOTATION_FASTA,
    ANNOTATION_MSA,
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE,
    ANNOTATION_CSV
)

from two_odd_annotator.utils.io import write_metadata


class State:
    """
    Manage the state of the two_odd_annotator pipeline. 


    This class is responsible for:
    - Validating the input (file or directory) and extracting valid FASTA files.
    - Inferring species names from input file names.
    - Tracking which pipeline steps have been completed based on the presence of expected output files.
    - Managing the output directory structure and metadata for each species.

    """
    def __init__(self, input_path: Path, output_base_dir: Path):
        self.input_path = Path(input_path)
        self.output_base_dir = Path(output_base_dir)

        self.input_fasta_paths = []  
        self._validate_fasta_input_exists()

        # store which steps have been completed for each species subdir,
        # which can be used to determine which steps need to be run / re-run
        # when resuming from an existing output directory
        self.subdir_steps_completed = {}

        # backward-compatible aggregate results dict used elsewhere in the
        # pipeline/tests
        self.results = None

        # store results of completed annotation steps at the top-level
        # output directory (one set per run, not per species)
        self.annotation_steps_completed = {
            "annotation_fasta": False,
            "annotation_msa": False,
            "annotation_msa_trim": False,
            "annotation_tree": False,
            "annotation_csv": False,
        }

        self.visualization_steps_completed = {
            "plot_1": False,
            "plot_2": False,
            "plot_3": False,
        }

        # store metadata for each species, 
        #  including inferred scientific name, tax id, original input file path, and any other relevant info
        self.metadata = None

        # Initialize output directory structure
        self._prepare_output_directory()
        self._initialize_species_subdirs()

    # ------------ HELPER METHODS ---------------
    def _get_base_name(self, path: Path) -> str:
        p = Path(path)
        name = p.name
        while "." in name:
            name = name.rsplit(".", 1)[0]
        return name
    
    def _infer_species_from_file_name(self, file_name: str) -> str:
        """Infer species name from an input FASTA file name.

        Assumes the file is named according to the Latin species name, e.g.
        "Solanum_tuberosum"-> "Solanum tuberosum".
        """

        # Replace underscores with spaces
        scientific_name = file_name.replace("_", " ").strip()

        parts = scientific_name.split()
        if len(parts) < 2:
            raise ValueError(f"Could not infer species name from file: {file_name}")

        return scientific_name
    
    def _fetch_completed_pipeline_steps_output_dir(self):
        """
        Check which annotation steps have been completed by checking the 
        output directory for the presence of expected output files. 
        This can be used to determine which steps need to be run / re-run when resuming from an existing output directory.

        """
        # check which annotation steps are completed
        if (self.output_base_dir / ANNOTATION_FASTA).exists():
            self.annotation_steps_completed["annotation_fasta"] = True
        if (self.output_base_dir / ANNOTATION_MSA).exists():
            self.annotation_steps_completed["annotation_msa"] = True
        if (self.output_base_dir / ANNOTATION_MSA_TRIM).exists():
            self.annotation_steps_completed["annotation_msa_trim"] = True
        if (self.output_base_dir / ANNOTATION_TREE).exists():
            self.annotation_steps_completed["annotation_tree"] = True
        if (self.output_base_dir / ANNOTATION_CSV).exists():
            self.annotation_steps_completed["annotation_csv"] = True   


    def _fetch_completed_pipeline_steps_for_subdir(self, subdir_path: Path) -> dict:
        """
        Check which steps have been completed based on the presence of expected output files.
        """
        pipeline_steps = {
            # can be populated with a list of the filtering methods
            # ["diamond", "hmmer", "blastp"]; annotation/visualize are
            # currently not tracked per-subdir but kept for API compatibility
            "seq_sim_filter": None,
            "annotate": None,
            "visualize": None,
        }

        expected_files_diamond = [DIAMOND_RESULTS, FILTERED_DIAMOND_HITS, FILTERED_DIAMOND_FASTA]
        expected_files_hmmer = [HMMER_RESULTS, FILTERED_HMMER_HITS, FILTERED_HMMER_FASTA]
        expected_files_blastp = [BLASTP_RESULTS, FILTERED_BLASTP_HITS, FILTERED_BLASTP_FASTA]

        # check for each filtering method if all expected files are present, if so, add to the list of completed filtering methods
        if all((subdir_path / f).exists() for f in expected_files_diamond):
            if pipeline_steps["seq_sim_filter"] is None:
                pipeline_steps["seq_sim_filter"] = []
            pipeline_steps["seq_sim_filter"].append("diamond")

        if all((subdir_path / f).exists() for f in expected_files_hmmer):
            if pipeline_steps["seq_sim_filter"] is None:
                pipeline_steps["seq_sim_filter"] = []
            pipeline_steps["seq_sim_filter"].append("hmmer")

        if all((subdir_path / f).exists() for f in expected_files_blastp):
            if pipeline_steps["seq_sim_filter"] is None:
                pipeline_steps["seq_sim_filter"] = []
            pipeline_steps["seq_sim_filter"].append("blastp")


        # append the dict of completed steps for the subdir to the overall results dict
        result = {subdir_path.name: pipeline_steps}

        return result
    
    
    # -----------INITIALIZATION METHODS ---------------
    def _validate_fasta_input_exists(self):

        # prevent input == output directory
        if self.input_path.resolve() == self.output_base_dir.resolve():
            raise ValueError(
                f"Input path and output directory must be different: "
                f"{self.input_path}"
            )

        if not self.input_path.exists():
            raise FileNotFoundError(f"{self.input_path} does not exist")

        # Single file
        if self.input_path.is_file():
            if self.input_path.suffix.lower() not in [".fasta", ".fa", ".faa"]:
                raise ValueError(f"Input file {self.input_path} is not a FASTA file")
            self.input_fasta_paths = [self.input_path]

        # Directory
        elif self.input_path.is_dir():
            fasta_files = [
                Path(root) / f
                for root, _, files in os.walk(self.input_path)
                for f in files
                if f.lower().endswith((".fasta", ".fa", ".faa")) and ".cds" not in f
                # ensure to not take any files from the output directory if it is a subdirectory of the input directory
                and not Path(root).resolve().is_relative_to(self.output_base_dir.resolve())
            ]

            if len(fasta_files) == 0:
                raise ValueError(
                    f"No FASTA files found in input directory {self.input_path}"
                )

            self.input_fasta_paths = fasta_files

            # detect duplicate species names
            species_names = [f.stem.split(".")[0] for f in fasta_files]
            duplicates = {s for s in species_names if species_names.count(s) > 1}

            if duplicates:
                raise ValueError(
                    f"Duplicate species names found: {', '.join(sorted(duplicates))}"
                )

        else:
            raise ValueError(f"Input path {self.input_path} is not a file or directory")

    def _prepare_output_directory(self):
        """
        Create the output base directory if it doesn't exist, 
        and check for any completed pipeline steps in the output directory to determine which 
        steps don't need to be re-run (if Runner.reuse_existing is set toTrue).

        The expected output directory files are:
        - annotation.fasta
        - annotation_msa.fasta
        - annotation_msa_trim.fasta
        - annotation_tree.nwk
        - annotation_results.csv
        """

        self.output_base_dir.mkdir(parents=True, exist_ok=True)

        # check if any annotation steps have been completed in the output
        # directory; this is used when reuse_existing=True to skip work.
        self._fetch_completed_pipeline_steps_output_dir()

    def _initialize_species_subdirs(self):

        for orig_file_path in self.input_fasta_paths:
            orig_file = orig_file_path.name


            print(f"Initializing subdir for {orig_file}")
            input_file_name = self._get_base_name(orig_file)
            inferred_sp_name = self._infer_species_from_file_name(input_file_name)

            # with inferred species name, get tax id and true scientific name 
            # (in case of misspellings in the filename)
            tax_id_dict = map_scientific_notation_to_tax_id(inferred_sp_name, raise_on_error=True)
            scientific_sp_name, tax_id = list(tax_id_dict.keys())[0], list(tax_id_dict.values())[0]

            subdir_folder_name = scientific_sp_name.replace(" ", "_")
            cleaned_fasta_path = self.output_base_dir / subdir_folder_name / f"clean_{subdir_folder_name}.fasta"

            # all species info is stored in metadata dictionary
            metadata_dict = {subdir_folder_name: {
                "scientific_sp_name": scientific_sp_name,
                "tax_id": tax_id,
                "original_file_path": str(orig_file_path.resolve()),
                "cleaned_fasta_path": str(cleaned_fasta_path.resolve())
            }}

            # create subdir for the species
            # first, check if subdir already exists (e.g. from a previous run)
            subdir_path = self.output_base_dir / subdir_folder_name
            if subdir_path.exists():
                # check if metadata file exists, otherwise create it
                if not (self.output_base_dir / subdir_folder_name / METADATA_YML).exists():
                    write_metadata(output_path=self.output_base_dir / subdir_folder_name , metadata=metadata_dict)

                # check if copied fasta file with cleaned headers exists, otherwise create it
                if not (self.output_base_dir / subdir_folder_name / CLEAN_FASTA_HEADERS_JSON).exists():
                    write_clean_fasta_with_taxid(
                        input_fasta_path=orig_file_path,
                        output_dir=self.output_base_dir / subdir_folder_name,
                        output_fasta_name=f"clean_{subdir_folder_name}",
                        scientific_sp_name=scientific_sp_name,
                        tax_info=tax_id,
                    )

                # check if any results files from previous runs exist, if so, keep them, otherwise initialize empty results dict
                result = self._fetch_completed_pipeline_steps_for_subdir(self.output_base_dir / subdir_folder_name)

            else:
                subdir_path.mkdir()
                write_metadata(output_path=subdir_path, metadata=metadata_dict)
                write_clean_fasta_with_taxid(
                    input_fasta_path=orig_file_path,
                    output_dir=subdir_path,
                    output_fasta_name=f"clean_{subdir_folder_name}",
                    scientific_sp_name=scientific_sp_name,
                    tax_info=tax_id,
                )

                result = {subdir_folder_name: {
                    "seq_sim_filter": None,
                    "annotate": None,
                    "visualize": None
                }}

            # append the result for the subdir to the overall results dict
            if self.results is None:
                self.results = result
            else:
                self.results.update(result)

            # append the metadata for the subdir to the overall metadata dict
            if self.metadata is None:
                self.metadata = metadata_dict
            else:
                self.metadata.update(metadata_dict)


