from pathlib import Path
import timeit

from two_odd_annotator.constants import DEFAULT_CONFIG_PATH
from two_odd_annotator.utils import io
from two_odd_annotator.pipeline import state 
from two_odd_annotator.services import (
    seq_sim_filter, 
    annotate, 
    vizualize
)


class Runner:
    """Run the two_odd_annotator pipeline. 
    Orchestrate the execution of each step in the pipeline, passing the necessary inputs and outputs between steps.
    Handle logging and error handling for the entire pipeline run.
    """

    def __init__(
            self, 
            # required arguments
            input_path: str, 
            output_base_dir: str, 

            # optional arguments
            config_path: str | None = None,
            reuse_existing: bool | None = None,
            sp_name_mapping: str | None = None,
            seq_sim_method: str | None = None,
            compute_plots: bool | None = None
        ):
        self.input_path = Path(input_path)
        self.output_base_dir = Path(output_base_dir)

        # If no config path provided, use default config path. 
        # If config path provided, load config from that path.
        if config_path is None:
            config_path = Path(__file__).parents[3] / DEFAULT_CONFIG_PATH
        self.config_path = config_path
        self.config = io.load_config(self.config_path)

        # override config values with any provided optional arguments
        if reuse_existing is not None:
            self.config["pipeline"]["reuse_existing"] = reuse_existing
        if sp_name_mapping is not None:
            self.config["pipeline"]["species_name_map"] = sp_name_mapping
        if seq_sim_method is not None:
            self.config["pipeline"]["seq_sim_method"] = seq_sim_method
        if compute_plots is not None:
            self.config["pipeline"]["compute_plots"] = compute_plots

        # initialize the pipeline run by creating a state object.
        # State tracks initializes output directory and subdirectories,
        # validates input files, 
        # and tracks the progress of each step in the pipeline for each input file.
        self.state = state.State(
            input_path=self.input_path,
            output_base_dir=self.output_base_dir
        )


    def run(self):
        """Run the full two_odd_annotator pipeline. 
        Steps:
        1. Sequence similarity filtering (Diamond, HMMER or BLASTP)
        2. 2ODD annotation (based on phylogenetic evidence)
        3. Visualization of results"""

        start = timeit.default_timer()

        for subdir, metadata in self.state.metadata.items():
            subdir = self.state.output_base_dir / subdir
            cleaned_fasta_path = metadata["cleaned_fasta_path"]

            seq_sim_filter.run(
                input_path = cleaned_fasta_path, 
                subdir = subdir, 
                config = self.config, 
                seq_sim_method = self.config["pipeline"]["seq_sim_method"]
                )
            
            # annotate.run(cleaned_fasta_path, subdir, self.config_path)
        
        # if self.config["pipeline"]["compute_plots"]:
        #     vizualize.run(self.state.output_base_dir, self.config_path)


        elapsed = (timeit.default_timer() - start) / 60
        print(f"Pipeline completed in {elapsed:.2f} minutes.")



