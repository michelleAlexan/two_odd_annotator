from pathlib import Path
import timeit

from two_odd_annotator.constants import DEFAULT_CONFIG_PATH
from two_odd_annotator.utils import io
from two_odd_annotator.utils.logging import init_log, log_line
from two_odd_annotator.pipeline import state
from two_odd_annotator.services import seq_sim_filter, annotate, vizualize

from two_odd_annotator.constants import (
    ANNOTATION_FASTA,
    ANNOTATION_MSA,
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE,
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
        threads: int | None = None,
        reuse_existing: bool | None = None,
        sp_name_mapping: str | None = None,
        seq_sim_method: str | None = None,
        compute_plots: bool | None = None,
        seq_len_thresh: str | None = None,
        delete_intermediate_files: bool | None = None,
        step: str | None = None,
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
        if threads is not None:
            if threads < 1:
                raise ValueError(f"threads must be >= 1, got {threads}")
            self.config.setdefault("parameters", {})
            self.config["parameters"]["threads"] = int(threads)
        if reuse_existing is not None:
            self.config["pipeline"]["reuse_existing"] = reuse_existing
        if sp_name_mapping is not None:
            self.config["pipeline"]["species_name_map"] = sp_name_mapping
        if seq_sim_method is not None:
            self.config["pipeline"]["seq_sim_method"] = seq_sim_method
        if compute_plots is not None:
            self.config["pipeline"]["compute_plots"] = compute_plots
        if seq_len_thresh is not None:
            self.config["pipeline"]["seq_len_thresh"] = int(seq_len_thresh)
        if delete_intermediate_files is not None:
            self.config["pipeline"][
                "delete_intermediate_files"
            ] = delete_intermediate_files

        # which part of the pipeline to run
        # "all" (default) runs the full pipeline; otherwise one of
        # "filter_seq_sim", "annotate", or "visualize"
        self.step = step or "all"

        valid_steps = {"all", "filter_seq_sim", "annotate", "visualize"}
        if self.step not in valid_steps:
            raise ValueError(
                f"Invalid step: {self.step!r}. Expected one of {sorted(valid_steps)}."
            )

        # initialize run log in the output base directory
        self.log_path = init_log(str(self.output_base_dir))

        # log basic run context and effective configuration
        log_line(self.log_path, "====== INITIALIZING PIPELINE ======")
        log_line(self.log_path, f"Input path: {self.input_path}")
        log_line(self.log_path, f"Output dir: {self.output_base_dir}")
        log_line(self.log_path, f"Config path: {self.config_path}")

        # log a compact representation of the effective config
        from pprint import pformat

        log_line(self.log_path, "Effective config:")
        for line in pformat(self.config).splitlines():
            log_line(self.log_path, f"  {line}")

        # resolve any relative paths in the config to be absolute
        # based on the project root so the pipeline can be run from
        # any current working directory.
        self._resolve_config_paths()

        # initialize the pipeline run by creating a state object.
        # State tracks initializes output directory and subdirectories,
        # validates input files,
        # and tracks the progress of each step in the pipeline for each input file.
        self.state = state.State(
            input_path=self.input_path,
            output_base_dir=self.output_base_dir,
            log_path=self.log_path,
        )

    def _resolve_config_paths(self) -> None:
        """Normalize path-like entries in the config to absolute paths.

        The YAML config stores paths like "data/2ODDs/dmnd_ref_db"; when the
        pipeline is executed from outside the project root, those relative
        paths would otherwise break. Here we resolve them relative to the
        repository root (three parents above this file).
        """

        project_root = Path(__file__).parents[3]

        def _abs(p: str | None) -> str | None:
            if p is None:
                return None
            path = Path(p)
            if path.is_absolute():
                return str(path)
            return str(project_root / path)

        pipeline_cfg = self.config.get("pipeline", {})
        for key in ("sp_name_mapping", "species_name_map"):
            if key in pipeline_cfg:
                pipeline_cfg[key] = _abs(pipeline_cfg[key])

        filter_tools = self.config.get("filter_tools", {})
        if "blastp" in filter_tools and "reference_db" in filter_tools["blastp"]:
            filter_tools["blastp"]["reference_db"] = _abs(
                filter_tools["blastp"]["reference_db"]
            )
        if "diamond" in filter_tools and "reference_db" in filter_tools["diamond"]:
            filter_tools["diamond"]["reference_db"] = _abs(
                filter_tools["diamond"]["reference_db"]
            )
        if "hmmer" in filter_tools and "domain_model" in filter_tools["hmmer"]:
            filter_tools["hmmer"]["domain_model"] = _abs(
                filter_tools["hmmer"]["domain_model"]
            )

        annotate_cfg = self.config.get("annotate", {})
        if "bait_sequence_collection" in annotate_cfg:
            annotate_cfg["bait_sequence_collection"] = _abs(
                annotate_cfg["bait_sequence_collection"]
            )
        if "major_2ODDs_functional_characterization" in annotate_cfg:
            annotate_cfg["major_2ODDs_functional_characterization"] = _abs(
                annotate_cfg["major_2ODDs_functional_characterization"]
            )
        if "major_minor_2ODD_ids" in annotate_cfg:
            annotate_cfg["major_minor_2ODD_ids"] = _abs(
                annotate_cfg["major_minor_2ODD_ids"]
            )

    def run(self):
        """Run the full two_odd_annotator pipeline.
        Steps:
        1. Sequence similarity filtering (HMMER, Diamond, or BLASTP)
        2. 2ODD annotation (based on phylogenetic evidence)
        3. Visualization of results"""

        start = timeit.default_timer()

        # 1) sequence similarity filtering
        if self.step in ("all", "filter_seq_sim"):
            log_line(self.log_path, "====== SEQUENCE SIMILARITY FILTER ======")
            for subdir, metadata in self.state.metadata.items():
                subdir_path = self.state.output_base_dir / subdir
                cleaned_fasta_path = metadata["cleaned_fasta_path"]

                method = self.config["pipeline"]["seq_sim_method"]

                # special mode: run all three methods sequentially
                if method == "all":
                    for m in ("diamond", "hmmer", "blastp"):
                        log_line(self.log_path, f"[SEQ_SIM] {m} for {subdir_path.name}")
                        seq_sim_filter.run(
                            input_path=cleaned_fasta_path,
                            subdir=subdir_path,
                            config=self.config,
                            seq_sim_method=m,
                        )
                else:
                    log_line(
                        self.log_path, f"[SEQ_SIM] {method} for {subdir_path.name}"
                    )
                    seq_sim_filter.run(
                        input_path=cleaned_fasta_path,
                        subdir=subdir_path,
                        config=self.config,
                        seq_sim_method=method,
                    )

        # 2) phylogenetic annotation
        if self.step in ("all", "annotate"):
            log_line(self.log_path, "====== ANNOTATION ======")
            if (
                self.config["pipeline"]["reuse_existing"]
                and self.state.annotation_steps_completed["annotation_csv"]
            ):
                print("Reusing existing annotation results...")
                log_line(
                    self.log_path,
                    "Reusing existing annotation results (annotation_csv present).",
                )
            else:
                log_line(self.log_path, "Running annotation service.")
                annotate.run(
                    result_dir=self.state.output_base_dir,
                    config=self.config,
                    seq_sim_method=self.config["pipeline"]["seq_sim_method"],
                    completed_annotation_steps=self.state.annotation_steps_completed,
                )

        # 3) visualization (optional / placeholder for now)
        if self.step in ("all", "visualize"):
            log_line(self.log_path, "====== VISUALIZATION ======")
            if self.config["pipeline"].get("compute_plots"):
                if hasattr(vizualize, "run"):
                    log_line(self.log_path, "Running visualization service.")
                    vizualize.run(self.state.output_base_dir, self.config)
                else:
                    log_line(
                        self.log_path, "Visualization step is not implemented yet."
                    )

        elapsed = (timeit.default_timer() - start) / 60

        # only delete intermediate annotation files if we actually ran
        # the annotation step in this invocation
        if self.step in ("all", "annotate") and self.config["pipeline"].get(
            "delete_intermediate_files"
        ):
            log_line(
                self.log_path,
                f"Deleting intermediate annotation files: {ANNOTATION_FASTA, ANNOTATION_MSA, ANNOTATION_MSA_TRIM, ANNOTATION_TREE}...",
            )
            # delete intermediate files (ANNOTATION_FASTA, ANNOTATION_MSA, ANNOTATION_MSA_TRIM, ANNOTATION_TREE)
            # that are found in the output directory (results folder)
            for file in self.state.output_base_dir.rglob("*"):
                if file.name in [
                    ANNOTATION_FASTA,
                    ANNOTATION_MSA,
                    ANNOTATION_MSA_TRIM,
                    ANNOTATION_TREE,
                ]:
                    file.unlink()

        log_line(self.log_path, f"Pipeline completed in {elapsed:.2f} minutes.")
