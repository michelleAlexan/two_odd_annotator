## 2ODD Annotator

The 2ODD Annotator is a bioinformatics helper package for identifying putative 2‑oxoglutarate/Fe(II)‑dependent dioxygenases (2ODD) in plant proteomes based on sequence similarity and phylogenetic grouping.


## Overview
The `two-odd-annotator` is a Python pipeline designed to identify and annotate 2ODD enzymes from plant protein FASTA sequences. The pipeline is configurable via a YAML file, and can be run as a Python library or command-line tool. It provides three services:

1. **Sequence‑similarity filtering** of input plant proteomes for gene sequences with significant 2ODD similarity. There are three methods for this:
        - DIAMOND-based filtering utilizing the [diamond reference database](data/char_2ODDs/dmnd_ref_db.dmnd) of characterized 2ODD sequences
                - This is the default method for sequence similarity filtering, as it is much faster than BLASTP
        - HMMER‑based filtering using the [2ODD domain HMM profile](data/2ODD_domain.hmm) constructed from characterized 2ODD sequences
        - BLASTP‑based filtering using the [BLAST reference database](data/blast_ref_db/2ODD_ref_db) of characterized 2ODD sequences. 


2. **Phylogenetic grouping / annotation** (planned): group candidate 2ODDs into functional clades based on phylogenetic trees. 

3. **Visualization** (planned): generate summary plots of the distribution of 2ODD candidates across species and clades.


---

## Installation

Requirements:

- Python ≥ 3.13
- A Unix‑like environment (Linux or macOS) where either DIAMOND, HMMER or BLAST+ are installed and available on the command line:
    - `diamond` (from DIAMOND) is used in `services.seq_sim_filter.run_diamond`.
    - `hmmsearch` (from HMMER 3) is used in `services.seq_sim_filter.run_hmmsearch`.        
    - `blastp` (from BLAST+) is used in `services.seq_sim_filter.run_blastp`.


Make sure these binaries are available on your `PATH`, e.g.:
```
which diamond
which blastp
which hmmsearch
```

Python dependencies are managed via the project’s `pyproject.toml`. At minimum you will need:

- `pandas`
- `pyyaml`
- `biopython` (for `Bio.SeqIO`)


### Using uv (recommended for development)

```
git clone https://github.com/michelleAlexan/two_odd_annotator.git
cd two_odd_annotator
uv sync
```

### Using pip

In a fresh virtual environment:

```
git clone https://github.com/michelleAlexan/two_odd_annotator.git
cd two_odd_annotator
pip install -e .
```


---

## Configuration

You can configure the pipeline parameters, reference databases and thresholds via a YAML file. By default, the pipeline looks for `configs/filter.yml` in the project directory, but you can specify a different config file path when running the pipeline.


Example config.yml:
```
pipeline:
  reuse_existing: true # whether to reuse existing intermediate files if they exist, 
                              # set to false to force re-running all steps of the pipeline

  sp_name_mapping: "configs/sp_name_correction.json" # path to JSON file mapping species names in input FASTA headers to corrected names.
                                                  # this is needed, e.g.,  when the latin species name does not map to a taxonomic ID in the NCBI taxonomy database

  seq_sim_method: "diamond"  # default method for sequence similarity filtering, 
                              # choices are "diamond", "hmmer", "blastp" 
  compute_plots: false # whether to compute summary plots at the end of the pipeline run



filter_tools:
  blastp:
    reference_db: "data/char_2ODDs/blast_ref_db/2ODD_ref_db"

  diamond:
    reference_db: "data/char_2ODDs/dmnd_ref_db"  

  hmmer:
    domain_model: "data/char_2ODDs/2ODD_domain.hmm"


parameters:
  threads: 8
  thresholds_alignment:
    evalue: 1e-5
    pident: 15.0
    length: 80
    bitscore: 40
    num_hits: 100
  
  thresholds_hmmer:
    full_Evalue: 1e-5
    bestdom_Evalue: 1e-5
    full_score: 50
    bestdom_score: 50
    N: 1
  
```

Optionally, you can override specific config parameters via command-line arguments when running the pipeline, e.g.:
```
annodd --input-path data/fasta \
       --output-dir results \
       --config-path configs/filter.yml \
       --reuse-existing false \
       --seq-sim-method hmmer \
       --compute-plots true
```

or when using the Python library:
```python
from two_odd_annotator.pipeline.runner import Runner   

runner = Runner(
    input_path = "data/fasta", 
    output_dir = "results", 
    config_path = "configs/filter.yml")

runner.run(
        reuse_existing = False,
        seq_sim_method = "hmmer",
        compute_plots = True)
```     


---

## Result folder structure

For each species, a subdirectory will be created containing filtered results:

```
results/
├── Arabidopsis_thaliana/
│   ├── metadata.yml
│   ├── clean_fasta_headers.json
│   ├── diamond_results.tsv / blastp_results.tsv / hmmer_results.out
│   ├── filtered_diamond_hits.csv / filtered_blastp_hits.csv / filtered_hmmer_hits.csv
│   └── filtered_diamond.fasta / filtered_blastp.fasta / filtered_hmmer.fasta
├── Oryza_sativa/
│   └── ...
└── filter.log
```
---




