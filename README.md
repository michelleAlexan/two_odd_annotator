## 2ODD Annotator

The 2ODD Annotator is a bioinformatics helper package for identifying putative 2‑oxoglutarate/Fe(II)‑dependent dioxygenases (2ODDs) in plant proteomes based on sequence similarity and phylogenetic grouping.

It provides two services:

1. **Sequence‑similarity filtering** of input plant proteomes for gene sequences with significant 2ODD similarity. There are two approaches for this:
        - BLASTP‑based filtering (`filter_seq_sim.blast`)
        - HMMER‑based filtering using an HMM profile (`filter_seq_sim.hmmer`)

        Both filters work on protein FASTA files and can process either a single file or recursively walk a directory of FASTA files.

2. **Phylogenetic grouping / annotation** (planned): group candidate 2ODDs into functional clades based on phylogenetic trees. 

---

## Installation

Requirements:

- Python ≥ 3.13
- A Unix‑like environment (Linux or macOS) where BLAST+ and HMMER are installed and available on the command line:
    - `blastp` (from BLAST+) is used in `filter_seq_sim.blast.run_blastp`.
    - `hmmsearch` (from HMMER 3) is used in `filter_seq_sim.hmmer.run_hmmsearch`.

Make sure these binaries are available on your `PATH`, e.g.:
```bash
which blastp
which hmmsearch
```

`blastp` relies on a BLAST protein database as a reference to characterized 2ODD sequences (provided in `data/blast_ref_db/`).
`hmmsearch` relies on an HMM profile constructed from characterized 2ODD sequences (provided as `data/2ODD_domain.hmm`).

Python dependencies are managed via the project’s `pyproject.toml`. At minimum you will need:

- `pandas`
- `pyyaml`
- `biopython` (for `Bio.SeqIO`)


### Using uv (recommended for development)

```
git clone https://github.com/michelleAlexan/2ODD_annotator.git
cd 2ODD_annotator
uv sync
```

### Using pip

In a fresh virtual environment:

```
git clone https://github.com/michelleAlexan/2ODD_annotator.git
cd 2ODD_annotator
pip install -e .
```


---

## Configuration

Both BLAST and HMMER workflows are driven by a YAML configuration file, for example `configs/filter_seq_sim.yml`.

### BLAST configuration

The BLAST‑based filter expects (at minimum) keys like:

```yaml
input: path/to/proteomes_or_fasta
output_dir: results/filter_seq_sim
reuse_existing: true

tools:
        blast:
                reference_db: data/blast_ref_db/2ODD_ref_db
                threads: 8
                thresholds:
                        evalue: 1e-5
                        pident: 40
                        length: 100
                        bitscore: 50
                        num_hits: 5
```

### HMMER configuration

The HMMER‑based filter expects (at minimum):

```yaml
input: path/to/proteomes_or_fasta
output_dir: results/filter_seq_sim_hmmer
reuse_existing: true

tools:
        hmmer:
                domain_model: data/2ODD_domain.hmm
```

You can reuse the same top‑level keys (`input`, `output_dir`, `reuse_existing`) across both tools.


---

## Usage

The package is currently used as a Python library rather than a finished command‑line tool.  
Typical usage is to load your YAML configuration and call either the BLAST or HMMER pipeline.

### BLASTP‑based sequence similarity filter

```python
from 2odd_annotator.filter_seq_sim import blast

config = blast.load_config("configs/filter_seq_sim.yml")

input_path = config["input"]
if os.path.isfile(input_path):
        blast.main(input_path, config)
else:
        blast.process_directory(input_path, config)
```

This will:

- Run `blastp` against the specified reference database.
- Filter hits by identity, alignment length, bitscore, and number of hits per query.
- Write:
    - `blast_results.tsv` – raw BLAST tabular output.
    - `blast_filtered_hits.csv` – filtered hits table.
    - `filter_seq_sim_blast.fasta` – FASTA sequences that passed the thresholds.

### HMMER‑based sequence similarity filter

```python
from 2odd_annotator.filter_seq_sim import hmmer

config = hmmer.load_config("configs/filter_seq_sim.yml")

input_path = config["input"]
if os.path.isfile(input_path):
        hmmer.main(input_path, config)
else:
        hmmer.process_directory(input_path, config)
```

This will:

- Run `hmmsearch` with your 2ODD HMM profile.
- Parse the "Scores for complete sequences" table into a pandas DataFrame.
- Optionally restrict to sequences with a single domain (`N == 1`).
- Write:
    - `hmmer2ODD_hits.out` – raw hmmsearch text output.
    - `hmmer2ODD_hits.csv` – parsed table of hits.
    - `filter_seq_sim_hmmer.fasta` – FASTA of sequences with accepted HMMER hits.

When `reuse_existing` is set to `true` in the config, both pipelines will skip recomputation for FASTA files where all expected output files already exist.

---




