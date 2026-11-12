# Reproducing the analysis

Tested on Ubuntu 24.04 (WSL2). Heavy steps are local scripts with checkpoints, so any step can be stopped and resumed by repeating the same command.

## 1. Environments

```bash
conda env create -f environment/environment.yml          # analysis environment (pinned versions)
conda create -n phanotate --file environment/phanotate-explicit.txt   # PHANOTATE 1.6.7
conda activate phageannoscope
python -m venv --system-site-packages .venv && .venv/bin/pip install --no-deps -e .
```

`configs/run_linux.json` contains the absolute path of the PHANOTATE executable (`phanotate_command`). Edit it to match your installation. Then run the tests:

```bash
PHAGEANNOSCOPE_TEST_CONFIG=configs/run_linux.json .venv/bin/python -m pytest -q
```

Pinned versions: Python 3.13.16, Biopython 1.88, pyarrow 25.0.0, duckdb 1.5.6, Pyrodigal-gv 0.3.2 / Pyrodigal 3.7.1, PHANOTATE 1.6.7, pyhmmer 0.12.3, statsmodels 0.15.0, Phold 1.2.6 / Foldseek 10.941cd33 (H3 only).

## 2. Cohort

1. Download the GenBank records listed in `manifests/accessions_full_cohort.txt` from NCBI (`phageannoscope fetch-ncbi`, or E-utilities with `efetch -format gb`).
   `scripts/pipeline/fetch_genbank_records.py` writes one file per accession in the layout of the validated manifest; the files are byte-identical to those of the study unless NCBI has updated a record. `scripts/pipeline/make_batch_manifests.py` then builds `manifests/batch1-4.csv` and reports any record whose file changed (its sequence checksum is still enforced).
2. Run validation and the protocol amendments:

```bash
.venv/bin/phageannoscope validate-cohort --manifest <frozen INPHARED2 manifest> --genbank-dir data/raw/ncbi/primary --out <validated.csv>
.venv/bin/python scripts/cohort/resolve_host_genera.py --manifest <validated.csv> --out data/raw/ncbi/host_taxonomy
.venv/bin/python scripts/cohort/apply_amendments_a1_a2.py --manifest <validated.csv> --hosts data/raw/ncbi/host_taxonomy/metadata_host_domains.csv --out <validated_amended.csv>
.venv/bin/python scripts/cohort/apply_amendment_a5.py --manifest <validated_amended.csv> --out <cohort_manifest_validated.csv> --samples <sample manifests>
```

The frozen INPHARED2 metadata is in `manifests/inphared2_metadata_frozen_2026-10-16.csv.gz`. The final validated manifest is in `manifests/cohort_manifest_validated.csv.gz`.

## 3. PHROGs HMM library (about 40 min, single thread)

```bash
.venv/bin/python scripts/phrogs/build_hmm_library.py --msa MSA_phrogs.tar.gz --out data/db/phrogs_v4/phrogs_v4_m50.h3m
```

Check the inputs against `manifests/phrogs_v4_SHA256SUMS.txt`.

## 4. Gene calling, matching and PHROGs search per batch

Build a run manifest for each batch by selecting the rows of the validated manifest whose accessions are listed in `manifests/accessions_batch*.txt`. Then:

```bash
scripts/pipeline/run_batches.sh manifests/batch1.csv results/batch1 manifests/batch2.csv results/batch2 \
                                manifests/batch3.csv results/batch3 manifests/batch4.csv results/batch4
```

Set `WORKERS` and `CORES` to suit your machine; the script header shows the settings for a 16-core workstation. A batch that was interrupted on another machine can be assembled from its checkpoints with `scripts/pipeline/rebuild_from_checkpoints.py`.

## 5. Merge, statistics and figures

```bash
.venv/bin/python scripts/pipeline/merge_runs.py --out results/full_cohort --runs results/batch1 results/batch2 results/batch3 results/batch4
REUSE= scripts/pipeline/post_run.sh results/full_cohort figures
```

For cohorts as large as this one, run `analyze_run.py` once to build `analysis/locus_table.parquet`. Later runs can then use `--reuse-locus-table` (set `REUSE=1` for `post_run.sh`), which needs much less memory.

## 6. Phold (H3)

```bash
.venv/bin/python scripts/analysis/select_phold_sample.py --run <run> --out <run>/phold_sample
phold proteins-predict -i <run>/phold_sample/phold_sample.faa -o <run>/phold/predict -t 3 --batch_size 1
phold proteins-compare -i <run>/phold_sample/phold_sample.faa --predictions_dir <run>/phold/predict -o <run>/phold/compare -t 3
.venv/bin/python scripts/analysis/analyze_phold.py --run <run> --phold <run>/phold
```

## 7. Cross-machine check

```bash
.venv/bin/phageannoscope run --manifest scripts/reproducibility/reference_genomes.csv --config configs/run_linux.json --workers 3 --out results/cross_machine_check --raw-root data/raw/callers_check --no-caller-cache
.venv/bin/python scripts/reproducibility/cross_machine_check.py --run results/cross_machine_check --reference scripts/reproducibility/reference_calls.parquet
```

The script exits with code 0 only if the CDS calls are identical to the reference.

## 8. Record age (exploratory)

```bash
P=.venv/bin/python
$P scripts/record_age/fetch_genbank_metadata.py --mode headers  --accessions results/descriptive/per_genome_summary.csv --out data/raw/ncbi/headers
$P scripts/record_age/fetch_genbank_metadata.py --mode features --accessions results/descriptive/per_genome_summary.csv --out data/raw/ncbi/features
$P scripts/record_age/parse_record_metadata.py data/raw/ncbi/headers data/raw/ncbi/features results/record_age/record_metadata.csv.gz
$P scripts/record_age/analyze_record_age.py results/descriptive/per_genome_summary.csv \
    results/record_age/record_metadata.csv.gz manifests/cohort_manifest_validated.csv.gz results/record_age
$P scripts/record_age/plot_record_age.py results/record_age/genome_table.csv.gz results/record_age
```

GenBank headers change when records are updated, so a later download can differ slightly from the snapshot of 2026-11-10 stored in `results/record_age/record_metadata.csv.gz`. The figure is written to `results/record_age/record_age_components.png` (`figures/fig12_record_age_components.png`).

## 9. Consensus re-annotation (exploratory)

Requires the merged run of step 5 with its PHROGs tables, and `results/record_age/genome_table.csv.gz`.

```bash
.venv/bin/python scripts/reannotation/reannotate.py results/full_cohort results/record_age/genome_table.csv.gz results/reannotation
.venv/bin/python scripts/reannotation/summarize_reannotation.py results/reannotation results/full_cohort
```

The run takes about a minute and up to 6 GB of memory. The figure is written to `results/reannotation/reannotation.png` (`figures/fig13_consensus_reannotation.png`).
