# Study protocol and amendments

Protocol v1.0 was registered on 2026-10-16, before any gene calling or hypothesis testing. Its machine-readable parameters are in `configs/protocol_v1.json`. Each later change is listed below as a dated amendment with its rationale, and every amendment was introduced before the affected results were interpreted.

## 1. Aim and scope of inference

The study asks three questions about DNA bacteriophages:

- How consistently do two algorithmically different gene callers predict CDS?
- How is single- versus two-caller support related to ORF length, overlaps and functional characterisation?
- What additional functional assignments does structure-informed annotation provide?

NCBI GenBank serves as the baseline annotation, not as an experimental gold standard. INPHARED2 serves for cohort selection, not as evidence of GenBank errors. A predicted ORF, even one found by both callers, does not prove that a protein exists.

**Hypotheses**

- **H1.** The share of proteins without a specific function differs between ORFs supported by both callers and ORFs supported by one caller.
- **H2.** Short and overlapping ORFs are more often supported by a single caller.
- **H3.** Structure-informed annotation (Phold) assigns functions to some proteins that lack one.
- **H4.** Some ORFs absent from GenBank are homologous to functionally characterised phage protein families.

## 2. Data sources

- **INPHARED2.** Official TSV export of 2026-10-16: 25 columns, 37,387 records. The SHA-256 of the raw file is recorded in `manifests/`.
- **NCBI GenBank.** Records were retrieved per accession.version, each with its retrieval date and the SHA-256 of the raw file.
- **PHROGs v4.** Annotation table and M50 multiple sequence alignments.
- **Phold.** Version 1.2.6 with its bundled database and the ProstT5 model.

## 3. Cohort

### 3.1 Metadata screen

A record entered the primary-candidate set (28,664 records) if it met all of the following:

- Baltimore group I or II;
- CheckV quality Complete or High-quality;
- not flagged as a MIUVIG genome fragment.

### 3.2 Eligibility checks on the GenBank record

1. **Identity.** The record has a unique accession.version, a non-empty sequence and a length that agrees with the metadata.
2. **Host.** The phage infects bacteria. Records with an archaeal or eukaryotic host are excluded; records with an unresolved host are held.
3. **Completeness.** The genome is complete or high quality. Segmented genomes form a separate layer, and records flagged as proviral are held.
4. **Ambiguous bases.** Sequences with N, or (amendment A5) other IUPAC ambiguity codes, enter sensitivity analyses only.
5. **Genetic code and recoding.** Genomes with a non-standard genetic code or genome-level recoding are held. The narrower treatment of ribosomal slippage is set by amendments A2 and A2c.
6. **Duplicates.** Exact sequence duplicates are represented by one deterministic accession.

### 3.3 Relatedness

The relatedness block of a genome is its TaxMyPhage genus, otherwise the NCBI genus from the metadata, otherwise the genome itself. These blocks define the clusters for GEE and for the cluster bootstrap. No ANI clustering was performed.

### 3.4 Sample order

Pilot (100), nested 1,000 and nested 5,000 samples were drawn by length strata, ordered by SHA-256(seed, accession.version) with seed 30102026. The full cohort was then processed in prespecified nested batches; see "Scale" below.

## 4. Structural comparison

Coordinates are stored as 0-based half-open segments, ordered in the transcription direction. Only CDS features are compared. Pairs of calls are classified as exact CDS, alternative start (same stop, different start), contained with a distinct stop, partial overlap, antisense overlap, or no correspondence.

Calls that cannot be compared in a simple pairwise analysis are flagged and kept: partial, pseudo, exception, transl_except, order or non-contiguous joins, and incompatible genetic codes. For circular genomes, each caller is additionally run on the sequence rotated by L/2. Origin-spanning models and models sensitive to the rotation are excluded from the primary analysis.

The primary unit of analysis is the **stop-locus**: the union of PHANOTATE and Pyrodigal-gv models that share a stop codon on the same strand. Its representative model is the longest complete model. GenBank is a baseline, not a vote. Prokka and Prodigal are not counted as independent of Pyrodigal-gv.

## 5. Functional annotation

Every distinct protein sequence is searched once, and the result is propagated to all models that encode it. The search uses HMMER (pyhmmer) against PHROGs v4 profiles. The primary filter, set before analysis, is:

- E ≤ 1e-5;
- query coverage ≥ 0.5;
- profile coverage ≥ 0.5.

Each protein receives one status: specific function, PHROG of unknown function, no informative hit, or ambiguous (near-tied hits with conflicting known categories). Sensitivity variants vary the thresholds: E 1e-3 and 1e-10; coverage 0.3 and 0.7.

**Phold (H3)** was run on a stratified probability sample of 600 unique proteins:

- strata: caller support (1 or 2) × length (< 100 aa or ≥ 100 aa);
- per stratum: 100 proteins without a specific function plus 50 controls with a specific function;
- at most 2 proteins per genome;
- inclusion probabilities recorded for weighting.

An informative Phold assignment requires E ≤ 1e-3; a confirmation layer additionally requires query and target coverage ≥ 0.7. Phold reuses PHROGs-derived reference information, so the two are not independent.

## 6. Statistics

- **H1.** Logistic GEE for "no specific function". Covariates: support (both callers vs one), a spline of log(length), GC content, log(genome length) and any overlap. Independence working correlation with cluster-robust SE by relatedness block.
- **H2.** Logistic GEE for single-caller support. Covariates: log(length), GC content, log(genome length), same-strand overlap and antisense overlap. A spline variant shows the shape of the length effect.
- **Proportions.** 95% CIs from 2,000 bootstrap replicates of relatedness blocks.
- **Multiplicity.** Holm correction over the two primary tests.
- **Sensitivity analyses.** Excluding slippage genomes; ORFs ≥ 30 aa; excluding genomes without GenBank CDS; length strata; one genome per block; overlap ≥ 30 nt; the PHROGs threshold variants.

## Amendments

**A1. Host domain (2026-10-18).** The host domain is taken from the INPHARED2 `metadata_host` genus, resolved with exact NCBI Taxonomy queries. All 304 names were queried; raw XML responses and their SHA-256 are kept. For names that occur in more than one kingdom, an additional query is restricted to the Bacteria and Archaea subtrees. Hosts in Archaea or Eukaryota are excluded. Placeholders and names absent from NCBI Taxonomy are held.

**A2. Ribosomal slippage (2026-10-18).** Genomes whose only recoding feature is `/ribosomal_slippage` are eligible, because a programmed frameshift does not change the stop-codon system. The affected CDS remain excluded from pairwise comparison. Genomes with `/transl_except` or `/exception` remain on hold. A sensitivity analysis excludes all slippage genomes.

**A2b / A2c. Code corrections (2026-10-18).**

- A2b: `/ribosomal_slippage` CDS are now treated as complex biology. Before this fix, the in-frame product of a frameshifted gene appeared as an artefactual "additional" ORF.
- A2c: such CDS no longer place the whole genome on hold.

Both changes have tests.

**A3. PHROGs profiles (2026-10-18).** The distributed PHROGs HMMs are in HH-suite format. HMMER3 profiles (38,880) were rebuilt from the official M50 alignments with `pyhmmer.plan7.Builder` defaults. The resulting E-values are not directly comparable with HH-suite or MMseqs2 searches.

**A4. H4 layers (2026-10-18).** Candidate ORFs are split into four layers:

- baseline_clean: genome with GenBank CDS, no complex or partial GenBank context;
- complex or partial context;
- sparse annotation (see A6);
- genome without any GenBank CDS.

**A5. Ambiguous bases (2026-10-18).** IUPAC ambiguity codes other than N follow the rule for N: such genomes enter sensitivity analyses only. This applied to 276 genomes. Pre-drawn samples were filtered, not re-drawn.

**A6. GenBank completeness (2026-10-18).** `genbank_completeness_ratio` is the number of GenBank CDS divided by the number of loci supported by both callers. The main H4 layer requires a ratio ≥ 0.8. This rule was set from the distribution of the ratio before the H4 list was interpreted.

**Scale (2026-10-18).** All 17,882 eligible genomes were processed in prespecified nested batches of about 5,000 genomes each (same seed). Gene calling is cached per genome and PHROGs results are cached per protein, so no work is repeated between batches.

**Computational changes (2026-10-21/22). None of these changes the estimands.**

- A quadratic `list_contains` join in result validation was replaced by an equivalent equi-join on unnested call IDs. Validation of the 4,909-genome batch gave the same verdict (passed).
- Batch 3 was assembled from its verified per-genome checkpoints. The run identity differed only in the WSL kernel string; code and tool versions were identical.
- The GEE for the full cohort uses a memory-lean IRLS implementation with a cluster sandwich estimator (`scripts/analysis/lowmem_gee.py`). It is identical to statsmodels GEE with independence working correlation and robust SE. Checked on the 4,909-genome batch: OR, 95% CI and p agree to 6 decimals.
- Cross-machine check: 17 genomes recomputed from scratch on a second machine gave 7,247/7,247 identical CDS calls.

## Limitations

- There is no experimental gold standard, so agreement is measured, not accuracy.
- Relatedness is approximated by taxonomy rather than ANI.
- The host domain comes from curated metadata.
- PHROGs and Phold share reference information.
- H1 cannot distinguish genuine but uncharacterised single-caller ORFs from false predictions.
- Additional ORFs are hypotheses that require experimental validation.
