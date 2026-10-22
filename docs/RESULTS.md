# Results (full cohort)

Generated from the analysis outputs in `results/`; numbers only, no interpretation.

## Cohort

- Genomes analysed: 17,882; primary analysis: 17,862; relatedness blocks: 4,490
- Primary stop-loci: 2,300,621; supported by both callers: 1,858,746; by one caller: 441,875 (PHANOTATE only 387,648, Pyrodigal-gv only 54,227)
- Median protein length (aa): one caller 45, both callers 137

## Pairwise agreement (pooled Jaccard)

| Pair | Exact CDS | Shared stop |
|---|---:|---:|
| genbank – phanotate | 0.614 | 0.798 |
| genbank – pyrodigal_gv | 0.749 | 0.880 |
| phanotate – pyrodigal_gv | 0.635 | 0.808 |

## H1: caller support and functional uncertainty

- No specific PHROGs function: one caller 95.3% (95% CI 94.8–95.8%), both callers 64.0% (62.5–65.5%); difference -31.3 pp (-32.8–-29.9)
- Adjusted OR (both vs one caller) = 0.486 (95% CI 0.426–0.555), p = 1.3e-26; n = 2,300,621 loci, 4,490 clusters

| Sensitivity analysis | OR | 95% CI |
|---|---:|---|
| Excluding genomes with ribosomal slippage | 0.505 | 0.437–0.583 |
| ORFs ≥ 30 aa | 0.487 | 0.426–0.556 |
| Excluding genomes without GenBank CDS | 0.491 | 0.429–0.561 |
| ORFs < 100 aa only | 0.453 | 0.357–0.573 |
| ORFs ≥ 100 aa only | 0.514 | 0.462–0.572 |
| One genome per relatedness block | 0.570 | 0.546–0.594 |

## H2: length and overlaps

- Single-caller support: proteins < 100 aa 38.9%, ≥ 100 aa 3.7%

Logistic GEE for single-caller support (cluster-robust 95% CI):

| Covariate | OR | 95% CI |
|---|---:|---|
| log(protein length), per unit | 0.026 | 0.021–0.031 |
| Genome GC fraction, per unit | 11.515 | 6.429–20.625 |
| log(genome length), per unit | 1.087 | 0.994–1.189 |
| Same-strand overlap with another locus (≥ 1 nt) | 0.605 | 0.569–0.644 |
| Antisense overlap with another locus (≥ 1 nt) | 43.915 | 34.985–55.125 |

- Same-strand overlap ≥ 30 nt: OR 6.37

## H3: Phold (probability sample of 600 proteins)

- Informative assignment among proteins without a specific function: 18.2% (95% CI 14.8–22.3%), n = 400; design-weighted 19.0%

| Caller support | Protein length | PHROGs function | n | Phold informative | Category agreement |
|---|---|---|---:|---:|---:|
| 1 caller | ≥ 100 aa | none (target) | 100 | 26% | — |
| 1 caller | ≥ 100 aa | specific (control) | 50 | 88% | 86% |
| 1 caller | < 100 aa | none (target) | 100 | 2% | — |
| 1 caller | < 100 aa | specific (control) | 50 | 52% | 50% |
| 2 callers | ≥ 100 aa | none (target) | 100 | 38% | — |
| 2 callers | ≥ 100 aa | specific (control) | 50 | 94% | 94% |
| 2 callers | < 100 aa | none (target) | 100 | 7% | — |
| 2 callers | < 100 aa | specific (control) | 50 | 56% | 56% |

## H4: ORFs absent from GenBank with a specific PHROGs function

| Layer | ORFs | Genomes |
|---|---:|---:|
| Main layer: complete GenBank annotation, no complex context | 6,399 | 3,986 |
| Sparse GenBank annotation (completeness ratio < 0.8) | 2,292 | 182 |
| Next to frameshift / partial GenBank CDS | 1,536 | 1,265 |
| Genomes without any GenBank CDS | 27,483 | 716 |

- Main layer by category: DNA, RNA and nucleotide metabolism 1,921, lysis 1,233, other 737, moron, auxiliary metabolic gene and host takeover 714, head and packaging 692, tail 591

## H1 robustness (adjusted OR, 95% CI)

| Section | Stratum | OR | 95% CI | n |
|---|---|---:|---|---:|
| threshold_variants | E≤1e-5, cov≥0.5 | 0.486 | 0.426–0.555 | 2,300,621 |
| threshold_variants | E≤1e-3 | 0.493 | 0.435–0.559 | 2,300,621 |
| threshold_variants | E≤1e-10 | 0.467 | 0.405–0.537 | 2,300,621 |
| threshold_variants | cov≥0.3 | 0.496 | 0.441–0.558 | 2,300,621 |
| threshold_variants | cov≥0.7 | 0.503 | 0.430–0.588 | 2,300,621 |
| length_classes | <50 aa | 0.170 | 0.096–0.301 | 306,643 |
| length_classes | 50–99 aa | 0.548 | 0.454–0.661 | 705,696 |
| length_classes | 100–199 aa | 0.635 | 0.543–0.742 | 674,845 |
| length_classes | ≥200 aa | 0.319 | 0.236–0.431 | 613,437 |
| host_genera | Mycobacterium | 0.436 | 0.323–0.588 | 283,380 |
| host_genera | Escherichia | 0.483 | 0.295–0.790 | 340,429 |
| host_genera | Pseudomonas | 0.477 | 0.404–0.563 | 167,363 |
| host_genera | Klebsiella | 0.529 | 0.367–0.764 | 161,455 |
| host_genera | Salmonella | 0.668 | 0.457–0.977 | 152,227 |
| host_genera | Vibrio | 0.558 | 0.374–0.835 | 111,186 |
| host_genera | Staphylococcus | 0.257 | 0.126–0.525 | 102,530 |
| host_genera | Gordonia | 0.804 | 0.633–1.021 | 69,578 |
| host_genera | Microbacterium | 0.256 | 0.040–1.647 | 44,199 |
| host_genera | Bacillus | 1.034 | 0.578–1.850 | 81,946 |

## Descriptive

- GenBank CDS: 1,921,788; named 654,825; 'hypothetical' 1,265,958 (65.9%); 'hypothetical' with a specific PHROGs function 176,230 (13.9%)
- Unique representative proteins: 1,122,903
- CDS per kb (median): genbank 1.505, phanotate 1.695, pyrodigal_gv 1.453
- Alternative starts PHANOTATE vs Pyrodigal-gv: median 16 codons (n = 243,706 pairs)
