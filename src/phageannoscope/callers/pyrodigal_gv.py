"""Adapter for the installed ViralGeneFinder, never ordinary Pyrodigal substitution."""

from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    SeqFeature,
    SimpleLocation,
)

from ..genbank_io import normalize_cds
from ..locations import map_rotated_segments


def normalize_gene(genome, gene, source_run_id, ordinal=0, rotation_shift=0):
    """Map one real native Gene; partial flags refer to genomic input edges."""
    start0, end0 = gene.begin - 1, gene.end
    if gene.strand not in (1, -1):
        raise ValueError("Native Pyrodigal returned an unknown strand")
    segments = map_rotated_segments(
        [{"start0": start0, "end0": end0}], rotation_shift, genome["length_bp"], gene.strand
    )
    partial5 = gene.partial_begin if gene.strand == 1 else gene.partial_end
    partial3 = gene.partial_end if gene.strand == 1 else gene.partial_begin
    parts = [SimpleLocation(s["start0"], s["end0"], strand=gene.strand) for s in segments]
    if partial5:
        p = parts[0]
        parts[0] = SimpleLocation(
            BeforePosition(p.start) if gene.strand == 1 else p.start,
            AfterPosition(p.end) if gene.strand == -1 else p.end,
            strand=gene.strand,
        )
    if partial3:
        p = parts[-1]
        parts[-1] = SimpleLocation(
            BeforePosition(p.start) if gene.strand == -1 else p.start,
            AfterPosition(p.end) if gene.strand == 1 else p.end,
            strand=gene.strand,
        )
    location = parts[0] if len(parts) == 1 else CompoundLocation(parts, operator="join")
    native_translation = gene.translate(include_stop=True)
    feature = SeqFeature(
        location,
        type="CDS",
        qualifiers={
            "transl_table": [str(gene.translation_table)],
            "translation": [native_translation],
            "note": [
                f"pyrodigal-gv native partial_begin={gene.partial_begin}; partial_end={gene.partial_end}"
            ],
        },
    )
    call = normalize_cds(
        genome, feature, "pyrodigal_gv", source_run_id, ordinal, rotation_shift=rotation_shift
    )
    call.update(
        native_begin=gene.begin,
        native_end=gene.end,
        native_partial_begin=bool(gene.partial_begin),
        native_partial_end=bool(gene.partial_end),
        native_translation=native_translation,
        native_nucleotide=gene.sequence(),
        stop_coordinate_convention="native-begin/end-1-based-inclusive-including-terminal-stop",
    )
    return call


def predict(genome, source_run_id, rotation_shift=0, params=None):
    import pyrodigal_gv

    params = {"meta": True, **(params or {})}
    if not params["meta"]:
        raise ValueError("The Pyrodigal-gv adapter requires native metagenomic mode")
    finder = pyrodigal_gv.ViralGeneFinder(**params)
    actual = {
        key: getattr(finder, key, None)
        for key in (
            "meta",
            "closed",
            "mask",
            "min_mask",
            "min_gene",
            "min_edge_gene",
            "max_overlap",
            "backend",
        )
    }
    actual["viral_only"] = params.get("viral_only", False)
    sequence = genome["sequence"]
    if not sequence:
        raise ValueError("Pyrodigal-gv requires a nonempty genome sequence")
    rotation_shift %= len(sequence)
    rotated = sequence[rotation_shift:] + sequence[:rotation_shift]
    genes = finder.find_genes(rotated.encode("ascii"))
    selected_model = str(getattr(genes, "metagenomic_bin", None))
    calls = []
    for ordinal, gene in enumerate(genes):
        call = normalize_gene(genome, gene, source_run_id, ordinal, rotation_shift)
        call["selected_model"] = selected_model
        call["automatic_code_selection"] = True
        call["code_selection_uncertainty"] = "not_quantified_by_adapter"
        calls.append(call)
    return {
        "calls": calls,
        "actual_params": actual,
        "selected_model": str(getattr(genes, "metagenomic_bin", None)),
        "selected_translation_table": getattr(
            getattr(genes, "training_info", None), "translation_table", None
        ),
    }
