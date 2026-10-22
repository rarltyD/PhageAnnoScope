from pathlib import Path

import pytest

from phageannoscope.fasta_io import read_fasta, resolve_fasta_identifier, validate_fasta_genbank
from phageannoscope.genbank_io import read_genbank

FIXTURES = Path(__file__).parent / "fixtures" / "synthetic"


def test_matching_exact_version_sequence_and_pipe_alias():
    genbank = read_genbank(FIXTURES / "synthetic.gb")["genomes"]
    fasta = read_fasta(FIXTURES / "synthetic.fna")
    pipe = read_fasta(FIXTURES / "pipe.fna")
    alias = read_fasta(FIXTURES / "alias.fna", {"local_fixture": "SYNTHETIC_A.7"})
    assert validate_fasta_genbank(fasta, genbank)
    assert validate_fasta_genbank(pipe, genbank)
    assert validate_fasta_genbank(alias, genbank)
    assert fasta[0]["genome_id"] == genbank[0]["genome_id"]


def test_version_and_sha_mismatch():
    genbank = read_genbank(FIXTURES / "synthetic.gb")["genomes"]
    fasta = read_fasta(FIXTURES / "synthetic.fna")
    changed = [{**fasta[0], "accession_version": "SYNTHETIC_A.8"}]
    with pytest.raises(ValueError, match="version mismatch"):
        validate_fasta_genbank(changed, genbank)
    changed = [{**fasta[0], "sequence_sha256": "0" * 64}]
    with pytest.raises(ValueError, match="SHA mismatch"):
        validate_fasta_genbank(changed, genbank)


def test_duplicate_header_is_an_error(tmp_path):
    text = (FIXTURES / "synthetic.fna").read_text()
    path = tmp_path / "duplicate.fna"
    path.write_text(text + text)
    with pytest.raises(ValueError, match="Duplicate FASTA header"):
        read_fasta(path)


@pytest.mark.parametrize(
    "header",
    [
        "lcl|SYNTHETIC_A.7",
        "anything_SYNTHETIC_A.7_extra",
        "SYNTHETIC_A",
        "ref|SYNTHETIC_A.7|gb|SYNTHETIC_OTHER.3|",
    ],
)
def test_substrings_unversioned_and_ambiguous_headers_rejected(header):
    with pytest.raises(ValueError):
        resolve_fasta_identifier(header)
