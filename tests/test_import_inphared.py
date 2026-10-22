"""Local synthetic archives only; no INPHARED download is performed."""

import hashlib
import io
import tarfile
import zipfile

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from phageannoscope.import_inphared import import_inphared, inventory_archive


def tar(path, name, data=b"SYNTHETIC TEST ONLY", *, link=False):
    with tarfile.open(path, "w") as handle:
        member = tarfile.TarInfo(name)
        if link:
            member.type = tarfile.SYMTYPE
            member.linkname = "../../outside"
            handle.addfile(member)
        else:
            member.size = len(data)
            handle.addfile(member, io.BytesIO(data))


@pytest.mark.parametrize(
    "name",
    [
        "../escape.fna",
        "/absolute.fna",
        "C:\\escape.fna",
        "safe/../../escape.fna",
        "safe\\..\\escape.fna",
    ],
)
def test_tar_traversal_refused(tmp_path, name):
    path = tmp_path / "bad.tar"
    tar(path, name)
    with pytest.raises(ValueError, match="Unsafe archive"):
        inventory_archive(path)


def test_links_refused_before_import_output(tmp_path):
    path = tmp_path / "bad.tar"
    tar(path, "safe/link.fna", link=True)
    with pytest.raises(ValueError, match="links"):
        import_inphared(path, [], tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_zip_traversal_and_budget_refused(tmp_path):
    path = tmp_path / "bad.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("../escape.fna", "bad")
    with pytest.raises(ValueError, match="Unsafe"):
        inventory_archive(path)
    tar(tmp_path / "safe.tar", "safe/file.fna", b"12345")
    with pytest.raises(ValueError, match="budget"):
        inventory_archive(tmp_path / "safe.tar", max_bytes=4)


def test_sha_mismatch_quarantined(tmp_path):
    path = tmp_path / "synthetic.tar"
    tar(path, "synthetic/test.fna", b">SYNTH001.2 SYNTHETIC TEST ONLY\nATGAAATAA\n")
    summary = import_inphared(
        path,
        [{"accession_version": "SYNTH001.2", "sequence_sha256": "incorrect"}],
        tmp_path / "out",
    )
    assert summary["matched_records"] == 0
    assert summary["imported_cds_calls"] == 0
    assert summary["quarantine"][0]["reason"] == "SEQUENCE_SHA_MISMATCH"


def test_pipe_header_exact_match_and_unversioned_rejected(tmp_path):
    path = tmp_path / "synthetic.tar"
    tar(
        path,
        "synthetic/test.fna",
        b">ref|SYNTH001.2| SYNTHETIC TEST ONLY\nATGAAATAA\n>SYNTH001 SYNTHETIC TEST ONLY\nATGAAATAA\n",
    )
    sha = hashlib.sha256(b"ATGAAATAA").hexdigest()
    summary = import_inphared(
        path, [{"accession_version": "SYNTH001.2", "sequence_sha256": sha}], tmp_path / "out"
    )
    assert summary["matched_records"] == 1 and summary["quarantined_records"] == 1
    assert summary["matches"][0]["accession_version"] == "SYNTH001.2"
    assert summary["quarantine"][0]["reason"] == "ACCESSION_VERSION_NOT_MATCHED"


def test_conflicting_duplicate_header_member_not_imported(tmp_path):
    path = tmp_path / "synthetic.tar"
    tar(
        path,
        "synthetic/test.fna",
        b">SYNTH001.2 SYNTHETIC TEST ONLY\nATGAAATAA\n>SYNTH001.2 SYNTHETIC TEST ONLY\nATGCCCTAA\n",
    )
    sha = hashlib.sha256(b"ATGAAATAA").hexdigest()
    summary = import_inphared(
        path, [{"accession_version": "SYNTH001.2", "sequence_sha256": sha}], tmp_path / "out"
    )
    assert summary["matched_records"] == 0
    assert summary["quarantine"][0]["reason"] == "DUPLICATE_RECORD_IDENTITY"


def test_genbank_multi_record_import_filters_unmatched_genomes(tmp_path):
    records = []
    for identity in ("SYNTH001.2", "SYNTH002.8"):
        record = SeqRecord(
            Seq("ATGAAATAA"),
            id=identity,
            name=identity.split(".")[0],
            description="SYNTHETIC TEST ONLY",
            annotations={"molecule_type": "DNA"},
        )
        record.features = [
            SeqFeature(
                SimpleLocation(0, 9, strand=1), type="CDS", qualifiers={"transl_table": ["11"]}
            )
        ]
        records.append(record)
    handle = io.StringIO()
    SeqIO.write(records, handle, "genbank")
    path = tmp_path / "synthetic.tar"
    tar(path, "synthetic/test.gbk", handle.getvalue().encode())
    sha = hashlib.sha256(b"ATGAAATAA").hexdigest()
    summary = import_inphared(
        path, [{"accession_version": "SYNTH001.2", "sequence_sha256": sha}], tmp_path / "out"
    )
    assert summary["matched_records"] == 1 and summary["imported_cds_calls"] == 1
    assert [g["accession_version"] for g in summary["annotations"]["genomes"]] == ["SYNTH001.2"]
    assert {f["source"] for f in summary["annotations"]["features"]} == {"inphared"}
