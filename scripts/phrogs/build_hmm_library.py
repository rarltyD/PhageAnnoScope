"""Build a HMMER3 profile library from the PHROGs v4 M50 MSAs.

The distributed PHROGs HMMs are HH-suite ``.hhm`` files, which HMMER/pyhmmer
cannot read.  Profiles are therefore rebuilt from the official MSAs with
``pyhmmer.plan7.Builder`` default parameters.  The archive is read as a
stream; nothing is extracted to disk and the source archive is never modified.

Usage (bounded test first):
    nice -n 10 python scripts/phrogs/build_hmm_library.py --msa data/db/phrogs_v4/MSA_phrogs.tar.gz \
        --out data/db/phrogs_v4/phrogs_v4_m50.h3m --limit 200
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import tarfile
import time
from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, MSAFile


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msa", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args(argv)
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f"Refusing to overwrite existing {out}")
    alphabet = Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    temporary = out.with_name(out.name + ".tmp")
    started = time.perf_counter()
    built, failed = 0, []
    with tarfile.open(args.msa, "r:gz") as archive, temporary.open("wb") as handle:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".fma"):
                continue
            name = re.sub(r"\.fma$", "", Path(member.name).name)
            data = archive.extractfile(member).read()
            try:
                with MSAFile(io.BytesIO(data), format="afa", digital=True, alphabet=alphabet) as f:
                    msa = f.read()
                msa.name = name.encode()
                hmm, _, _ = builder.build_msa(msa, background)
                hmm.write(handle, binary=True)
                built += 1
            except Exception as error:  # keep going; failures are reported
                failed.append({"phrog": name, "error": f"{type(error).__name__}: {error}"})
            if args.limit and built + len(failed) >= args.limit:
                break
    temporary.replace(out)
    report = {
        "source_archive": str(Path(args.msa).resolve()),
        "source_sha256": hashlib.sha256(Path(args.msa).read_bytes()).hexdigest(),
        "builder": "pyhmmer.plan7.Builder defaults",
        "pyhmmer_version": pyhmmer.__version__,
        "hmmer_version": pyhmmer.LIBHMMER_VERSION if hasattr(pyhmmer, "LIBHMMER_VERSION") else None,
        "profiles_built": built,
        "failed": failed,
        "limit": args.limit,
        "wall_seconds": round(time.perf_counter() - started, 2),
        "output": str(out.resolve()),
        "output_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
    }
    out.with_suffix(".build.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "failed"} | {"n_failed": len(failed)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
