#!/usr/bin/env python3
"""Rebuild the preliminary manifest from the frozen local source, offline."""
import json
from pathlib import Path
from freeze_manifest import build

root = Path(__file__).resolve().parents[1]
raw = root / "data/raw/metadata_2026-10-16"
meta = json.loads((raw / "download_metadata.json").read_text())
build(raw / "inphared2_metadata.tsv", root,
      "INPHARED2_online_2026-10-03_export_2026-10-16",
      meta["finished_at_utc"], meta["url"])
