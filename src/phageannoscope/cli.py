"""Public CLI contract for the structural pipeline."""

import argparse
import json
import logging
import sys
from pathlib import Path

from . import __version__
from .config import load_config, resolve_external_tools
from .util import hardware, versions, write_json


def doctor(config=None, out=None):
    from .callers.phanotate import command

    settings = resolve_external_tools(load_config(config))
    installed = versions()
    from .callers.phanotate import command_version

    installed["phanotate"] = command_version(settings)
    report = {
        "package_version": __version__,
        "hardware": hardware(),
        "versions": installed,
        "phanotate_command": command(settings.get("phanotate_command")),
        "core_ready": all(
            installed[n] for n in ("biopython", "pyarrow", "duckdb", "pyrodigal-gv", "psutil")
        ),
        "independent_comparison_ready": bool(
            all(
                installed[n]
                for n in ("biopython", "pyarrow", "duckdb", "pyrodigal-gv", "psutil", "phanotate")
            )
            and command(settings.get("phanotate_command"))
        ),
        "config": settings,
    }
    if out:
        write_json(out, report)
    return report


def parser():
    root = argparse.ArgumentParser(description="PhageAnnoScope structural analysis")
    root.add_argument("--version", action="version", version=__version__)
    commands = root.add_subparsers(dest="command", required=True)
    p = commands.add_parser("doctor")
    p.add_argument("--config")
    p.add_argument("--out")
    p = commands.add_parser("fetch-ncbi")
    inputs = p.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--manifest")
    inputs.add_argument(
        "--accessions", nargs="+", help="Exact identifiers or a text file containing them"
    )
    p.add_argument("--out", required=True)
    p.add_argument("--events", required=True)
    p.add_argument("--limit", type=int, default=3, help="Safe initial batch limit (default:3)")
    p.add_argument("--all", action="store_true", help="Explicitly retrieve the full manifest")
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--no-resolve-hosts", action="store_true")
    p.add_argument("--batch-size", type=int, default=20)
    p.add_argument("--progress")
    p = commands.add_parser("resolve-hosts")
    p.add_argument("--manifest", required=True)
    p.add_argument("--genbank-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--data-root")
    p.add_argument("--refresh", action="store_true")
    p = commands.add_parser("validate-cohort")
    p.add_argument("--manifest", required=True)
    p.add_argument("--genbank-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--host-evidence")
    p.add_argument("--data-root")
    p = commands.add_parser("select")
    p.add_argument("--manifest", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--seed", type=int, default=30102026)
    p.add_argument("--strategy", default="protocol-length-strata")
    p.add_argument("--include")
    p.add_argument("--out", required=True)
    p = commands.add_parser("run")
    p.add_argument("--manifest", required=True)
    p.add_argument("--config")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--out", required=True)
    p.add_argument("--raw-root", default="data/raw/callers")
    p.add_argument(
        "--no-caller-cache",
        action="store_true",
        help="Execute fresh native attempts, preserving previous raw logs",
    )
    p.add_argument(
        "--diagnostic", action="store_true", help="Technical sample outside primary research N"
    )
    p = commands.add_parser("replace-selection")
    p.add_argument("--manifest", required=True)
    p.add_argument("--previous", required=True)
    p.add_argument("--run", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=30102026)
    p.add_argument("--attempted", nargs="*")
    for name in ("validate-results", "benchmark"):
        p = commands.add_parser(name)
        p.add_argument("--run", required=True)
        p.add_argument("--out")
    p = commands.add_parser("resource-report")
    p.add_argument("--run", required=True)
    p.add_argument("--population-manifest", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--available-hours", type=float)
    p = commands.add_parser("build-db")
    p.add_argument("--data-root", required=True)
    p.add_argument("--out", required=True)
    p = commands.add_parser("import-inphared")
    p.add_argument("--archive", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", required=True)
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
    try:
        if args.command == "doctor":
            result = doctor(args.config, args.out)
        elif args.command == "fetch-ncbi":
            from .acquire import fetch_ncbi

            accessions = args.accessions
            if accessions and len(accessions) == 1 and Path(accessions[0]).is_file():
                accessions = Path(accessions[0]).read_text(encoding="utf-8").split()
            result = fetch_ncbi(
                manifest=args.manifest,
                accessions=accessions,
                out=args.out,
                events=args.events,
                resume=not args.no_resume,
                limit=None if args.all else args.limit,
                resolve_hosts=not args.no_resolve_hosts,
                batch_size=args.batch_size,
                progress_path=args.progress,
            )
        elif args.command == "validate-cohort":
            from .manifest import validate_cohort

            result = validate_cohort(
                args.manifest,
                args.genbank_dir,
                args.out,
                host_evidence=args.host_evidence,
                data_root=args.data_root,
            )
        elif args.command == "resolve-hosts":
            from .hosts import resolve_hosts

            result = resolve_hosts(
                args.manifest,
                args.genbank_dir,
                args.out,
                refresh=args.refresh,
                data_root=args.data_root,
                include_records=False,
            )
        elif args.command == "select":
            from .manifest import select_manifest

            result = select_manifest(
                args.manifest, args.n, args.seed, args.strategy, args.out, args.include
            )
        elif args.command == "run":
            from .benchmark import export_benchmark
            from .pipeline import run_pipeline
            from .validation import validate_results

            log_dir = Path(args.out)
            log_dir.mkdir(parents=True, exist_ok=True)
            logging.getLogger().addHandler(
                logging.FileHandler(log_dir / "pipeline.log", encoding="utf-8")
            )
            config = load_config(args.config)
            if args.no_caller_cache:
                config["use_caller_cache"] = False
            result = run_pipeline(
                args.manifest,
                config,
                args.out,
                args.workers,
                args.diagnostic,
                args.raw_root,
            )
            result["validation"] = validate_results(args.out)["status"]
            export_benchmark(args.out)
        elif args.command == "replace-selection":
            from .replacement import replace_selection

            result = replace_selection(
                args.manifest,
                args.previous,
                args.run,
                args.out,
                seed=args.seed,
                attempted=args.attempted,
            )
        elif args.command == "validate-results":
            from .validation import validate_results

            result = validate_results(args.run, args.out)
        elif args.command == "benchmark":
            from .benchmark import export_benchmark

            result = export_benchmark(args.run, args.out)
        elif args.command == "resource-report":
            from .benchmark import resource_report

            result = resource_report(
                args.run, args.population_manifest, args.out, args.available_hours
            )
        elif args.command == "build-db":
            from .storage import build_db

            result = {"database": build_db(args.data_root, args.out)}
        else:
            from .import_inphared import import_inphared

            result = import_inphared(args.archive, args.manifest, args.out)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        if isinstance(result, dict) and (
            result.get("status") in {"failed", "incomplete"}
            or result.get("validation") == "failed"
            or result.get("failed_imports")
            or result.get("pilot_status") == "incomplete"
        ):
            raise SystemExit(1)
    except (ValueError, OSError) as error:
        logging.error("%s", error)
        raise SystemExit(2) from error


if __name__ == "__main__":
    main(sys.argv[1:])
