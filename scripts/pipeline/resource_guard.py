"""Run a command under a resource guard (stops it safely on low memory, swap growth or heat).

Samples /proc/meminfo, swap and (if exposed) thermal zones every --interval s.
On MemAvailable < --min-avail-mb, swap used > --max-swap-mb, or temperature >
--max-temp-c it sends SIGTERM to the whole process group (then SIGKILL after
60 s).  The pipeline checkpoints per genome, so a stopped run resumes cleanly.
Samples are appended to --log as JSON lines.

Example:
    python scripts/pipeline/resource_guard.py --log logs/run.guard.jsonl -- \
        nice -n 10 .venv/bin/phageannoscope run ...
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import signal
import subprocess
import sys
import time


def meminfo():
    values = {}
    with open("/proc/meminfo") as handle:
        for line in handle:
            key, value = line.split(":", 1)
            values[key] = int(value.split()[0]) // 1024
    return values


def temperatures():
    temps = []
    for path in glob.glob("/sys/class/thermal/thermal_zone*/temp"):
        try:
            temps.append(int(open(path).read()) / 1000)
        except (OSError, ValueError):
            pass
    return temps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", required=True)
    parser.add_argument("--interval", type=float, default=10)
    parser.add_argument("--min-avail-mb", type=int, default=900)
    parser.add_argument("--max-swap-mb", type=int, default=1024)
    parser.add_argument("--max-temp-c", type=float, default=90)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    swap0 = None
    child = subprocess.Popen(command, start_new_session=True)
    reason = None
    with open(args.log, "a") as log:
        while child.poll() is None:
            m = meminfo()
            swap = m["SwapTotal"] - m["SwapFree"]
            swap0 = swap if swap0 is None else swap0
            temps = temperatures()
            sample = {"t": time.strftime("%Y-%m-%dT%H:%M:%S"), "avail_mb": m["MemAvailable"],
                      "swap_used_mb": swap, "load1": os.getloadavg()[0],
                      "max_temp_c": max(temps) if temps else None}
            log.write(json.dumps(sample) + "\n")
            log.flush()
            if m["MemAvailable"] < args.min_avail_mb:
                reason = f"MemAvailable {m['MemAvailable']} MB < {args.min_avail_mb}"
            elif swap - swap0 > args.max_swap_mb:
                reason = f"swap grew by {swap - swap0} MB"
            elif temps and max(temps) > args.max_temp_c:
                reason = f"temperature {max(temps)} C"
            if reason:
                log.write(json.dumps({"t": sample["t"], "guard_stop": reason}) + "\n")
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(60)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                break
            time.sleep(args.interval)
    code = child.wait()
    print(json.dumps({"exit_code": code, "guard_stop": reason}), file=sys.stderr)
    sys.exit(3 if reason else code)


if __name__ == "__main__":
    main()
