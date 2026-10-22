"""One native Pyrodigal-gv worker process per genome/rotation."""

import json
import sys
import time
from pathlib import Path

from ..util import write_json
from .pyrodigal_gv import predict


def main():
    request = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    started = time.process_time()
    result = predict(**request)
    result["worker_cpu_seconds"] = time.process_time() - started
    write_json(sys.argv[2], result)


if __name__ == "__main__":
    main()
