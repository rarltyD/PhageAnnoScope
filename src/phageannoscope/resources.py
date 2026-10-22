"""Best-effort resource probes, independent of biological job success.

An unavailable or interrupted probe never becomes a measured zero. A failed
sample invalidates that metric for this observation, rather than reporting an
earlier sample as if monitoring had remained available.
"""

import shutil
import threading

import psutil


def probe(producer, name):
    """Return a measured value and its absence reason without raising probe errors."""
    try:
        value = producer()
    except Exception as error:
        # Platform probes can raise OSError/NotImplementedError as well as psutil.Error.
        return None, f"{name}: {type(error).__name__}: {error}"
    if value is None:
        return None, f"{name}: no value available"
    return value, None


class ProcessTreeSampler:
    """Sample a process tree; failures affect measurements, never subprocess control."""

    def __init__(self, pid=None):
        self.values = {"cpu_seconds": None, "peak_rss_bytes": None}
        self.reasons = dict.fromkeys(self.values)
        self.process, reason = probe(lambda: psutil.Process(pid), "process_init")
        if reason:
            self.invalidate(reason)

    def invalidate(self, reason, metrics=None):
        for metric in metrics or self.values:
            self.values[metric] = None
            self.reasons[metric] = self.reasons[metric] or reason

    def sample(self):
        if self.process is None:
            return
        members, reason = probe(
            lambda: [self.process] + self.process.children(recursive=True), "process_tree"
        )
        if reason:
            self.invalidate(reason)
            return
        producers = {
            "peak_rss_bytes": lambda: sum(member.memory_info().rss for member in members),
            "cpu_seconds": lambda: sum(
                (times := member.cpu_times()).user + times.system for member in members
            ),
        }
        for metric, producer in producers.items():
            if self.reasons[metric]:
                continue
            value, reason = probe(producer, metric)
            if reason:
                self.invalidate(reason, [metric])
            else:
                previous = self.values[metric]
                self.values[metric] = value if previous is None else max(previous, value)

    def result(self):
        reasons = {
            metric: self.reasons[metric]
            or ("no_resource_samples" if self.values[metric] is None else None)
            for metric in self.values
        }
        observed = sum(value is not None for value in self.values.values())
        return {
            **self.values,
            "cpu_unavailable_reason": reasons["cpu_seconds"],
            "peak_rss_unavailable_reason": reasons["peak_rss_bytes"],
            "resource_monitor_status": ("available" if observed == 2 else "partial")
            if observed
            else "unavailable",
        }

    def kill_descendants(self):
        if self.process is None:
            return
        members, _ = probe(lambda: self.process.children(recursive=True), "timeout_process_tree")
        for member in members or []:
            probe(member.kill, "timeout_descendant_kill")


class BackgroundResourceMonitor:
    """Guard both monitor thread setup and its complete sampling body."""

    def __init__(self, interval=0.05, disk_path=None):
        self.sampler = ProcessTreeSampler()
        self.interval = interval
        self.stop_event = threading.Event()
        self.thread = None
        self.started = False
        self.disk_path = disk_path
        self.disk_start = self.disk_min = None
        self.disk_reason = "disk_probe_not_requested" if disk_path is None else None
        if disk_path is not None:
            self._sample_disk()

    def _sample_disk(self):
        if self.disk_path is None or self.disk_reason:
            return
        free, reason = probe(lambda: shutil.disk_usage(self.disk_path).free, "disk_free_bytes")
        if reason:
            self.disk_start = self.disk_min = None
            self.disk_reason = reason
        else:
            self.disk_start = free if self.disk_start is None else self.disk_start
            self.disk_min = free if self.disk_min is None else min(self.disk_min, free)

    def _run(self):
        try:
            self.sampler.sample()
            self._sample_disk()
            while not self.stop_event.wait(self.interval):
                self.sampler.sample()
                self._sample_disk()
        except Exception as error:
            self.sampler.invalidate(f"monitor_thread: {type(error).__name__}: {error}")
            if self.disk_path is not None:
                self.disk_reason = (
                    self.disk_reason or f"monitor_thread: {type(error).__name__}: {error}"
                )

    def start(self):
        try:
            self.thread = threading.Thread(target=self._run, daemon=True)
            self.thread.start()
            self.started = True
        except Exception as error:
            self.sampler.invalidate(f"monitor_thread_start: {type(error).__name__}: {error}")
            if self.disk_path is not None:
                self.disk_reason = (
                    self.disk_reason or f"monitor_thread_start: {type(error).__name__}: {error}"
                )

    def stop(self):
        self.stop_event.set()
        if self.started:
            try:
                self.thread.join()
            except Exception as error:
                self.sampler.invalidate(f"monitor_thread_join: {type(error).__name__}: {error}")
                if self.disk_path is not None:
                    self.disk_reason = (
                        self.disk_reason or f"monitor_thread_join: {type(error).__name__}: {error}"
                    )
        self._sample_disk()
        return {
            **self.sampler.result(),
            "disk_free_start_bytes": self.disk_start,
            "disk_free_min_bytes": self.disk_min,
            "disk_peak_delta_bytes": max(0, self.disk_start - self.disk_min)
            if self.disk_start is not None and self.disk_min is not None and not self.disk_reason
            else None,
            "disk_unavailable_reason": self.disk_reason,
            "disk_measurement_scope": "system filesystem free space sampled; unrelated process activity may contribute",
        }
