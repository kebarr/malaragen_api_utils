Cell 1:

# Parameters cell — Papermill overrides these at runtime
dataset_size = 100_000
batch_size = 512
algorithm = "sort"
output_csv = "memory_benchmarks.csv"


Cell 2:

import time

def workload(size, batch, algo):
    # Simulated memory-intensive task
    data = [i % 1000 for i in range(size)]
    if algo == "sort":
        return sorted(data)
    elif algo == "set_lookup":
        s = set(data)
        return [x for x in data if x in s]
    return data


Cell2:
import itertools
from pathlib import Path
import papermill as pm

Call 3:
import os
import csv
import filelock  # pip install filelock (prevents race conditions if running concurrently)
from memory_profiler import memory_usage

# memory_usage returns peak RAM in MiB during execution
start_time = time.perf_counter()
peak_memory_mib, result = memory_usage(
    (workload, (dataset_size, batch_size, algorithm)),
    max_usage=True,
    retval=True
)
wall_time_sec = time.perf_counter() - start_time

# Record dictionary
record = {
    "dataset_size": dataset_size,
    "batch_size": batch_size,
    "algorithm": algorithm,
    "peak_memory_mib": round(peak_memory_mib, 2),
    "runtime_sec": round(wall_time_sec, 4),
}

# Append safely using a file lock
lock = filelock.FileLock(f"{output_csv}.lock")
with lock:
    file_exists = os.path.isfile(output_csv)
    with open(output_csv, mode="a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=record.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(record)

print(f"Recorded: {record}")


external_Script:
# 1. Define parameter grid
matrix = {
    "dataset_size": [100_000, 500_000, 1_000_000],
    "batch_size": [256, 1024],
    "algorithm": ["sort", "set_lookup"],
}

# Generate all parameter combinations
keys, values = zip(*matrix.items())
param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

# 2. Directory to store executed notebook artifacts
runs_dir = Path("executed_runs")
runs_dir.mkdir(exist_ok=True)

# 3. Execute runs systematically
template_nb = "benchmark_template.ipynb"

for idx, params in enumerate(param_combinations, start=1):
    output_nb = runs_dir / f"run_{idx}_{params['algorithm']}_{params['dataset_size']}.ipynb"
    print(f"[{idx}/{len(param_combinations)}] Running with: {params}")

    pm.execute_notebook(
        input_path=template_nb,
        output_path=str(output_nb),
        parameters={
            **params,
            "output_csv": "memory_benchmarks.csv",
        },
        progress_bar=False,
        log_output=False,
    )

print("\nBenchmark complete! Results written to memory_benchmarks.csv.")
