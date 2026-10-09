import itertools
from pathlib import Path
import papermill as pm

# 1. Define parameter grid
matrix = {
    "dataset_size": 2,
    "start_base": 240000,
    "size_increment": 10000,
}

# Generate all parameter combinations
keys, values = zip(*matrix.items())
param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

# 2. Directory to store executed notebook artifacts
runs_dir = Path("executed_runs")
runs_dir.mkdir(exist_ok=True)

# 3. Execute runs systematically
template_nb = "test_capture_mem_output.ipynb"

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