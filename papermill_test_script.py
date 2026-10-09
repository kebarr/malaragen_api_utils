import itertools
from pathlib import Path
import papermill as pm

# 1. Define parameter grid
# 10 sample sets + 10 size increments for each = 100 runs
matrix = {
    "dataset_size": [10],
    "start_base": [2440000],
    "size_increment": [100000],
    "sample_set":["1326-VO-UG-KAYONDO-KAJO-UG-2203", "1315-VO-NG-OMITOLA-OMOL-NG-2008",
                  "1288-VO-UG-DONNELLY-VMF00219", "1273-VO-ZM-MULEBA-VMF00176",
                  "1236-VO-TZ-OKUMU-VMF00261", "1354-VO-KE-DONNELLY-VMF00281",
                  "1231-VO-MULTI-WONDJI-VMF00043", "1235-VO-MZ-PAAIJMANS-VMF00094",
                  "1354-VO-KE-DONNELLY-VMF00281", "1236-VO-TZ-OKUMU-VMF00261"]
}


# Generate all parameter combinations
keys, values = zip(*matrix.items())
param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

# 2. Directory to store executed notebook artifacts
runs_dir = Path("executed_runs")
runs_dir.mkdir(exist_ok=True)

# 3. Execute runs systematically
template_nb = "test_capture_mem_output.ipynb"

run_name = "10_samples_10_increments"
for idx, params in enumerate(param_combinations, start=1):
    output_nb = runs_dir / f"run_{idx}_{params['size_increment']}_{params['dataset_size']}.ipynb"
    print(f"[{idx}/{len(param_combinations)}] Running with: {params}")

    pm.execute_notebook(
        input_path=template_nb,
        output_path=str(output_nb),
        parameters={
            **params,
            "output_csv_prefix": run_name,
        },
        progress_bar=False,
        log_output=False,
    )

print("\nBenchmark complete!")