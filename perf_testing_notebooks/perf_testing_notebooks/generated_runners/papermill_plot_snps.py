#!/usr/bin/env python3
"""CSV-driven Papermill runner generated for one notebook."""

import argparse
import csv
import json
from pathlib import Path

import papermill as pm

NOTEBOOK_PATH = Path('/Users/katie.barr/malariagen_api_utils/perf_testing_notebooks/notebooks/plot_snps.ipynb').resolve()
PARAMETER_NAMES = ['API_KEEP_CALLSITE_VALUE', 'API_CALL_PARAMETERS']
OMIT = object()


def parse_value(value):
    value = value.strip()
    if not value:
        return OMIT
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return value


def main():
    parser = argparse.ArgumentParser(
        description=f"Run {NOTEBOOK_PATH.name} once for each CSV row."
    )
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=NOTEBOOK_PATH.parent / "executed_runs",
    )
    args = parser.parse_args()
    config_path = args.config.resolve()
    output_dir = args.output_dir.resolve()

    with config_path.open(newline="") as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []
        if not columns:
            raise ValueError(f"{config_path} is missing a CSV header.")
        if len(columns) != len(set(columns)):
            raise ValueError(f"{config_path} contains duplicate column names.")
        unknown = set(columns) - set(PARAMETER_NAMES)
        if unknown:
            raise ValueError(
                f"Unsupported notebook parameters in {config_path}: "
                f"{', '.join(sorted(unknown))}. Valid parameters: "
                f"{', '.join(PARAMETER_NAMES)}"
            )

        runs = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(
                    f"Line {line_number} of {config_path} has more values than "
                    "there are CSV columns."
                )
            parameters = {}
            for name, value in row.items():
                if value is None:
                    continue
                parsed = parse_value(value)
                if parsed is not OMIT:
                    parameters[name] = parsed
            if not parameters:
                raise ValueError(
                    f"Line {line_number} of {config_path} has no parameter values."
                )
            runs.append(parameters)

    if not runs:
        raise ValueError(f"No configurations found in {config_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for index, parameters in enumerate(runs, start=1):
        output_path = output_dir / f"{NOTEBOOK_PATH.stem}_{index:03d}.ipynb"
        print(f"[{index}/{len(runs)}] {parameters}")
        pm.execute_notebook(
            input_path=str(NOTEBOOK_PATH),
            output_path=str(output_path),
            parameters=parameters,
            cwd=str(NOTEBOOK_PATH.parent),
            progress_bar=False,
            log_output=False,
        )

    print(f"Completed {len(runs)} runs in {output_dir}")


if __name__ == "__main__":
    main()
