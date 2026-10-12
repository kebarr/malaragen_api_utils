"""Run snp_calls_to_vcf benchmarks from CSV rows.

Use method parameter names as CSV headers. Plain text cells are strings;
JSON-formatted cells preserve lists, mappings, numbers, booleans, and nulls.
Leave optional inputs blank to use the method defaults.
"""

import argparse
import csv
import inspect
import json
from pathlib import Path
from typing import Any

import malariagen_data
import papermill as pm


SCRIPT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = (
    SCRIPT_DIR / "perf_testing_notebooks" / "test_capture_mem_output_from_config.ipynb"
)
DEFAULT_CONFIG = SCRIPT_DIR / "benchmark_inputs.csv"
_OMIT = object()


def parse_cell(value: str) -> Any:
    value = value.strip()
    if not value:
        return _OMIT
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return value


def read_config(
    config_path: Path, method_signature: inspect.Signature
) -> list[dict[str, Any]]:
    with config_path.open(newline="") as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []
        if not columns:
            raise ValueError(f"{config_path} is missing a CSV header.")
        if len(columns) != len(set(columns)):
            raise ValueError(f"{config_path} contains duplicate column names.")

        accepted_parameters = set(method_signature.parameters)
        unknown_columns = set(columns) - accepted_parameters
        if unknown_columns:
            raise ValueError(
                f"{config_path} contains unsupported snp_calls_to_vcf inputs: "
                f"{', '.join(sorted(unknown_columns))}. Supported inputs are: "
                f"{', '.join(method_signature.parameters)}"
            )
        if "output_path" not in columns:
            raise ValueError(f"{config_path} must include an output_path column.")

        runs = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(
                    f"Line {line_number} of {config_path} has more values than "
                    "there are CSV columns."
                )
            kwargs = {}
            for name, value in row.items():
                if value is None:
                    continue
                parsed = parse_cell(value)
                if parsed is not _OMIT:
                    kwargs[name] = parsed
            output_path = kwargs.get("output_path")
            if not isinstance(output_path, str) or not output_path.strip():
                raise ValueError(
                    f"Line {line_number} of {config_path} must have a non-empty "
                    "output_path."
                )
            output_path = Path(output_path)
            if not output_path.is_absolute():
                output_path = config_path.parent / output_path
            kwargs["output_path"] = str(output_path.resolve())

            try:
                method_signature.bind(**kwargs)
            except TypeError as exc:
                raise ValueError(
                    f"Invalid snp_calls_to_vcf inputs on line {line_number} of "
                    f"{config_path}: {exc}"
                ) from exc
            runs.append(kwargs)

    if not runs:
        raise ValueError(f"No benchmark configurations found in {config_path}")
    return runs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark snp_calls_to_vcf calls defined in a CSV file."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"CSV configuration path (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=SCRIPT_DIR / "snp_calls_to_vcf_benchmarks.csv",
        help="CSV file for per-call timing, memory, and SNP-count results.",
    )
    args = parser.parse_args()

    af1 = malariagen_data.Af1()
    method_signature = inspect.signature(af1.snp_calls_to_vcf)
    runs = read_config(args.config.resolve(), method_signature)
    output_csv = args.output_csv.resolve()

    runs_dir = NOTEBOOK_PATH.parent / "executed_runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    for index, kwargs in enumerate(runs, start=1):
        output_nb = runs_dir / f"run_{index:03d}_snp_calls_to_vcf.ipynb"
        print(
            f"[{index}/{len(runs)}] Exporting {kwargs['output_path']} "
            f"for region {kwargs['region']}"
        )
        pm.execute_notebook(
            input_path=str(NOTEBOOK_PATH),
            output_path=str(output_nb),
            parameters={
                "snp_calls_to_vcf_kwargs": kwargs,
                "result_csv": str(output_csv),
            },
            cwd=str(NOTEBOOK_PATH.parent),
            progress_bar=False,
            log_output=False,
        )

    print(f"\nBenchmark complete. Results written to {output_csv}.")


if __name__ == "__main__":
    main()
