import argparse
import csv
import random
from pathlib import Path

import malariagen_data
import papermill as pm


SCRIPT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = SCRIPT_DIR / "perf_testing_notebooks" / "test_capture_mem_output.ipynb"
DEFAULT_CONFIG = SCRIPT_DIR / "benchmark_inputs.csv"
REQUIRED_COLUMNS = {"dataset_size", "start_base", "size_increment", "contig"}


def read_config(config_path: Path) -> list[dict[str, int | str]]:
    with config_path.open(newline="") as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames or [])
        missing_columns = REQUIRED_COLUMNS - columns
        if missing_columns:
            raise ValueError(
                f"{config_path} is missing required columns: "
                f"{', '.join(sorted(missing_columns))}"
            )

        runs = []
        for line_number, row in enumerate(reader, start=2):
            try:
                dataset_size = int(row["dataset_size"])
                start_base = int(row["start_base"])
                size_increment = int(row["size_increment"])
                contig = row["contig"].strip()
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid benchmark settings on line {line_number} of {config_path}"
                ) from exc

            if dataset_size < 1 or start_base < 0 or size_increment < 1 or not contig:
                raise ValueError(
                    f"Invalid benchmark settings on line {line_number}: "
                    "dataset_size and size_increment must be positive, "
                    "start_base must be non-negative, and contig must not be empty"
                )

            runs.append(
                {
                    "dataset_size": dataset_size,
                    "start_base": start_base,
                    "size_increment": size_increment,
                    "contig": contig,
                }
            )

    if not runs:
        raise ValueError(f"No benchmark configurations found in {config_path}")
    return runs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run per-region benchmarks from a CSV configuration."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"CSV configuration path (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for repeatable sample-set and sample selection.",
    )
    parser.add_argument(
        "--output-prefix",
        default="benchmarks",
        help="Prefix for the SNP-call and IGV-view result CSV files.",
    )
    args = parser.parse_args()

    runs = read_config(args.config)
    rng = random.Random(args.seed)
    af1 = malariagen_data.Af1()
    available_sample_sets = (
        af1.sample_sets()["sample_set"].dropna().astype(str).unique().tolist()
    )
    if not available_sample_sets:
        raise ValueError("Af1.sample_sets() returned no sample sets.")

    runs_dir = NOTEBOOK_PATH.parent / "executed_runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    for idx, config in enumerate(runs, start=1):
        sample_set = rng.choice(available_sample_sets)
        sample_metadata = af1.sample_metadata(sample_sets=sample_set)
        sample_ids = (
            sample_metadata["sample_id"].dropna().astype(str).unique().tolist()
        )
        if not sample_ids:
            raise ValueError(f"No samples found for selected sample set {sample_set!r}.")
        sample_id = rng.choice(sample_ids)

        params = {
            **config,
            "sample_set": sample_set,
            "sample_id": sample_id,
            "output_csv_prefix": args.output_prefix,
        }
        output_nb = runs_dir / (
            f"run_{idx}_{params['size_increment']}_{params['dataset_size']}.ipynb"
        )
        print(
            f"[{idx}/{len(runs)}] Running with sample_set={sample_set!r}, "
            f"sample_id={sample_id!r}, settings={config}"
        )

        pm.execute_notebook(
            input_path=str(NOTEBOOK_PATH),
            output_path=str(output_nb),
            parameters=params,
            cwd=str(NOTEBOOK_PATH.parent),
            progress_bar=False,
            log_output=False,
        )

    print(f"\nBenchmark complete. Result files use prefix {args.output_prefix!r}.")


if __name__ == "__main__":
    main()
