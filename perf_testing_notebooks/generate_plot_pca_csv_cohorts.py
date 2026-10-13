import ast
import csv
import json
from copy import deepcopy
from pathlib import Path

import malariagen_data


SCRIPT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = SCRIPT_DIR / "notebooks" / "plot_pca.ipynb"
OUTPUT_CSV_PATH = SCRIPT_DIR / "plot_pca_config_cohorts.csv"
N_COMPONENTS_VALUES = [5, 10, 20]
COHORT_VALUES = ["country", "admin1_year", "taxon"]
MAX_COHORT_SIZE = 20
AG3_URL = "simplecache::gs://vo_agam_release_master_us_central1"


def load_notebook_defaults():
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    tree = ast.parse("".join(notebook["cells"][0]["source"]))
    return next(
        ast.literal_eval(node.value)
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name)
            and target.id == "API_CALL_PARAMETERS"
            for target in node.targets
        )
    )


def find_valid_cohort_columns(api, pca_parameters):
    metadata = api.sample_metadata(
        sample_sets=pca_parameters.get("sample_sets"),
        sample_query=pca_parameters.get("sample_query"),
        sample_query_options=pca_parameters.get("sample_query_options"),
        sample_indices=pca_parameters.get("sample_indices"),
    )

    valid_columns = {}
    for cohort in COHORT_VALUES:
        column = cohort
        if column not in metadata.columns and f"cohort_{cohort}" in metadata.columns:
            column = f"cohort_{cohort}"
        if column in metadata.columns and metadata[column].notna().any():
            valid_columns[cohort] = column
    return valid_columns


def main():
    defaults = load_notebook_defaults()
    pca_call_keys = [key for key in defaults if key.endswith("_ag3_pca")]
    if not pca_call_keys:
        raise ValueError(f"No Ag3 PCA calls found in {NOTEBOOK_PATH}")

    api = malariagen_data.Ag3(
        AG3_URL,
        simplecache=dict(cache_storage=str(SCRIPT_DIR / "gcs_cache")),
    )
    valid_columns_by_call = {
        key: find_valid_cohort_columns(api, defaults[key])
        for key in pca_call_keys
    }
    supported_cohorts = [
        cohort
        for cohort in COHORT_VALUES
        if any(cohort in columns for columns in valid_columns_by_call.values())
    ]
    if not supported_cohorts:
        raise ValueError(
            "None of the configured cohort columns are present with values in "
            "sample metadata for the notebook's Ag3 PCA calls."
        )

    for cohort in COHORT_VALUES:
        supported_calls = [
            key for key in pca_call_keys
            if cohort in valid_columns_by_call[key]
        ]
        if supported_calls:
            print(f"{cohort!r} is valid for: {', '.join(supported_calls)}")
        else:
            print(f"Skipping {cohort!r}: no matching populated metadata column.")

    with OUTPUT_CSV_PATH.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=["API_CALL_PARAMETERS"],
            lineterminator="\n",
        )
        writer.writeheader()

        for n_components in N_COMPONENTS_VALUES:
            for cohorts in supported_cohorts:
                settings = deepcopy(defaults)
                for key in pca_call_keys:
                    settings[key]["n_components"] = n_components
                    cohort_column = valid_columns_by_call[key].get(cohorts)
                    settings[key]["cohorts"] = cohort_column
                    settings[key]["max_cohort_size"] = (
                        MAX_COHORT_SIZE if cohort_column is not None else None
                    )

                writer.writerow(
                    {"API_CALL_PARAMETERS": json.dumps(settings)}
                )

    print(
        f"Wrote {len(N_COMPONENTS_VALUES) * len(supported_cohorts)} configurations "
        f"to {OUTPUT_CSV_PATH}"
    )


if __name__ == "__main__":
    main()
