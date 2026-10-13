import ast
import csv
import json
from copy import deepcopy
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = SCRIPT_DIR / "notebooks" / "plot_pca.ipynb"
OUTPUT_CSV_PATH = SCRIPT_DIR / "plot_pca_config_cohorts.csv"
N_COMPONENTS_VALUES = [5, 10, 20]
COHORT_VALUES = ["country", "admin1_year", "taxon"]
MAX_COHORT_SIZE = 20


def main():
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    tree = ast.parse("".join(notebook["cells"][0]["source"]))
    defaults = next(
        ast.literal_eval(node.value)
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name)
            and target.id == "API_CALL_PARAMETERS"
            for target in node.targets
        )
    )

    pca_call_keys = [key for key in defaults if key.endswith("_ag3_pca")]
    if not pca_call_keys:
        raise ValueError(f"No Ag3 PCA calls found in {NOTEBOOK_PATH}")

    with OUTPUT_CSV_PATH.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["API_CALL_PARAMETERS"])
        writer.writeheader()

        for n_components in N_COMPONENTS_VALUES:
            for cohorts in COHORT_VALUES:
                settings = deepcopy(defaults)
                for key in pca_call_keys:
                    settings[key]["n_components"] = n_components
                    settings[key]["cohorts"] = cohorts
                    settings[key]["max_cohort_size"] = MAX_COHORT_SIZE

                writer.writerow(
                    {"API_CALL_PARAMETERS": json.dumps(settings)}
                )

    print(
        f"Wrote {len(N_COMPONENTS_VALUES) * len(COHORT_VALUES)} configurations "
        f"to {OUTPUT_CSV_PATH}"
    )


if __name__ == "__main__":
    main()
