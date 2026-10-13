import ast
import csv
import json
from pathlib import Path

notebook_path = Path("notebooks/plot_snps.ipynb")
notebook = json.loads(notebook_path.read_text())

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

plot_call_keys = [
    key for key in defaults
    if key.endswith("_ag3_plot_snps_track") or key.endswith("_ag3_plot_snps")
]

with Path("perf_testing_notebooks/plot_snps_config.csv").open(
    "w", newline=""
) as f:
    writer = csv.DictWriter(f, fieldnames=["API_CALL_PARAMETERS"])
    writer.writeheader()

    for max_snps in [50_000, 100_000, 200_000]:
        settings = json.loads(json.dumps(defaults))  # fresh copy for each row

        for key in plot_call_keys:
            if "max_snps" in settings[key]:
                settings[key]["max_snps"] = max_snps

        writer.writerow({"API_CALL_PARAMETERS": json.dumps(settings)})