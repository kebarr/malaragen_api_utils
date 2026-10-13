import ast
import csv
import json
from pathlib import Path

notebook_path = Path("notebooks/plot_pca.ipynb")
notebook = json.loads(notebook_path.read_text())

parameter_source = "".join(notebook["cells"][0]["source"])
tree = ast.parse(parameter_source)
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

for name, options in defaults.items():
    if name.endswith("_ag3_pca"):
        options["n_components"] = 10
        options["random_seed"] = 123

with Path("perf_testing_notebooks/plot_pca_config.csv").open(
    "w", newline=""
) as f:
    writer = csv.DictWriter(f, fieldnames=["API_CALL_PARAMETERS"])
    writer.writeheader()
    writer.writerow({"API_CALL_PARAMETERS": json.dumps(defaults)})
