"""Generate a CSV-driven Papermill runner for a malariagen-data notebook."""

import argparse
import ast
import json
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_NOTEBOOKS_DIR = SCRIPT_DIR.parent / "malariagen-data-python" / "notebooks"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "generated_runners"


def _parameter_names(notebook_path: Path) -> list[str]:
    notebook = json.loads(notebook_path.read_text())
    names = []
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        if "parameters" not in cell.get("metadata", {}).get("tags", []):
            continue
        try:
            module = ast.parse("".join(cell.get("source", [])))
        except SyntaxError as exc:
            raise ValueError(
                f"Cannot parse Papermill parameter cell in {notebook_path}: {exc}"
            ) from exc

        for node in module.body:
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for target in targets:
                    if not isinstance(target, ast.Name):
                        raise ValueError(
                            f"Unsupported parameter assignment in {notebook_path}; "
                            "parameter cells must assign simple variable names."
                        )
                    if target.id not in names:
                        names.append(target.id)

    if not names:
        raise ValueError(
            f"No parameters-tagged code cell with simple assignments found in "
            f"{notebook_path}."
        )
    return names


def _resolve_notebook(name: str, notebooks_dir: Path) -> Path:
    notebook_name = name if name.endswith(".ipynb") else f"{name}.ipynb"
    matches = list(notebooks_dir.rglob(notebook_name))
    if not matches:
        raise FileNotFoundError(
            f"Notebook {notebook_name!r} was not found under {notebooks_dir}."
        )
    if len(matches) > 1:
        match_list = "\n".join(f"  {path}" for path in matches)
        raise ValueError(
            f"Notebook name {notebook_name!r} is ambiguous; matches:\n{match_list}"
        )
    return matches[0].resolve()


def _render_runner(notebook_path: Path, parameter_names: list[str]) -> str:
    template = r'''#!/usr/bin/env python3
"""CSV-driven Papermill runner generated for one notebook."""

import argparse
import csv
import json
from pathlib import Path

import papermill as pm

NOTEBOOK_PATH = Path(__NOTEBOOK_PATH__).resolve()
PARAMETER_NAMES = __PARAMETER_NAMES__
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
'''
    return (
        template.replace("__NOTEBOOK_PATH__", repr(str(notebook_path)))
        .replace("__PARAMETER_NAMES__", repr(parameter_names))
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a CSV-driven Papermill runner for a notebook under "
            "malariagen-data-python/notebooks."
        )
    )
    parser.add_argument("notebook", help="Notebook filename or unique name without .ipynb.")
    parser.add_argument(
        "--notebooks-dir",
        type=Path,
        default=DEFAULT_NOTEBOOKS_DIR,
        help=f"Notebook directory (default: {DEFAULT_NOTEBOOKS_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for generated runner scripts (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args()

    notebook_path = _resolve_notebook(args.notebook, args.notebooks_dir.resolve())
    parameter_names = _parameter_names(notebook_path)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    runner_path = args.output_dir / f"papermill_{notebook_path.stem}.py"
    runner_path.write_text(_render_runner(notebook_path, parameter_names))
    print(f"Generated {runner_path}")
    print(f"Notebook parameters: {', '.join(parameter_names)}")


if __name__ == "__main__":
    main()
