#!/usr/bin/env python3
"""Inspect parameterized notebooks and create CSV configurations for their API calls.

Examples:
    python papermill_csv_config.py list
    python papermill_csv_config.py describe plot_pca
    python papermill_csv_config.py create plot_pca \\
        --set 'cell_8_call_1_ag3_pca.n_components=10' \\
        --sweep 'cell_8_call_1_ag3_pca.site_mask=["default","gamb_colu"]'
    python papermill_csv_config.py create plot_pca \\
        --validate-metadata --output plot_pca_config.csv
"""

import argparse
import ast
import csv
import fnmatch
import itertools
import json
import re
from copy import deepcopy
from pathlib import Path

# Notes from Copilot:
# Use --set 'CALL_KEY.PARAMETER=JSON_VALUE' for a fixed override; use --sweep for values that should create separate runs. If you want to set an API argument to a list, use --set so the list is treated as one value.

#Optional --validate-metadata checks sample selections and cohort columns against Malariagen metadata; you can provide a dataset URL with --api-url ag3=URL. For PCA, it resolves a short cohort name such as admin1_year to cohort_admin1_year when metadata confirms that column is populated. It doesn’t validate every kind of API argument, such as regions or masks.

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_NOTEBOOKS_DIR = SCRIPT_DIR / "perf_testing_notebooks" / "notebooks"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "perf_testing_notebooks"
API_PARAMETERS_NAME = "API_CALL_PARAMETERS"


def _load_api_parameters(notebook_path: Path) -> dict:
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        if "parameters" not in cell.get("metadata", {}).get("tags", []):
            continue

        try:
            module = ast.parse("".join(cell.get("source", [])))
        except SyntaxError as exc:
            raise ValueError(
                f"Cannot parse the parameters cell in {notebook_path}: {exc}"
            ) from exc

        for node in module.body:
            if not isinstance(node, (ast.Assign, ast.AnnAssign)):
                continue
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if not any(
                isinstance(target, ast.Name)
                and target.id == API_PARAMETERS_NAME
                for target in targets
            ):
                continue
            try:
                value = ast.literal_eval(node.value)
            except (ValueError, TypeError) as exc:
                raise ValueError(
                    f"{API_PARAMETERS_NAME} in {notebook_path} must be a "
                    "literal mapping."
                ) from exc
            if not isinstance(value, dict) or not all(
                isinstance(key, str) and isinstance(params, dict)
                for key, params in value.items()
            ):
                raise ValueError(
                    f"{API_PARAMETERS_NAME} in {notebook_path} must map "
                    "call keys to argument mappings."
                )
            return value

    raise ValueError(
        f"No parameters-tagged {API_PARAMETERS_NAME} mapping found in "
        f"{notebook_path}."
    )


def _iter_notebooks(notebooks_dir: Path):
    for path in sorted(notebooks_dir.rglob("*.ipynb")):
        if ".ipynb_checkpoints" in path.parts:
            continue
        try:
            parameters = _load_api_parameters(path)
        except (ValueError, json.JSONDecodeError):
            continue
        if parameters:
            yield path, parameters


def _resolve_notebook(name: str, notebooks_dir: Path) -> Path:
    candidate = Path(name)
    if candidate.suffix != ".ipynb":
        candidate = candidate.with_suffix(".ipynb")

    if candidate.is_absolute():
        matches = [candidate] if candidate.is_file() else []
    else:
        direct = notebooks_dir / candidate
        if direct.is_file():
            matches = [direct]
        else:
            matches = [
                path for path, _ in _iter_notebooks(notebooks_dir)
                if path.name == candidate.name
            ]

    if not matches:
        raise FileNotFoundError(
            f"Notebook {name!r} was not found under {notebooks_dir}. "
            "Use the path shown by the 'list' command."
        )
    if len(matches) > 1:
        choices = "\n".join(
            f"  {path.relative_to(notebooks_dir)}" for path in matches
        )
        raise ValueError(
            f"Notebook name {name!r} is ambiguous. Use one of:\n{choices}"
        )
    return matches[0].resolve()


def _parse_json_value(value: str, option: str):
    try:
        return json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON value for {option}: {value!r}. "
            'Quote string values, e.g. --set "call.param=\\"value\\"".'
        ) from exc


def _parse_assignment(text: str, option: str):
    if "=" not in text:
        raise ValueError(
            f"{option} expects CALL_KEY.PARAMETER=JSON_VALUE; got {text!r}."
        )
    target, value = text.split("=", 1)
    call_pattern, separator, parameter_name = target.rpartition(".")
    if not separator or not call_pattern or not parameter_name:
        raise ValueError(
            f"{option} expects CALL_KEY.PARAMETER=JSON_VALUE; got {text!r}."
        )
    return call_pattern, parameter_name, _parse_json_value(value, option)


def _resolve_target(call_pattern: str, parameter_name: str, defaults: dict):
    matching_calls = [
        call_key
        for call_key in defaults
        if fnmatch.fnmatchcase(call_key, call_pattern)
    ]
    if not matching_calls:
        raise ValueError(
            f"No API call matches {call_pattern!r}. Use 'describe' to list call keys."
        )

    invalid = [
        call_key
        for call_key in matching_calls
        if parameter_name not in defaults[call_key]
    ]
    if invalid:
        raise ValueError(
            f"Parameter {parameter_name!r} is not available for "
            f"{', '.join(invalid)}."
        )
    return [(call_key, parameter_name) for call_key in matching_calls]


def _prepare_overrides(defaults: dict, set_values: list[str], sweep_values: list[str]):
    constants = {}
    sweeps = []

    for text in set_values:
        pattern, parameter_name, value = _parse_assignment(text, "--set")
        targets = _resolve_target(pattern, parameter_name, defaults)
        for target in targets:
            if target in constants:
                raise ValueError(f"Duplicate --set override for {target[0]}.{target[1]}.")
            constants[target] = value

    for text in sweep_values:
        pattern, parameter_name, values = _parse_assignment(text, "--sweep")
        if not isinstance(values, list) or not values:
            raise ValueError(
                f"--sweep value must be a non-empty JSON array; got {values!r}."
            )
        targets = _resolve_target(pattern, parameter_name, defaults)
        for target in targets:
            if target in constants or any(target in targets for targets, _ in sweeps):
                raise ValueError(
                    f"Parameter {target[0]}.{target[1]} has more than one override."
                )
        sweeps.append((targets, values))

    return constants, sweeps


def _configuration_rows(defaults: dict, constants: dict, sweeps: list):
    sweep_combinations = itertools.product(
        *(values for _, values in sweeps)
    ) if sweeps else [()]

    for combination in sweep_combinations:
        settings = deepcopy(defaults)
        for (call_key, parameter_name), value in constants.items():
            settings[call_key][parameter_name] = value
        for (targets, _), value in zip(sweeps, combination):
            for call_key, parameter_name in targets:
                settings[call_key][parameter_name] = value
        yield settings


def _describe(notebook_path: Path, defaults: dict):
    print(f"Notebook: {notebook_path}")
    print(f"API calls: {len(defaults)}")
    for call_key, parameters in defaults.items():
        print(f"\n{call_key}")
        for name, value in parameters.items():
            rendered = json.dumps(value, ensure_ascii=False, default=repr)
            print(f"  {name} = {rendered}")


def _parse_dataset_urls(entries: list[str]) -> dict:
    urls = {}
    for entry in entries:
        alias, separator, url = entry.partition("=")
        if not separator or not alias or not url:
            raise ValueError(
                f"Invalid --api-url {entry!r}; expected CLIENT=URL, "
                "for example ag3=https://example/data."
            )
        urls[alias.lower()] = url
    return urls


def _validate_metadata(settings: dict, dataset_urls: dict):
    import inspect

    import malariagen_data

    clients = {}
    metadata_cache = {}
    validated_calls = 0
    api_classes = {
        name.lower(): value
        for name, value in inspect.getmembers(malariagen_data, inspect.isclass)
        if re.fullmatch(r"[A-Z][a-z]+\d+", name)
        and callable(getattr(value, "sample_metadata", None))
    }

    for call_key, parameters in settings.items():
        if not any(
            name in parameters
            for name in ("sample_sets", "sample_query", "sample_indices", "cohorts")
        ):
            continue
        aliases = [
            alias
            for alias in api_classes
            if re.search(rf"(?:^|_){re.escape(alias)}(?:_|$)", call_key)
        ]
        if not aliases:
            raise ValueError(
                f"Cannot identify a malariagen-data client for metadata "
                f"validation of {call_key!r}."
            )
        if len(aliases) != 1:
            raise ValueError(
                f"Cannot determine the data client for API call {call_key!r}."
            )
        alias = aliases[0]
        if alias not in clients:
            url = dataset_urls.get(alias)
            clients[alias] = (
                api_classes[alias](url=url) if url else api_classes[alias]()
            )

        metadata_parameters = {
            name: parameters.get(name)
            for name in (
                "sample_sets",
                "sample_query",
                "sample_query_options",
                "sample_indices",
            )
            if name in parameters
        }
        cache_key = (alias, json.dumps(metadata_parameters, sort_keys=True, default=repr))
        if cache_key not in metadata_cache:
            metadata_cache[cache_key] = clients[alias].sample_metadata(
                **metadata_parameters
            )
        metadata = metadata_cache[cache_key]
        validated_calls += 1

        cohorts = parameters.get("cohorts")
        if not isinstance(cohorts, str):
            continue
        possible_columns = [cohorts]
        if not cohorts.startswith("cohort_"):
            possible_columns.append(f"cohort_{cohorts}")
        valid_column = next(
            (
                column
                for column in possible_columns
                if column in metadata.columns and metadata[column].notna().any()
            ),
            None,
        )
        if valid_column is None:
            raise ValueError(
                f"{call_key}: cohorts={cohorts!r} is not a populated column "
                "in sample metadata for the configured sample selection."
            )
        if call_key.endswith("_pca"):
            parameters["cohorts"] = valid_column

    if not validated_calls:
        raise ValueError(
            "Metadata validation was requested, but no API calls in this "
            "notebook expose sample-selection or cohort parameters."
        )


def _create_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Create Papermill CSV configurations from notebook API-call "
            "parameter cells."
        ),
        epilog=(
            "Examples:\n"
            "  papermill_csv_config.py list\n"
            "  papermill_csv_config.py describe plot_pca\n"
            "  papermill_csv_config.py create plot_pca "
            "--set 'cell_8_call_1_ag3_pca.n_components=10'\n"
            "  papermill_csv_config.py create plot_pca "
            "--sweep 'cell_8_call_1_ag3_pca.n_components=[5,10,20]'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--notebooks-dir",
        type=Path,
        default=DEFAULT_NOTEBOOKS_DIR,
        help=f"Directory containing notebooks (default: {DEFAULT_NOTEBOOKS_DIR})",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="List notebooks with API parameter cells.")

    describe = commands.add_parser(
        "describe", help="Show API call keys, available arguments, and defaults."
    )
    describe.add_argument("notebook", help="Notebook name or path relative to notebooks-dir.")

    create = commands.add_parser(
        "create", help="Write one CSV row per configuration."
    )
    create.add_argument("notebook", help="Notebook name or path relative to notebooks-dir.")
    create.add_argument(
        "--output",
        type=Path,
        help="Output CSV (default: <notebooks-dir>/../<notebook>_config.csv).",
    )
    create.add_argument(
        "--set",
        dest="set_values",
        action="append",
        default=[],
        metavar="CALL.PARAM=JSON",
        help="Set one argument value; may be repeated. Supports * in CALL.",
    )
    create.add_argument(
        "--sweep",
        dest="sweep_values",
        action="append",
        default=[],
        metavar="CALL.PARAM=[JSON,...]",
        help="Create a sweep over a non-empty JSON array; may be repeated.",
    )
    create.add_argument(
        "--validate-metadata",
        action="store_true",
        help="Check sample selections/cohort columns against malariagen-data metadata.",
    )
    create.add_argument(
        "--api-url",
        dest="api_urls",
        action="append",
        default=[],
        metavar="CLIENT=URL",
        help="Override the data URL during metadata checks; may be repeated.",
    )
    return parser


def main(argv=None):
    parser = _create_parser()
    args = parser.parse_args(argv)
    notebooks_dir = args.notebooks_dir.resolve()

    if args.command == "list":
        for path, parameters in _iter_notebooks(notebooks_dir):
            print(
                f"{path.relative_to(notebooks_dir)} "
                f"({len(parameters)} API calls)"
            )
        return

    notebook_path = _resolve_notebook(args.notebook, notebooks_dir)
    defaults = _load_api_parameters(notebook_path)

    if args.command == "describe":
        _describe(notebook_path, defaults)
        return

    constants, sweeps = _prepare_overrides(
        defaults,
        args.set_values,
        args.sweep_values,
    )
    dataset_urls = _parse_dataset_urls(args.api_urls)
    output_path = args.output
    if output_path is None:
        output_path = DEFAULT_OUTPUT_DIR / f"{notebook_path.stem}_config.csv"
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = list(_configuration_rows(defaults, constants, sweeps))
    if args.validate_metadata:
        for settings in rows:
            _validate_metadata(settings, dataset_urls)

    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=[API_PARAMETERS_NAME],
            lineterminator="\n",
        )
        writer.writeheader()
        for settings in rows:
            writer.writerow(
                {API_PARAMETERS_NAME: json.dumps(settings, allow_nan=False)}
            )

    print(f"Wrote {len(rows)} configurations to {output_path}")


if __name__ == "__main__":
    main()
