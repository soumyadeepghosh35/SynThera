#!/usr/bin/env python3
"""Generate doranet_config.yaml from a DORAnet Python input script."""

import argparse
import ast
import json
from pathlib import Path
import textwrap
from typing import Any


def literal_value(node: ast.AST, variables: dict[str, Any]) -> Any:
    """Read safe literal values without executing the input Python file."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        if node.id not in variables:
            raise ValueError(f"Unknown variable: {node.id}")
        return variables[node.id]
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        return [literal_value(item, variables) for item in node.elts]
    if isinstance(node, ast.Dict):
        return {
            literal_value(key, variables): literal_value(value, variables)
            for key, value in zip(node.keys, node.values)
        }
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -literal_value(node.operand, variables)
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "Path"
        and len(node.args) == 1
    ):
        return str(literal_value(node.args[0], variables))
    raise ValueError(f"Unsupported non-literal expression: {ast.unparse(node)}")


def call_name(node: ast.Call) -> str:
    if isinstance(node.func, ast.Name):
        return node.func.id
    if isinstance(node.func, ast.Attribute):
        return node.func.attr
    return ""


def parse_doranet_script(source: Path) -> dict[str, Any]:
    tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    variables: dict[str, Any] = {}

    # Resolve top-level literal assignments such as user_starters and job_name.
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            try:
                value = literal_value(statement.value, variables)
            except ValueError:
                continue
            for target in statement.targets:
                if isinstance(target, ast.Name):
                    variables[target.id] = value

    generate_args: dict[str, Any] = {}
    pathway_args: dict[str, Any] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = call_name(node)
        destination = generate_args if name == "generate_network" else pathway_args if name in {
            "one_step", "pretreat_networks"
        } else None
        if destination is None:
            continue
        for keyword in node.keywords:
            if keyword.arg is None:
                continue
            try:
                destination[keyword.arg] = literal_value(keyword.value, variables)
            except ValueError:
                pass

    required = {
        "doranetPath": variables.get("DORANET_PATH"),
        "jobName": generate_args.get("job_name", variables.get("job_name")),
        "starters": generate_args.get("starters", variables.get("user_starters")),
        "generations": generate_args.get("gen"),
        "direction": generate_args.get("direction"),
        "ruleset": generate_args.get("ruleset"),
        "maxAtoms": generate_args.get("max_atoms"),
    }
    missing = [key for key, value in required.items() if value is None]
    if missing:
        raise ValueError("Could not find required setting(s): " + ", ".join(missing))

    config = dict(required)
    config["totalGenerations"] = pathway_args.get(
        "total_generations", required["generations"]
    )
    helpers = pathway_args.get("helpers", variables.get("user_helpers"))
    if helpers is not None:
        config["helpers"] = helpers
    return config


def render_yaml(config: dict[str, Any]) -> str:
    """Render the configuration in the same layout as doranet_config.yaml."""
    starters = "\n".join(f"  - {json.dumps(smi)}" for smi in config["starters"])
    max_atoms = ", ".join(f"{element}: {limit}" for element, limit in config["maxAtoms"].items())
    helpers = textwrap.fill(
        "helpers: [" + ", ".join(json.dumps(smi) for smi in config.get("helpers", [])) + "]",
        width=88,
        subsequent_indent="          ",
        break_long_words=False,
        break_on_hyphens=False,
    )
    return f'''# DORAnet job configuration.
# Edit this file only; runDORAnet.py is reusable across starters.

# Path to the local DORAnet checkout (added to sys.path at runtime).
doranetPath: {config["doranetPath"]}

# Prefix for all output files (molecule CSV and pathway files).
jobName: {config["jobName"]}

# Starter molecule(s), as SMILES. Provide one or several.
starters:
{starters}

# Network generation settings.
generations: {config["generations"]}
direction: {config["direction"]}
ruleset: {config["ruleset"]}
maxAtoms: {{{max_atoms}}}

# Pathway generation depth. Defaults to 'generations' if omitted.
totalGenerations: {config["totalGenerations"]}

# Cofactor / helper molecules. Optional: remove to use the default set in runDORAnet.py.
{helpers}
'''


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract DORAnet settings from a Python script into YAML."
    )
    parser.add_argument("python_file", type=Path, help="DORAnet .py input file")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("doranet_config.yaml"),
        help="output path (default: doranet_config.yaml)",
    )
    args = parser.parse_args()

    if not args.python_file.is_file():
        parser.error(f"input file not found: {args.python_file}")

    try:
        config = parse_doranet_script(args.python_file)
    except (OSError, SyntaxError, ValueError) as exc:
        parser.error(str(exc))

    args.output.write_text(render_yaml(config), encoding="utf-8")
    print(f"Created {args.output}")


if __name__ == "__main__":
    main()
