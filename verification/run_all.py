"""Run the project's mathematical verification programs.

Use --write-results DIRECTORY to save fresh JSON outputs. This never
updates the committed expected/ records unless that directory is chosen
explicitly. Assertions in the programs enforce the mathematical checks.
"""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys


PROGRAMS = (
    "verify_code_parameters.py",
    "verify_concatenated_family.py",
    "verify_physical_channel.py",
)


def compare(actual, expected, path: str, tolerance: float = 0.0) -> None:
    """Check the full recorded result, ignoring only runtime-version metadata."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict), f"{path}: expected an object"
        assert set(actual) == set(expected), f"{path}: result keys changed"
        for key, value in expected.items():
            if key != "environment":
                compare(actual[key], value, f"{path}.{key}", tolerance)
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), f"{path}: list length changed"
        for index, (got, wanted) in enumerate(zip(actual, expected)):
            compare(got, wanted, f"{path}[{index}]", tolerance)
    elif isinstance(expected, bool):
        assert actual is expected, f"{path}: {actual!r} != {expected!r}"
    elif isinstance(expected, (int, float)):
        assert isinstance(actual, (int, float)) and not isinstance(actual, bool), f"{path}: expected a number"
        assert math.isfinite(actual) and abs(actual - expected) <= tolerance, (
            f"{path}: {actual!r} differs from {expected!r} by more than {tolerance}"
        )
    else:
        assert actual == expected, f"{path}: {actual!r} != {expected!r}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    if args.write_results:
        args.write_results.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for program in PROGRAMS:
        completed = subprocess.run([sys.executable, str(root / program)],
                                   check=True, capture_output=True, text=True)
        result = json.loads(completed.stdout)
        assert result["checks_passed"] is True
        expected_path = root / "expected" / f"{Path(program).stem}.json"
        expected = json.loads(expected_path.read_text())
        tolerance = 2e-12 if program == "verify_physical_channel.py" else 0.0
        compare(result, expected, program, tolerance)
        outputs[program] = result
        if args.write_results:
            (args.write_results / f"{Path(program).stem}.json").write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n"
            )
    print(json.dumps({"expected_records_matched": True, "results": outputs}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
