"""CLI entrypoint for stateid."""

from __future__ import annotations

import argparse
import json
import sys

from .abacus import load_abacus_evidence
from .core import identify_states


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evidence-based electronic-state identification (ABACUS)")
    parser.add_argument("input", help="ABACUS evidence JSON file or directory containing one")
    parser.add_argument("--energy-window-ev", type=float, default=1.0, help="Energy window around Fermi level")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    evidence = load_abacus_evidence(args.input)
    states = identify_states(evidence, energy_window_ev=args.energy_window_ev)
    payload = {
        "source": evidence["source"],
        "fermi_energy_ev": evidence["fermi_energy_ev"],
        "identified_states": states,
    }
    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
