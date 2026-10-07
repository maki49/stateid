"""ABACUS input loading for state identification."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List


DEFAULT_ABACUS_FILENAMES = (
    "abacus_state_data.json",
    "stateid_abacus.json",
    "abacus.json",
)


def _as_float(value: Any, *, field: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {field}: {value!r}") from exc


def _normalized_projections(raw: Dict[str, Any]) -> Dict[str, float]:
    projections = {k: _as_float(v, field=f"projection[{k}]") for k, v in raw.items()}
    total = sum(v for v in projections.values() if v > 0.0)
    if total <= 0.0:
        return {}
    return {k: max(v, 0.0) / total for k, v in projections.items()}


def _coerce_record(record: Dict[str, Any]) -> Dict[str, Any]:
    projections = record.get("projections", {})
    if not isinstance(projections, dict):
        raise ValueError("projections must be a mapping")

    return {
        "kpoint_index": int(record.get("kpoint_index", 0)),
        "band_index": int(record.get("band_index", 0)),
        "spin": str(record.get("spin", "none")),
        "energy_ev": _as_float(record.get("energy_ev"), field="energy_ev"),
        "occupation": _as_float(record.get("occupation"), field="occupation"),
        "projections": _normalized_projections(projections),
    }


def load_abacus_evidence(path: str | Path) -> Dict[str, Any]:
    """Load normalized ABACUS-derived state evidence.

    Expected JSON shape:
      {
        "fermi_energy_ev": 0.0,
        "bands": [
          {
            "kpoint_index": 0,
            "band_index": 1,
            "spin": "up",
            "energy_ev": -0.2,
            "occupation": 1.0,
            "projections": {"s": 0.1, "p": 0.9}
          }
        ]
      }
    """

    p = Path(path)
    if p.is_dir():
        matches: Iterable[Path] = (p / name for name in DEFAULT_ABACUS_FILENAMES)
        p = next((candidate for candidate in matches if candidate.exists()), None)  # type: ignore[assignment]
        if p is None:
            raise FileNotFoundError(
                f"No ABACUS evidence JSON found in {path!s}; tried {DEFAULT_ABACUS_FILENAMES}"
            )

    raw = json.loads(Path(p).read_text(encoding="utf-8"))
    if "bands" not in raw or not isinstance(raw["bands"], list):
        raise ValueError("ABACUS evidence must include list field: bands")

    fermi = _as_float(raw.get("fermi_energy_ev", 0.0), field="fermi_energy_ev")
    bands: List[Dict[str, Any]] = [_coerce_record(entry) for entry in raw["bands"]]

    return {
        "source": str(p),
        "fermi_energy_ev": fermi,
        "bands": bands,
    }
