"""Core evidence-based state identification logic."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List


@dataclass(frozen=True)
class IdentifiedState:
    kpoint_index: int
    band_index: int
    spin: str
    label: str
    confidence: float
    evidence: List[str]


def _dominant_orbital(projections: Dict[str, float]) -> tuple[str, float]:
    if not projections:
        return "unknown", 0.0
    orbital, weight = max(projections.items(), key=lambda x: x[1])
    return orbital, weight


def _state_label(energy: float, fermi: float, occupation: float) -> str:
    if occupation >= 0.95 and energy <= fermi:
        return "valence"
    if occupation <= 0.05 and energy >= fermi:
        return "conduction"
    return "frontier"


def _confidence(energy: float, fermi: float, occupation: float, dom_weight: float, label: str) -> float:
    occ_signal = occupation if label == "valence" else (1.0 - occupation if label == "conduction" else 0.5)
    proximity = max(0.0, 1.0 - min(abs(energy - fermi), 2.0) / 2.0)
    score = 0.45 * occ_signal + 0.35 * dom_weight + 0.20 * proximity
    return max(0.0, min(1.0, round(score, 3)))


def identify_states(abacus_evidence: Dict[str, Any], *, energy_window_ev: float = 1.0) -> List[Dict[str, Any]]:
    """Identify electronic states near the Fermi level with explicit evidence."""

    fermi = float(abacus_evidence["fermi_energy_ev"])
    identified: List[IdentifiedState] = []

    for band in abacus_evidence["bands"]:
        energy = float(band["energy_ev"])
        if abs(energy - fermi) > energy_window_ev:
            continue

        occupation = float(band["occupation"])
        label = _state_label(energy, fermi, occupation)
        dom_orbital, dom_weight = _dominant_orbital(dict(band.get("projections", {})))
        confidence = _confidence(energy, fermi, occupation, dom_weight, label)

        evidence = [
            f"energy_ev={energy:.6f}",
            f"fermi_energy_ev={fermi:.6f}",
            f"occupation={occupation:.4f}",
            f"dominant_orbital={dom_orbital}",
            f"dominant_weight={dom_weight:.4f}",
        ]

        identified.append(
            IdentifiedState(
                kpoint_index=int(band["kpoint_index"]),
                band_index=int(band["band_index"]),
                spin=str(band.get("spin", "none")),
                label=f"{label} ({dom_orbital}-dominated)",
                confidence=confidence,
                evidence=evidence,
            )
        )

    return [
        {
            "kpoint_index": state.kpoint_index,
            "band_index": state.band_index,
            "spin": state.spin,
            "label": state.label,
            "confidence": state.confidence,
            "evidence": state.evidence,
        }
        for state in sorted(identified, key=lambda s: (s.kpoint_index, s.band_index, s.spin))
    ]
