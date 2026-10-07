import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from stateid.abacus import load_abacus_evidence
from stateid.core import identify_states


SAMPLE = {
    "fermi_energy_ev": 0.0,
    "bands": [
        {
            "kpoint_index": 0,
            "band_index": 1,
            "spin": "up",
            "energy_ev": -0.1,
            "occupation": 1.0,
            "projections": {"p": 0.8, "s": 0.2},
        },
        {
            "kpoint_index": 0,
            "band_index": 2,
            "spin": "up",
            "energy_ev": 0.15,
            "occupation": 0.0,
            "projections": {"d": 0.7, "p": 0.3},
        },
        {
            "kpoint_index": 1,
            "band_index": 3,
            "spin": "down",
            "energy_ev": 2.5,
            "occupation": 0.0,
            "projections": {"d": 1.0},
        },
    ],
}


class StateIdTests(unittest.TestCase):
    def test_load_and_normalize_abacus_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "abacus_state_data.json"
            path.write_text(json.dumps(SAMPLE), encoding="utf-8")

            evidence = load_abacus_evidence(td)

        self.assertEqual(evidence["fermi_energy_ev"], 0.0)
        self.assertEqual(len(evidence["bands"]), 3)
        self.assertAlmostEqual(evidence["bands"][0]["projections"]["p"], 0.8)

    def test_identify_states_near_fermi(self) -> None:
        states = identify_states(
            {
                "source": "sample",
                "fermi_energy_ev": SAMPLE["fermi_energy_ev"],
                "bands": SAMPLE["bands"],
            },
            energy_window_ev=0.5,
        )

        self.assertEqual(len(states), 2)
        self.assertTrue(states[0]["label"].startswith("valence"))
        self.assertTrue(states[1]["label"].startswith("conduction"))
        self.assertGreaterEqual(states[0]["confidence"], 0.0)
        self.assertIn("dominant_orbital=", " ".join(states[0]["evidence"]))

    def test_cli_outputs_json(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "abacus_state_data.json"
            path.write_text(json.dumps(SAMPLE), encoding="utf-8")

            result = subprocess.run(
                [sys.executable, "-m", "stateid", str(path), "--energy-window-ev", "0.5"],
                capture_output=True,
                text=True,
                check=True,
            )

        payload = json.loads(result.stdout)
        self.assertIn("identified_states", payload)
        self.assertEqual(len(payload["identified_states"]), 2)


if __name__ == "__main__":
    unittest.main()
