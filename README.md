# stateid

Evidence-based electronic-state identification from first-principles calculations, starting with ABACUS.

## Current scope

This repository currently provides a minimal ABACUS-first workflow:

1. Load ABACUS-derived electronic-state evidence from JSON.
2. Identify states near the Fermi level.
3. Emit labels with explicit evidence traces and confidence scores.

## Expected ABACUS evidence format

```json
{
  "fermi_energy_ev": 0.0,
  "bands": [
    {
      "kpoint_index": 0,
      "band_index": 1,
      "spin": "up",
      "energy_ev": -0.12,
      "occupation": 1.0,
      "projections": { "p": 0.76, "s": 0.24 }
    }
  ]
}
```

## Usage

```bash
python -m stateid /path/to/abacus_state_data.json --energy-window-ev 1.0
```

If a directory is provided, the loader will try these file names in order:

- `abacus_state_data.json`
- `stateid_abacus.json`
- `abacus.json`
