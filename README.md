# RhG_shift_MIC

MATLAB workflow for nonlinear conductivity in RhG model, including:

- shift current (`sigma_abc`)
- magnetic injection current (`eta_abc`, metric form)

## Included folders

- `+MTB`
- `+tbHFMF`
- `+tbNLC` (new toolbox wrapper for nonlinear conductivity)

## Main scripts

- `BLG_shift_current_MIC.m` (legacy full workflow, now wired to `+tbNLC`)
- `BLG_shift_current_MIC_refactor.m` (simplified and documented workflow)

## Data policy

Large local datasets are intentionally **not tracked** in this repository.
If you need to run with your own data, place files under `data/` locally.

## Run from terminal (no MATLAB desktop)

```bash
matlab -batch "run('/absolute/path/to/BLG_shift_current_MIC_refactor.m')"
```

