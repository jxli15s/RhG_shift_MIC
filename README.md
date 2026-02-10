# RhG_shift_MIC

MATLAB workflow for nonlinear conductivity in RhG model, including:

- shift current (`sigma_abc`)
- magnetic injection current (`eta_abc`, metric form)

Also includes a Python + Numba package (`tbnlc`) for dense k-mesh runs.

## Included folders

- `+MTB`
- `+tbHFMF`
- `+tbNLC` (new toolbox wrapper for nonlinear conductivity)

## Main scripts

- `BLG_shift_current_MIC.m` (legacy full workflow, now wired to `+tbNLC`)
- `BLG_shift_current_MIC_refactor.m` (simplified and documented workflow)

## Python package (`tbnlc`)

Install (editable mode):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Quick import check:

```bash
python -c "import tbnlc; print('tbnlc import ok')"
```

Run U-scan script (minimal output arrays):

```bash
python scripts/hfmf_u_scan.py \
  --hfmf-data-dir /path/to/matlab_data \
  --file-pattern "ne=0.0000e12_U={u:.3f}data.mat" \
  --u-start 1 --u-stop 31 --u-step 1 \
  --band-list 1,2
```

Saved `.npz` fields include:

- `sigma_shift_u`
- `eta_mic_u` (unless `--disable-mic`)
- `Kx`, `Ky`, `Eph_list`, `U_list`, `success_mask`, `Ef_list`

## Data policy

Large local datasets are intentionally **not tracked** in this repository.
If you need to run with your own data, place files under `data/` locally.

## Run from terminal (no MATLAB desktop)

```bash
matlab -batch "run('/absolute/path/to/BLG_shift_current_MIC_refactor.m')"
```
