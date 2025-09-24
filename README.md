# Abundatron: Batch INSPECT Fetcher

Automated tool for fetching NLTE abundance corrections from the  [INSPECT database](https://www.inspect-stars.com).

## What it does

INSPECT provides online calculators to:
- Convert LTE abundances into **3D NLTE abundances**.
- Convert **equivalent widths (mÅ)** into abundances.

This script wraps those calculators into a **batch mode**.  
Give it a list of input values (LTE abundances or equivalent widths) plus stellar parameters,  
and it will query INSPECT for each case and return a clean CSV table.

---

## Features

- **Two modes**:
  - `--mode ew` : Equivalent width → LTE/3D NLTE abundance.
  - `--mode lte`: LTE abundance → 3D NLTE abundance.
- **Flexible input**: pass values as a comma list, from a file, or piped via stdin.
- **Line selection**: choose the spectral line either by wavelength (`--wavelength`) or INSPECT’s internal index (`--wi`).
- **Robust**: retries on network hiccups, checks ranges, throttles requests (`--sleep`).
- **Progress output**: shows `[1/3] … DONE` with the key numbers.
- **Student-friendly**: errors explain what went wrong (e.g. “check if input is within parameter space”).

---

## Installation

Requires **Python 3.8+** and a few common packages:

```bash
pip install requests beautifulsoup4 urllib3
```

Clone this repo, then run:

```bash
python3 abundatron.py --help
```

---

## Usage

### Oxygen example (EW mode)

```bash
python3 abundatron.py   --element O --mode ew --values 65,80,100   --teff 5777 --logg 4.44 --feh 0.0 --vt 1.0   --wavelength 7771.957 --out o7771_from_ew.csv
```

### LTE mode from a file

```bash
python3 abundatron.py   --element O --mode lte --values-file abundances.txt   --teff 5777 --logg 4.44 --feh 0.0 --vt 1.0   --wi 3
```

### Pipe values via stdin

```bash
cat ew_list.txt | python3 abundatron.py   --element O --mode ew   --teff 5777 --logg 4.44 --feh 0.0 --vt 1.0   --wavelength 7774.156
```

---

## Flags

**Required:**
- `--element, -e` : element symbol (e.g. `O`, `Li`, `Na`).
- `--mode` : `"ew"` (EW→abundance) or `"lte"` (LTE→NLTE).
- `--teff` : effective temperature [K].
- `--logg` : surface gravity log g [cgs].
- `--feh` : metallicity [Fe/H].
- `--vt` : microturbulence [km/s].
- **Line choice**: either `--wavelength` (Å) or `--wi` (INSPECT index).

**Input values:**
- `--values` : comma-separated list (EW in mÅ if `--mode ew`, or LTE abundance if `--mode lte`).
- `--values-file` : text/CSV file with one value per line.
- *stdin* : pipe values directly.

**Output & control:**
- `--out` : output CSV (default = stdout).
- `--sleep` : pause (s) between requests, default 0.2.
- `--quiet` : suppress progress prints.
- `--clip` : clip out-of-range inputs to INSPECT’s allowed parameter ranges.

---

## Output

- Progress for each query:
  ```
  [1/3] mode=ew val=65 → A_LTE=8.778 A_NLTE=8.582 Δ=-0.196 [O/Fe]_NLTE=-0.118 … DONE
  ```
- Final CSV written to `--out` (or stdout).

---

## Notes

- INSPECT enforces parameter ranges that depend on element and line (e.g., Teff 5000–6500 K for O I).  
  Out-of-range values will error out with a clear message.
- For big batches, be polite: use `--sleep 0.5` or so.
- Microturbulence ranges differ by element! (e.g. Li requires 1.0–5.0 km/s, O allows 0.5–2.0 km/s).

---

## Citation

If you use INSPECT data in publications, please cite the original sources  
(e.g. [Amarsi et al. (2015) MNRAS, 454,11](https://ui.adsabs.harvard.edu/abs/2015MNRAS.454L..11A/abstract)) as well as Pietrow et al. (in prep).
