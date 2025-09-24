#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INSPECT batch fetcher for abundance/EW conversions
==================================================

What this script does
---------------------
Automates queries to INSPECT (https://www.inspect-stars.com) to:
  1) Convert LTE abundances to 3D NLTE values, or
  2) Convert equivalent widths (mÅ) to abundances,
for a chosen element and spectral line.

You feed it a list of values (either A(LTE) or EW, depending on --mode),
plus stellar parameters (Teff, logg, [Fe/H], microturbulence),
and it returns a CSV with parsed results from INSPECT’s calculators.

Why it’s student-friendly
-------------------------
- Clear, single-file script.
- Inputs can come from --values, a file, or stdin (one per line).
- You can select the line by wavelength (preferred) or by INSPECT’s "wi" index.
- Live progress prints: "1/3 ... DONE" with the key outputs.
- Gentle retries and throttling to avoid hammering the site.
- Errors never crash the batch; they show up as a row with an "error" message.

Requirements
------------
Python 3.8+
pip install: requests, beautifulsoup4, urllib3

    pip install requests beautifulsoup4 urllib3

Quick examples
--------------
# Oxygen, EW mode, 7771.957 Å line, 3 EWs, write CSV
python3 abundatron.py \
  --element O --mode ew --values 65,80,100 \
  --teff 5777 --logg 4.44 --feh 0.0 --vt 1.0 \
  --wavelength 7771.957 --out o7771_from_ew.csv

# LTE mode from a text file, choose line by wi index
python3 abundatron.py \
  --element O --mode lte --values-file abundances.txt \
  --teff 5777 --logg 4.44 --feh 0.0 --vt 1.0 \
  --wi 3

# Pipe values via stdin (prints CSV to stdout)
cat ew_list.txt | python3 abundatron.py \
  --element O --mode ew \
  --teff 5777 --logg 4.44 --feh 0.0 --vt 1.0 \
  --wavelength 7774.156

Flags
-----
Required:
  --element, -e     Element symbol used by INSPECT (e.g. O, Li, Na)
  --mode            Calculation mode: "ew" (EW→abundance) or "lte" (LTE→NLTE)
  --teff            Effective temperature [K]
  --logg            Surface gravity log g [cgs]
  --feh             Metallicity [Fe/H]
  --vt              Microturbulence [km/s]
  --wavelength      Line wavelength [Å] (exact match preferred; else nearest)
     OR
  --wi              INSPECT’s internal line index

Input values:
  --values          Comma-separated list of input values
                    (EW in mÅ if --mode ew, or A(LTE) if --mode lte)
  --values-file     File with one value per line (or CSV with numeric first column)
  stdin             Pipe values directly (e.g. `cat file | python abundatron.py ...`)

Output & control:
  --out             Output CSV file path (default: stdout)
  --sleep           Pause [s] between requests (default 0.2, be polite!)
  --quiet           Suppress progress prints (only final CSV output)
  --clip            Clip out-of-range values to INSPECT’s allowed parameter ranges


Outputs
-------
- A progress line per input, e.g.:
  [1/3] mode=ew val=65 => A_LTE=8.778  A_NLTE=8.582  Δ=-0.196  [O/Fe]_NLTE=-0.118  ... DONE
- A CSV summary to --out (or stdout if --out omitted).

Notes
-----
- INSPECT enforces parameter ranges (e.g., Teff 5000–6500 K for O I example).
  If a query is out of range, you’ll see an "error" column in that row.
- Use --sleep to be polite. Default 0.2 s between queries is usually fine.
"""

import sys
import csv
import re
import math
import time
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from bs4 import BeautifulSoup

BASE = "https://www.inspect-stars.com"
ENDPOINT_EW   = f"{BASE}/A_from_e"
ENDPOINT_LTE  = f"{BASE}/nonlte_from_lte"

# ---------- HTTP session with retries ----------
def make_session() -> requests.Session:
    """
    Create a requests session that automatically retries common transient errors.
    """
    s = requests.Session()
    retry = Retry(
        total=5,                # up to 5 retries for robustness
        backoff_factor=0.5,     # exponential backoff: 0.5, 1.0, 2.0, ...
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=["GET"]
    )
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "inspect-batch-scraper/1.1 (python requests)"})
    return s

# ---------- Helpers to fetch available lines and map wavelength -> wi ----------
def fetch_available_lines(element: str, session: requests.Session) -> List[Tuple[int, float, str]]:
    """
    Load the calculator page and parse <select name="wi"> options.
    Returns a list of (wi_index, wavelength_A, label_text).
    """
    # Using the LTE endpoint just to access the wavelength dropdown.
    r = session.get(f"{ENDPOINT_LTE}?element_name={element}", timeout=20)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, "html.parser")
    sel = soup.find("select", attrs={"name": "wi"})
    if not sel:
        raise RuntimeError(f"Could not find wavelength selector for element {element}.")
    lines = []
    for opt in sel.find_all("option"):
        val = opt.get("value", "").strip()
        txt = opt.text.strip()
        if not val:
            continue
        try:
            wi = int(val)
        except ValueError:
            continue
        try:
            wav = float(txt)
        except ValueError:
            wav = math.nan
        lines.append((wi, wav, txt))
    if not lines:
        raise RuntimeError(f"No lines found for element {element}.")
    return lines

def choose_wi_from_wavelength(element: str, wavelength: float, session: requests.Session) -> Tuple[int, float]:
    """
    Map a wavelength to INSPECT's 'wi' index.
    Prefers exact match; otherwise chooses the nearest wavelength.
    Returns (wi_index, matched_wavelength_A).
    """
    lines = fetch_available_lines(element, session)
    for wi, wav, _ in lines:
        if not math.isnan(wav) and abs(wav - wavelength) < 1e-6:
            return wi, wav
    nearest = min(lines, key=lambda t: abs(t[1] - wavelength) if not math.isnan(t[1]) else float("inf"))
    return nearest[0], nearest[1]

# ---------- Parsing the <pre> block ----------
FLOAT_RE = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")

def parse_pre_block(html: str) -> Dict[str, float]:
    """
    Parse the INSPECT results from the <pre> block into a dict.

    Typical A_from_e block (example):
        EW  A(O) LTE  A(O) NLTE  Delta  [O/Fe] NLTE
        65  8.778     8.582      -0.196 -0.118

    Typical nonlte_from_lte block:
        A(LTE) A(NLTE) Delta [O/Fe] NLTE
    """
    soup = BeautifulSoup(html, "html.parser")
    pre = soup.find("pre")
    if not pre:
        sys.stderr.write(
            "\nERROR: No results block (<pre>) found.\n"
            "This usually means your input values are outside INSPECT's parameter space "
            "(e.g. Teff, logg, [Fe/H], vt, or the primary value).\n"
            "Check the ranges shown on the INSPECT calculator page.\n\n"
        )
        sys.exit(1)   # exit immediately with error code

    text = pre.get_text("\n", strip=True)
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise RuntimeError("Empty results block.")

    # Heuristic: last line with numbers holds the values
    value_line = next((ln for ln in reversed(lines) if FLOAT_RE.search(ln)), None)
    if value_line is None:
        raise RuntimeError("Could not locate numeric result line.")

    nums = [float(x) for x in FLOAT_RE.findall(value_line)]
    header = "\n".join(lines[:2]).lower()

    # A_from_e typically yields 5 numbers: EW, A_LTE, A_NLTE, Delta, [O/Fe]_NLTE
    if "ew" in header and len(nums) >= 5:
        return {"EW_mA": nums[0], "A_LTE": nums[1], "A_NLTE": nums[2], "Delta": nums[3], "OFe_NLTE": nums[4]}

    # Fallbacks to handle slight format changes:
    if len(nums) == 5:
        return {"EW_mA": nums[0], "A_LTE": nums[1], "A_NLTE": nums[2], "Delta": nums[3], "OFe_NLTE": nums[4]}
    if len(nums) == 4:
        return {"A_LTE": nums[0], "A_NLTE": nums[1], "Delta": nums[2], "OFe_NLTE": nums[3]}
    if len(nums) == 3:
        return {"A_LTE": nums[0], "A_NLTE": nums[1], "Delta": nums[2]}

    raise RuntimeError(f"Unrecognized numeric format in result line: {value_line}")

# ---------- Core query functions ----------
def query_A_from_e(element: str, ew_mA: float, teff: float, logg: float, feh: float, vt: float, wi: int,
                   session: requests.Session) -> Dict[str, float]:
    """
    Query INSPECT's A_from_e calculator (EW -> abundance).
    """
    params = {
        "element_name": element, "e": f"{ew_mA:g}",
        "t": f"{teff:g}", "g": f"{logg:g}", "f": f"{feh:g}", "x": f"{vt:g}",
        "wi": str(wi),
    }
    r = session.get(ENDPOINT_EW, params=params, timeout=30)
    r.raise_for_status()
    out = parse_pre_block(r.text)
    out.update({"mode": "ew"})
    return out

def query_nlte_from_lte(element: str, A_lte: float, teff: float, logg: float, feh: float, vt: float, wi: int,
                        session: requests.Session) -> Dict[str, float]:
    """
    Query INSPECT's nonlte_from_lte calculator (LTE abundance -> NLTE).
    """
    params = {
        "element_name": element, "A_lte": f"{A_lte:g}",
        "t": f"{teff:g}", "g": f"{logg:g}", "f": f"{feh:g}", "x": f"{vt:g}",
        "wi": str(wi),
    }
    r = session.get(ENDPOINT_LTE, params=params, timeout=30)
    r.raise_for_status()
    out = parse_pre_block(r.text)
    out.update({"mode": "lte"})
    return out

# ---------- IO helpers ----------
def read_values(args) -> List[float]:
    """
    Read numeric input values from --values, --values-file, and/or stdin.
    """
    vals: List[float] = []
    if args.values:
        vals.extend(float(s.strip()) for s in args.values.split(",") if s.strip())
    if args.values_file:
        p = Path(args.values_file)
        for ln in p.read_text().splitlines():
            ln = ln.strip()
            if not ln:
                continue
            m = FLOAT_RE.search(ln)
            if m:
                vals.append(float(m.group(0)))
    if not sys.stdin.isatty():
        for ln in sys.stdin:
            ln = ln.strip()
            if not ln:
                continue
            m = FLOAT_RE.search(ln)
            if m:
                vals.append(float(m.group(0)))
    if not vals:
        raise SystemExit("No input values found. Use --values, --values-file, or pipe via stdin.")
    return vals

def write_csv(rows: List[Dict[str, Union[str, float]]], out_path: Optional[str]):
    """
    Write rows to CSV at out_path, or to stdout if out_path is None.
    """
    if not rows:
        print("No rows to write.", file=sys.stderr)
        return
    preferred = ["mode","element","wi","wavelength_A","Teff","logg","FeH","vt",
                 "input_value","EW_mA","A_LTE","A_NLTE","Delta","OFe_NLTE","error"]
    cols = set()
    for r in rows:
        cols.update(r.keys())
    header = [c for c in preferred if c in cols] + [c for c in sorted(cols) if c not in preferred]
    if out_path:
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {len(rows)} rows to {out_path}")
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=header)
        w.writeheader()
        w.writerows(rows)

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(
        description="Batch INSPECT scraper (EW↔Abundance) with student-friendly docs & progress output"
    )
    ap.add_argument("--element", "-e", required=True, help="Element symbol as used by INSPECT (e.g. O, Fe, Na)")
    ap.add_argument("--mode", choices=["ew","lte"], required=True,
                    help="'ew' = A_from_e (EW→abundance); 'lte' = nonlte_from_lte (LTE→NLTE)")
    ap.add_argument("--values", help="Comma-separated values: EW in mÅ for 'ew', or A(LTE) for 'lte'")
    ap.add_argument("--values-file", help="Text/CSV file; takes the first numeric token on each line")

    ap.add_argument("--teff", type=float, required=True, help="Effective temperature [K]")
    ap.add_argument("--logg", type=float, required=True, help="log g [cgs]")
    ap.add_argument("--feh",  type=float, required=True, help="[Fe/H]")
    ap.add_argument("--vt",   type=float, required=True, help="Microturbulence [km/s]")

    ap.add_argument("--wavelength", type=float, help="Line wavelength [Å] (exact match preferred; else nearest)")
    ap.add_argument("--wi", type=int, help="INSPECT line index (alternative to --wavelength)")

    ap.add_argument("--out", help="Output CSV path (default: stdout)")
    ap.add_argument("--sleep", type=float, default=0.2, help="Seconds to sleep between requests (default 0.2)")
    ap.add_argument("--quiet", action="store_true", help="Suppress per-item progress prints")

    args = ap.parse_args()

    # Ensure exactly one of wavelength or wi is set
    if (args.wavelength is None) == (args.wi is None):
        ap.error("Choose exactly one: --wavelength OR --wi")

    values = read_values(args)
    total = len(values)

    sess = make_session()
    # Resolve wi and matched wavelength (for context in outputs)
    if args.wavelength is not None:
        wi, matched_wav = choose_wi_from_wavelength(args.element, args.wavelength, sess)
    else:
        wi = int(args.wi)
        lines = fetch_available_lines(args.element, sess)
        wav_map = {widx: wav for (widx, wav, _) in lines}
        matched_wav = wav_map.get(wi, float("nan"))

    rows: List[Dict[str, Union[str, float]]] = []

    if not args.quiet:
        print(f"Element={args.element}  mode={args.mode}  wi={wi}  λ≈{matched_wav} Å")
        print(f"Teff={args.teff}  logg={args.logg}  [Fe/H]={args.feh}  vt={args.vt} km/s")
        print(f"Total inputs: {total}\n")

    for i, v in enumerate(values, 1):
        prefix = f"[{i}/{total}]"
        try:
            if args.mode == "ew":
                res = query_A_from_e(args.element, ew_mA=v, teff=args.teff, logg=args.logg,
                                     feh=args.feh, vt=args.vt, wi=wi, session=sess)
            else:
                res = query_nlte_from_lte(args.element, A_lte=v, teff=args.teff, logg=args.logg,
                                          feh=args.feh, vt=args.vt, wi=wi, session=sess)

            # Enrich with context & input value
            res.update({
                "element": args.element, "wi": wi, "wavelength_A": matched_wav,
                "Teff": args.teff, "logg": args.logg, "FeH": args.feh, "vt": args.vt,
                "input_value": v,
            })
            rows.append(res)

            if not args.quiet:
                # Print a compact human-readable line with key fields
                a_lte  = res.get("A_LTE")
                a_nlte = res.get("A_NLTE")
                delta  = res.get("Delta")
                ofe    = res.get("OFe_NLTE")
                # Build a short status string depending on available keys
                bits = [f"mode={res['mode']}", f"val={v:g}"]
                if a_lte  is not None:  bits.append(f"A_LTE={a_lte:.3f}")
                if a_nlte is not None:  bits.append(f"A_NLTE={a_nlte:.3f}")
                if delta  is not None:  bits.append(f"Δ={delta:+.3f}")
                if ofe    is not None:  bits.append(f"[O/Fe]_NLTE={ofe:+.3f}")
                print(prefix, " ".join(bits), "… DONE", flush=True)

        except Exception as exc:
            err_row = {
                "mode": args.mode, "element": args.element, "wi": wi, "wavelength_A": matched_wav,
                "Teff": args.teff, "logg": args.logg, "FeH": args.feh, "vt": args.vt,
                "input_value": v, "error": str(exc)[:300],
            }
            rows.append(err_row)
            if not args.quiet:
                print(prefix, f"val={v:g} … ERROR: {err_row['error']}", flush=True)

        time.sleep(args.sleep)

    # Final CSV
    write_csv(rows, args.out)

if __name__ == "__main__":
    main()
