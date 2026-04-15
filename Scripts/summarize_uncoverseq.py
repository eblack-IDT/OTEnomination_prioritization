#!/usr/bin/env python3
"""
summarize_uncoverseq.py

Reads one or more UNCOVERseq CSV files (comma-separated list) and writes
a filter-count summary CSV for each file.

Filter steps (applied in order):
  1. Levenshtein Distance < 7  (base filter; all subsequent counts use this subset)
  2. Count per Tier value
  3. Count where Tier in {1, 2, 3}
  4. Count where Tier in {1, 2, 3}  AND  Reproducibility > 1
  5. Count where Tier in {1, 2, 3}  AND  Reproducibility > 2
  6. Count where Tier in {1, 2, 3}  AND  Reproducibility > 2  AND  Avg. % On-Target UMI > 0.5

Usage:
    python summarize_uncoverseq.py file1.csv[,file2.csv,...] [--outdir OUTPUT_DIR]
    python summarize_uncoverseq.py file1.csv,file2.csv --combine combined_summary.csv [--outdir OUTPUT_DIR]
"""

from __future__ import annotations
import argparse
import os
import sys
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_table_auto(path: str) -> pd.DataFrame:
    """Read CSV (auto-detect separator) and strip BOM/whitespace from headers."""
    df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    df.columns = [c.encode("utf-8-sig").decode("utf-8-sig").strip() for c in df.columns]
    return df


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first matching column name (comparison ignores case and spaces)."""
    token_map = {"".join(c.lower().split()): c for c in df.columns}
    for cand in candidates:
        key = "".join(cand.lower().split())
        if key in token_map:
            return token_map[key]
    return None


def _safe_float(v) -> float | None:
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _safe_int(v) -> int | None:
    try:
        return int(v)
    except (ValueError, TypeError):
        try:
            return int(float(v))
        except (ValueError, TypeError):
            return None


# ---------------------------------------------------------------------------
# Core summarisation
# ---------------------------------------------------------------------------

def summarize_file(path: str) -> pd.DataFrame:
    """Build and return a filter-count summary DataFrame for one UNCOVERseq CSV."""
    df = _read_table_auto(path)

    # Locate columns (flexible name matching)
    ld_col = _find_col(df, ["Levenshtein Distance", "levenshtein_distance", "LevenshteinDistance"])
    tier_col = _find_col(df, ["Tier"])
    repro_col = _find_col(df, ["Reproducibility", "reproducibility"])
    ontar_col = _find_col(df, [
        "Avg. % On-Target UMI",
        "Avg. % On-Target UMI ",
        "Avg_%_On-Target_UMI",
        "Avg.%On-TargetUMI",
    ])

    total_input = len(df)

    # ------------------------------------------------------------------
    # Step 1 – filter by Levenshtein Distance < 7
    # ------------------------------------------------------------------
    if ld_col is None:
        print(
            f"WARNING: 'Levenshtein Distance' column not found in {path}; "
            "no LD filter applied.",
            file=sys.stderr,
        )
        df_filt = df.copy()
    else:
        df["_ld"] = df[ld_col].apply(_safe_float)
        df_filt = df[df["_ld"].notna() & (df["_ld"] < 7)].copy()

    total_after_ld = len(df_filt)

    # ------------------------------------------------------------------
    # Pre-compute typed helper columns on the filtered frame
    # ------------------------------------------------------------------
    df_filt = df_filt.copy()

    if tier_col:
        df_filt["_tier"] = df_filt[tier_col].apply(_safe_int)
    else:
        df_filt["_tier"] = pd.NA

    if repro_col:
        df_filt["_repro"] = df_filt[repro_col].apply(_safe_int)
    else:
        df_filt["_repro"] = pd.NA

    if ontar_col:
        df_filt["_ontar"] = df_filt[ontar_col].apply(_safe_float)
    else:
        df_filt["_ontar"] = pd.NA

    # ------------------------------------------------------------------
    # Build summary rows
    # ------------------------------------------------------------------
    rows: list[dict] = []

    rows.append({"Filter": "Total rows in input file", "Count": total_input})
    rows.append({"Filter": "Levenshtein Distance < 7", "Count": total_after_ld})

    # Per-Tier counts (all Tier values present after LD filter)
    if tier_col:
        tier_counts = df_filt["_tier"].value_counts(dropna=False)
        ordered_tier_values = sorted(
            tier_counts.index,
            key=lambda value: float("inf") if pd.isna(value) else int(value),
        )
        for tier_val in ordered_tier_values:
            cnt = tier_counts[tier_val]
            if pd.isna(tier_val):
                label = "Tier (unknown/blank)"
            else:
                label = f"Tier {int(tier_val)}"
            rows.append({"Filter": label, "Count": int(cnt)})
    else:
        rows.append({"Filter": "Tier (column not found)", "Count": "N/A"})

    # Tier 1, 2, or 3 combined
    if tier_col:
        mask_tier123 = df_filt["_tier"].isin([1, 2, 3])
    else:
        mask_tier123 = pd.Series(False, index=df_filt.index)

    rows.append({"Filter": "Tier = 1, 2, or 3", "Count": int(mask_tier123.sum())})

    # Tier 1–3  &  Reproducibility > 1
    if repro_col:
        mask_repro_gt1 = df_filt["_repro"].notna() & (df_filt["_repro"] > 1)
    else:
        mask_repro_gt1 = pd.Series(False, index=df_filt.index)

    rows.append({
        "Filter": "Tier = 1, 2, or 3 & Reproducibility > 1",
        "Count": int((mask_tier123 & mask_repro_gt1).sum()),
    })

    # Tier 1–3  &  Reproducibility > 2
    if repro_col:
        mask_repro_gt2 = df_filt["_repro"].notna() & (df_filt["_repro"] > 2)
    else:
        mask_repro_gt2 = pd.Series(False, index=df_filt.index)

    rows.append({
        "Filter": "Tier = 1, 2, or 3 & Reproducibility > 2",
        "Count": int((mask_tier123 & mask_repro_gt2).sum()),
    })

    # Tier 1–3  &  Reproducibility > 2  &  Avg. % On-Target UMI > 0.5
    if ontar_col:
        mask_ontar_gt05 = df_filt["_ontar"].notna() & (df_filt["_ontar"] > 0.5)
    else:
        mask_ontar_gt05 = pd.Series(False, index=df_filt.index)

    rows.append({
        "Filter": "Tier = 1, 2, or 3 & Reproducibility > 2 & Avg. % On-Target UMI > 0.5",
        "Count": int((mask_tier123 & mask_repro_gt2 & mask_ontar_gt05).sum()),
    })

    # ------------------------------------------------------------------
    # Emit warnings for any missing expected columns
    # ------------------------------------------------------------------
    missing = []
    if ld_col is None:
        missing.append("Levenshtein Distance")
    if tier_col is None:
        missing.append("Tier")
    if repro_col is None:
        missing.append("Reproducibility")
    if ontar_col is None:
        missing.append("Avg. % On-Target UMI")
    if missing:
        print(
            f"WARNING [{os.path.basename(path)}]: column(s) not found – "
            + ", ".join(f"'{c}'" for c in missing)
            + ". Affected counts set to 0 or N/A.",
            file=sys.stderr,
        )

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------

    return pd.DataFrame(rows, columns=["Filter", "Count"])


# ---------------------------------------------------------------------------
# CLI entry-point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize UNCOVERseq CSV files with progressive filter counts. "
            "Pass a comma-separated list of file paths as the first positional argument."
        )
    )
    parser.add_argument(
        "inputs",
        help="Comma-separated list of UNCOVERseq CSV file paths.",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory for output summary CSVs (default: current working directory).",
    )
    parser.add_argument(
        "--combine",
        metavar="FILENAME",
        default=None,
        help=(
            "Write all summaries into a single CSV file instead of one file per input. "
            "The combined file has one shared 'Filter' column and one count column per input file. "
            "The output is placed in --outdir unless FILENAME is an absolute path."
        ),
    )
    args = parser.parse_args()

    files = [f.strip() for f in args.inputs.split(",") if f.strip()]
    if not files:
        parser.error("No input files provided.")

    os.makedirs(args.outdir, exist_ok=True)

    for f in files:
        if not os.path.isfile(f):
            print(f"ERROR: File not found: {f}", file=sys.stderr)
            sys.exit(1)

    if args.combine:
        # ---- combined mode ------------------------------------------------
        combined: pd.DataFrame | None = None
        row_order: list[str] = []
        for f in files:
            summary = summarize_file(f)
            column_name = os.path.splitext(os.path.basename(f))[0]
            summary = summary.rename(columns={"Count": column_name})
            if combined is None:
                combined = summary.copy()
                row_order = combined["Filter"].tolist()
            else:
                for filter_name in summary["Filter"]:
                    if filter_name not in row_order:
                        row_order.append(filter_name)
                combined = combined.merge(summary, on="Filter", how="outer", sort=False)
        combined = combined.set_index("Filter").reindex(row_order).reset_index()
        if os.path.isabs(args.combine):
            out_path = args.combine
        else:
            out_path = os.path.join(args.outdir, args.combine)
        combined.to_csv(out_path, index=False)
        print(f"Wrote combined summary: {out_path}")
    else:
        # ---- per-file mode ------------------------------------------------
        for f in files:
            summary = summarize_file(f)
            basename = os.path.splitext(os.path.basename(f))[0]
            out_name = basename + "_filter_summary.csv"
            out_path = os.path.join(args.outdir, out_name)
            summary.to_csv(out_path, index=False)
            print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
