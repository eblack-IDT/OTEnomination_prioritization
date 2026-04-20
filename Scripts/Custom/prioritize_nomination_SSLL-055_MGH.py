# -*- coding: utf-8 -*-
"""
prioritize_nomination_SSLL-055_MGH.py

Multi-gRNA UNCOVERseq prioritization — all input files carry equal weight.

Naming convention
-----------------
    {gRNA}_Tier{N}_{X}
        gRNA  – label supplied via --grna_names (one per file, same order as --uncoverseq)
        N     – Tier number from the source file
        X     – 1-based sequential off-target index across all off-target rows in that
                        input file after LD < 7 filtering, matching the base pipeline convention
  ON-target rows are named {gRNA}_OnT.
  Sites detected by more than one gRNA (same chr/start/stop/strand) are merged:
    Name  = {gRNA1}_...__{gRNA2}_... (capped at 50 chars → coord fallback)
    Tier  = min across sources
    LD    = min across sources
    Repro = max across sources
    Avg % On-Target UMI = max across sources

Priority passes  (applied after Levenshtein Distance < 7 filter)
--------------------------------------------------------------
  Pass 0: ON-target sites                                (always included, no LD restriction)
  Pass 1: Tier ∈ {1, 2}                                 (only if total selected < MAX_SITES)
  Pass 2: Tier = 3  AND  Reproducibility > 2  AND  Avg. % On-Target UMI > 0.5
                                                         (only if total selected < MAX_SITES)
  Pass 3: Tier = 3  AND  Reproducibility > 1  AND  Avg. % On-Target UMI > 0.5
                                                         (only if total selected < MAX_SITES)
  Pass 4: all remaining LD < 7 sites                    (only if total selected < MAX_SITES)

Within each pass sites are ordered by: LD ↑  →  original row-order ↑
"Row order" is the minimum row-index of the site across all source files
(no file-label bias — the gRNA list order on the command line does not
determine priority).

Proximity merging of adjacent targets applies as in the base pipeline.

Key outputs (identical structure to prioritize_nomination_sites.py):
  *_prioritized.bed
  *_excluded.bed
  *_all_candidates_full.bed
  *_merge_mapping.tsv  (component_names column uses {gRNA}_Tier{N}_{X} tokens)
  *_dedup_nonmerged_prioritized_sites.bed
  *_prioritized_assay_mapping.csv
  *_excluded_assay_mapping.csv
  *_all_candidates_assay_mapping.csv

Usage
-----
  python prioritize_nomination_SSLL-055_MGH.py \\
      --uncoverseq file1.csv,file2.csv,file3.csv \\
      --grna_names  sgA,sgB,sgC \\
      [--out_prefix results/MySample | --outdir results/] \\
      [--max_sites 200] [--merge_targets proximity] [--verbose]
"""

from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MAX_SITES = 200
MAX_MERGED_NAME_LEN = 50


# ---------------------------------------------------------------------------
# Shared helpers (inlined to keep this script independent)
# ---------------------------------------------------------------------------
def _read_table_auto(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    df.columns = [c.encode("utf-8-sig").decode("utf-8-sig").strip() for c in df.columns]
    return df


def _canonicalize_columns(df: pd.DataFrame, rename_map: dict) -> pd.DataFrame:
    token_to_original = {}
    for c in df.columns:
        token = "".join(c.lower().split())
        token_to_original[token] = c
    for token, canonical in rename_map.items():
        if token in token_to_original:
            orig = token_to_original[token]
            if orig != canonical:
                df.rename(columns={orig: canonical}, inplace=True)
    return df


def _safe_int(v, default=None):
    try:
        if pd.isna(v):
            return default
        return int(v)
    except Exception:
        try:
            return int(float(v))
        except Exception:
            return default


def _safe_float(v, default=None):
    try:
        if pd.isna(v):
            return default
        return float(v)
    except Exception:
        return default


def _ensure_str_col(df: pd.DataFrame, col: str, default: str = "") -> None:
    if col in df.columns:
        df[col] = df[col].astype(str)
    else:
        df[col] = pd.Series([default] * len(df), index=df.index)


def _parse_repro(v) -> Optional[int]:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    try:
        n = int(str(v).strip())
        return n if 1 <= n <= 3 else None
    except Exception:
        return None


def _build_location(chr_: str, start: int, stop: int, strand: str) -> str:
    s = (strand or ".")
    try:
        return f"{chr_}:{int(start)}-{int(stop)}({s})"
    except Exception:
        return f"{chr_}:{start}-{stop}({s})"


@dataclass
class Component:
    source: str
    site_type: str
    raw_name: str
    cust_name: Optional[str] = None
    cust_nonon_index: Optional[int] = None
    cust_input_order: Optional[int] = None
    uncov_tier: Optional[int] = None
    uncov_nonon_index: Optional[int] = None
    uncov_input_order: Optional[int] = None
    uncov_ld: Optional[float] = None
    uncov_repro: Optional[int] = None
    ins_index: Optional[int] = None
    ins_ld: Optional[float] = None
    ins_tier: Optional[int] = None
    ins_mar_rank: Optional[int] = None
    gam_region: Optional[str] = None
    orig_chr: Optional[str] = None
    orig_start: Optional[int] = None
    orig_stop: Optional[int] = None
    orig_strand: Optional[str] = None
    name: Optional[str] = None
    location: Optional[str] = None
    lev_dist: Optional[str] = None
    gene: Optional[str] = None
    annot: Optional[str] = None
    pool: Optional[str] = None
    mar_str: Optional[str] = None
    tier_str: Optional[str] = None
    avg_pct_on_target_umi: Optional[str] = None

    def display_token(self) -> str:
        s = self.source
        if s == "customer":
            if self.cust_name:
                return f"Customer_{self.cust_name}"
            if self.cust_nonon_index:
                return f"Customer_OT{self.cust_nonon_index}"
            return "Customer"
        if s == "UNCOVERseq":
            y = self.uncov_tier if self.uncov_tier is not None else 9
            x = self.uncov_nonon_index if (self.uncov_nonon_index is not None and self.uncov_nonon_index > 0) else 1
            return f"Tier{y}_{x}"
        if s == "insilico_pop":
            return f"Insilico_pop{self.ins_index or 1}"
        if s == "insilico_ref":
            return f"Insilico_ref{self.ins_index or 1}"
        return s


@dataclass
class Record:
    chr: str
    start: int
    stop: int
    strand: str
    name_tokens: List[str]
    name_list: List[str]
    sources_present: Set[str]
    components: List[Component]
    is_on: bool = False
    uncov_present: bool = False
    uncov_is_on: bool = False
    uncov_tier: Optional[int] = None
    uncov_nonon_index: Optional[int] = None
    uncov_min_order: Optional[int] = None
    uncov_ld: Optional[float] = None
    uncov_repro: Optional[int] = None
    customer_present: bool = False
    customer_is_on: bool = False
    customer_nonon_index: Optional[int] = None
    customer_name: Optional[str] = None
    cust_min_order: Optional[int] = None
    ins_pop_present: bool = False
    ins_pop_index: Optional[int] = None
    ins_ref_present: bool = False
    ins_ref_index: Optional[int] = None
    ins_ld_min: Optional[float] = None
    ins_tier_min: Optional[int] = None
    ins_mar_rank_min: Optional[int] = None
    ins_ld3_exon: bool = False

    def bed_row(self) -> List[str]:
        name = "_".join([t for t in self.name_tokens if t])
        return [self.chr, str(self.start), str(self.stop), name, "0", (self.strand or ".")]


def _final_name_from_record(r: Record) -> str:
    return "_".join([t for t in r.name_tokens if t])


def _ordered_raw_names_from_components(comps: List[Component]) -> List[str]:
    ordered = sorted(
        comps,
        key=lambda c: (
            0 if c.source == "UNCOVERseq" else 1,
            c.uncov_input_order if c.uncov_input_order is not None else 10 ** 9,
            c.orig_chr or "",
            c.orig_start if c.orig_start is not None else 10 ** 9,
            c.orig_stop if c.orig_stop is not None else 10 ** 9,
            c.raw_name or "",
        ),
    )
    seen_parts: Set[str] = set()
    parts: List[str] = []
    for c in ordered:
        n = (c.raw_name or "").strip()
        if n and n not in seen_parts:
            parts.append(n)
            seen_parts.add(n)
    return parts


def _joined_name_from_components(comps: List[Component], chr_: str, start: int, stop: int, name_delim: str) -> Tuple[str, List[str]]:
    parts = _ordered_raw_names_from_components(comps)
    joined = name_delim.join(parts) if parts else "unknown"
    if len(joined) > MAX_MERGED_NAME_LEN:
        joined = f"{chr_}:{start}-{stop}"
    return joined, parts


def _recompute_aggregates_from_components(r: Record) -> None:
    comps = r.components or []
    r.sources_present = set(c.source for c in comps)
    r.is_on = any(bool(re.search(r"\bon\b", (c.site_type or "").lower())) for c in comps)
    r.uncov_present = any(c.source == "UNCOVERseq" for c in comps)
    r.uncov_is_on = any((c.source == "UNCOVERseq") and bool(re.search(r"\bon\b", (c.site_type or "").lower())) for c in comps)
    uncov_tiers = [c.uncov_tier for c in comps if c.uncov_tier is not None]
    r.uncov_tier = min(uncov_tiers) if uncov_tiers else None
    uncov_nonon = [c.uncov_nonon_index for c in comps if c.uncov_nonon_index]
    r.uncov_nonon_index = min(uncov_nonon) if uncov_nonon else None
    uncov_orders = [c.uncov_input_order for c in comps if c.uncov_input_order is not None]
    r.uncov_min_order = min(uncov_orders) if uncov_orders else None
    uncov_ld = [c.uncov_ld for c in comps if c.uncov_ld is not None]
    r.uncov_ld = min(uncov_ld) if uncov_ld else None
    uncov_repro = [c.uncov_repro for c in comps if c.uncov_repro is not None]
    r.uncov_repro = max(uncov_repro) if uncov_repro else None
    r.customer_present = any(c.source == "customer" for c in comps)
    r.customer_is_on = any((c.source == "customer") and bool(re.search(r"\bon\b", (c.site_type or "").lower())) for c in comps)
    cust_nonon = [c.cust_nonon_index for c in comps if c.cust_nonon_index]
    r.customer_nonon_index = min(cust_nonon) if cust_nonon else None
    cust_orders = [c.cust_input_order for c in comps if c.cust_input_order is not None]
    r.cust_min_order = min(cust_orders) if cust_orders else None
    cust_names = [c.cust_name for c in comps if c.cust_name]
    r.customer_name = cust_names[0] if cust_names else None
    r.ins_pop_present = any(c.source == "insilico_pop" for c in comps)
    r.ins_ref_present = any(c.source == "insilico_ref" for c in comps)
    pop_idx = [c.ins_index for c in comps if c.source == "insilico_pop" and c.ins_index is not None]
    ref_idx = [c.ins_index for c in comps if c.source == "insilico_ref" and c.ins_index is not None]
    r.ins_pop_index = min(pop_idx) if pop_idx else None
    r.ins_ref_index = min(ref_idx) if ref_idx else None
    ins_ld = [c.ins_ld for c in comps if c.ins_ld is not None]
    r.ins_ld_min = min(ins_ld) if ins_ld else None
    ins_tier = [c.ins_tier for c in comps if c.ins_tier is not None]
    r.ins_tier_min = min(ins_tier) if ins_tier else None
    ins_mar = [c.ins_mar_rank for c in comps if c.ins_mar_rank is not None]
    r.ins_mar_rank_min = min(ins_mar) if ins_mar else None
    r.ins_ld3_exon = any((c.ins_ld == 3) and (str(getattr(c, "gam_region", "")).strip().lower() == "exon") for c in comps)


def _merge_two_records(prev: Record, curr: Record, name_delim: str) -> Record:
    chr_ = prev.chr
    start = min(prev.start, curr.start)
    has_customer = prev.customer_present or curr.customer_present
    has_insilico = (prev.ins_pop_present or prev.ins_ref_present or curr.ins_pop_present or curr.ins_ref_present)
    if prev.start == curr.start and has_customer and has_insilico:
        stop = min(prev.stop, curr.stop)
    else:
        stop = curr.stop
    strand = prev.strand if prev.strand == curr.strand else "."
    merged_components = prev.components + curr.components
    merged_name, merged_name_list = _joined_name_from_components(
        merged_components, chr_, start, stop, name_delim
    )
    merged = Record(chr_, start, stop, strand,
                    name_tokens=[merged_name],
                    name_list=merged_name_list,
                    sources_present=set(),
                    components=merged_components)
    _recompute_aggregates_from_components(merged)
    return merged


def merge_adjacent_records(all_records: List[Record],
                           target_distance: int,
                           max_insert: int = 200,
                           target_flank: int = 40,
                           max_merge_size: Optional[int] = None,
                           name_delim: str = "__") -> List[Record]:
    if not all_records:
        return all_records
    if max_merge_size is None:
        max_merge_size = max_insert - target_flank - 15
    recs = sorted(all_records, key=lambda r: (r.chr, r.start, r.stop))
    out: List[Record] = []
    prev: Optional[Record] = None
    for curr in recs:
        if prev is None:
            prev = curr
            continue
        if curr.chr == prev.chr:
            gap = int(curr.start) - int(prev.stop)
            merged_len = int(curr.stop) - int(prev.start)
            if gap < target_distance and merged_len <= max_merge_size:
                prev = _merge_two_records(prev, curr, name_delim)
                continue
            out.append(prev)
            prev = curr
        else:
            out.append(prev)
            prev = curr
    if prev is not None:
        out.append(prev)
    return out


def merge_duplicate_coords(all_records: List[Record], name_delim: str = ";") -> List[Record]:
    if not all_records:
        return all_records
    from collections import defaultdict
    groups: Dict[Tuple[str, int, int, str], List[Record]] = defaultdict(list)
    for r in all_records:
        groups[(r.chr, r.start, r.stop, r.strand or ".")].append(r)
    out: List[Record] = []
    for (chr_, start, stop, strand), items in groups.items():
        if len(items) == 1:
            out.append(items[0])
            continue
        base = items[0]
        for nxt in items[1:]:
            base = _merge_two_records(base, nxt, name_delim)
            base.chr, base.start, base.stop, base.strand = chr_, start, stop, strand
        out.append(base)
    return out


def write_bed(path: str, recs: List[Record]) -> int:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    count = 0
    with open(path, "w") as f:
        for r in recs:
            f.write("\t".join(r.bed_row()) + "\n")
            count += 1
    return count


def build_dedup_nonmerged_prioritized_rows(prioritized: list) -> list:
    rows = []
    seen = set()

    def comp_is_on(c):
        st = (c.site_type or "").lower()
        return bool(re.search(r"\bon\b", st))

    for r in prioritized:
        for c in (r.components or []):
            chr_ = c.orig_chr or r.chr
            start = c.orig_start if c.orig_start is not None else r.start
            stop = c.orig_stop if c.orig_stop is not None else r.stop
            strand = (c.orig_strand or r.strand or ".")

            key = (chr_, int(start), int(stop), strand)
            if key in seen:
                continue
            seen.add(key)

            if comp_is_on(c):
                orig_name = "OnT"
            elif c.source == "UNCOVERseq":
                t = c.uncov_tier if c.uncov_tier is not None else 9
                x = c.uncov_nonon_index if (c.uncov_nonon_index is not None and c.uncov_nonon_index > 0) else 1
                orig_name = f"Tier{t}_{x}"
            elif c.source == "customer":
                if c.cust_name and str(c.cust_name).strip() != "":
                    orig_name = f"Customer_{c.cust_name}"
                else:
                    x = c.cust_nonon_index if c.cust_nonon_index is not None else 1
                    orig_name = f"Customer_OT{x}"
            elif c.source == "insilico_pop":
                orig_name = f"Insilico_pop{c.ins_index or 1}"
            elif c.source == "insilico_ref":
                orig_name = f"Insilico_ref{c.ins_index or 1}"
            else:
                try:
                    orig_name = c.display_token()
                except Exception:
                    orig_name = c.raw_name or "orig"

            rows.append([chr_, str(int(start)), str(int(stop)), orig_name, "0", strand])

    return rows


def write_bed6(path: str, rows: list) -> int:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    n = 0
    with open(path, "w") as f:
        for r in rows:
            line = "\t".join([str(r[0]), str(r[1]), str(r[2]), str(r[3]), str(r[4]), str(r[5])])
            f.write(line + "\n")
            n += 1
    return n


def _semi_join(values: List[Optional[str]]) -> str:
    normed = []
    for v in values:
        s = "" if v is None else str(v)
        s = s.strip()
        normed.append("NONE" if s == "" else s)
    uniq = set(normed)
    if len(uniq) == 1:
        return next(iter(uniq))
    return ";".join(normed)


def _origin_display(source: str) -> str:
    if source.startswith("insilico"):
        return "insilico"
    return source


def _component_sort_key(c: Component) -> int:
    s = c.source or ""
    if s == "UNCOVERseq":
        return 0
    if s.startswith("insilico"):
        return 1
    if s == "customer":
        return 2
    return 9


def _tier_final_from_components(comps: List[Component]) -> str:
    vals = []
    for c in comps:
        if c.tier_str is None:
            continue
        s = str(c.tier_str).strip()
        try:
            n = int(s)
            vals.append(n)
        except Exception:
            try:
                n = int(s.split()[-1])
                vals.append(n)
            except Exception:
                continue
    if not vals:
        return "NONE"
    return str(min(vals))


def write_assay_mapping(path: str, recs: List[Record]) -> int:
    import csv

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    headers = [
        "chr", "start", "stop", "strand",
        "Name", "Location", "Lev. Dist.", "Gene", "Annot.", "Pool",
        "MAR", "Tier", "Tier_Final", "Site Type", "Origin", "Avg. % On-Target UMI",
        "orig_chr", "orig_start", "orig_stop", "orig_strand", "orig_Name",
    ]
    n = 0
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in recs:
            comps = sorted(r.components, key=_component_sort_key)
            final_name = "_".join([t for t in r.name_tokens if t])
            name_uncapped = ";".join([s for s in (r.name_list or []) if s is not None and str(s).strip() != ""])

            locs = [c.location for c in comps]
            levd = [c.lev_dist for c in comps]
            genes = [c.gene for c in comps]
            annots = [c.annot for c in comps]
            pools = [c.pool for c in comps]
            mars = [c.mar_str for c in comps]
            tiers = [c.tier_str for c in comps]
            stypes = [c.site_type for c in comps]
            origins = [_origin_display(c.source or "") for c in comps]
            avg_pct = [c.avg_pct_on_target_umi for c in comps]
            o_chr = [c.orig_chr for c in comps]
            o_start = [str(c.orig_start) if c.orig_start is not None else None for c in comps]
            o_stop = [str(c.orig_stop) if c.orig_stop is not None else None for c in comps]
            o_strand = [c.orig_strand for c in comps]

            row = [
                r.chr, r.start, r.stop, (r.strand or "."),
                final_name,
                _semi_join(locs),
                _semi_join(levd),
                _semi_join(genes),
                _semi_join(annots),
                _semi_join(pools),
                _semi_join(mars),
                _semi_join(tiers),
                _tier_final_from_components(comps),
                _semi_join(stypes),
                _semi_join(origins),
                _semi_join(avg_pct),
                _semi_join(o_chr),
                _semi_join(o_start),
                _semi_join(o_stop),
                _semi_join(o_strand),
                name_uncapped if name_uncapped else final_name,
            ]
            w.writerow(row)
            n += 1
    return n


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _str_or_none(v) -> Optional[str]:
    if v is None:
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    s = str(v).strip()
    return s if s else None


def _get_max_avg_umi(r: Record) -> Optional[float]:
    """Return the maximum Avg. % On-Target UMI across all components (numeric)."""
    vals: List[float] = []
    for c in r.components:
        v = c.avg_pct_on_target_umi
        if v is not None:
            try:
                vals.append(float(str(v).strip()))
            except (ValueError, TypeError):
                pass
    return max(vals) if vals else None


# ---------------------------------------------------------------------------
# Per-file loader
# ---------------------------------------------------------------------------

def _load_one_uncov(path: str, grna_name: str) -> pd.DataFrame:
    """
    Load one UNCOVERseq CSV.
    Assigns canonical site names {gRNA}_Tier{N}_{X} (or {gRNA}_OnT).
    All parent column-normalisation conventions are preserved.
    """
    df = _read_table_auto(path)
    df["__input_order__"] = range(len(df))

    df = _canonicalize_columns(df, {
        "chr":                  "chr",
        "start":                "start",
        "stop":                 "stop",
        "strand":               "strand",
        "tier":                 "Tier",
        "reproducibility":      "Reproducibility",
        "sitetype":             "Site Type",
        "levenshteindistance":  "Levenshtein Distance",
        "annotation":           "Annotation",
        "location":             "Location",
        "gene":                 "Gene",
        "mar":                  "MAR",
        "avg.%on-targetumi":    "Avg. % On-Target UMI",
    })

    for col in ["chr", "start", "stop", "strand"]:
        if col not in df.columns:
            raise ValueError(
                f"[{grna_name}] missing required column '{col}' after normalisation"
            )

    df["chr"]    = df["chr"].astype(str).str.strip()
    df["start"]  = df["start"].apply(_safe_int)
    df["stop"]   = df["stop"].apply(_safe_int)
    df["strand"] = df["strand"].astype(str).str.strip()

    _ensure_str_col(df, "Site Type", "")

    if "Tier" in df.columns:
        df["Tier_num"] = df["Tier"].apply(_safe_int)
    else:
        df["Tier_num"] = None
        df["Tier"] = ""

    if "Reproducibility" not in df.columns:
        df["Reproducibility"] = None
    df["uncov_repro"] = df["Reproducibility"].apply(_parse_repro)

    if "Levenshtein Distance" in df.columns:
        df["UNCOV_ld"] = df["Levenshtein Distance"].apply(_safe_float)
    else:
        df["UNCOV_ld"] = None

    # Hard filter UNCOVERseq to LD < 7 so LD>=7/missing rows never appear downstream.
    # Keep the original surviving file order and use that order for naming/prioritization.
    pass_ld = df["UNCOV_ld"].notna() & (df["UNCOV_ld"] < 7)
    df.attrs["uncov_ld_lt7_count"] = int(pass_ld.sum())
    df.attrs["uncov_ld_ge7_or_missing_count"] = int((~pass_ld).sum())
    df = df.loc[pass_ld].copy().reset_index(drop=True)
    df["__input_order__"] = range(1, len(df) + 1)

    if "Avg. % On-Target UMI" not in df.columns:
        df["Avg. % On-Target UMI"] = None

    # X index: 1-based count across all off-target rows in LD-filtered file order,
    # matching the base pipeline UNCOVERseq convention on a per-input-file basis.
    is_on = df["Site Type"].str.contains(r"\bon\b", case=False, na=False)
    df["UNCOV_nonon_index"] = 0
    nonon_idx = 1
    for i in df.index:
        if is_on.at[i]:
            continue
        df.at[i, "UNCOV_nonon_index"] = nonon_idx
        nonon_idx += 1

    # Canonical site name
    def _site_name(row) -> str:
        if bool(re.search(r"\bon\b", str(row.get("Site Type", "")).lower())):
            return f"{grna_name}_OnT"
        tier = _safe_int(row.get("Tier_num"))
        n = tier if tier is not None else 9
        x_raw = _safe_int(row.get("UNCOV_nonon_index"), 1)
        x = x_raw if (x_raw is not None and x_raw > 0) else 1
        return f"{grna_name}_Tier{n}_{x}"

    df["site_name"]  = df.apply(_site_name, axis=1)
    df["raw_name"]   = df["site_name"]
    df["source"]     = "UNCOVERseq"  # kept as "UNCOVERseq" for parent function compatibility
    df["grna_label"] = grna_name
    df["site_type"]  = df["Site Type"]

    return df


# ---------------------------------------------------------------------------
# Record builder  (multi-gRNA, UNCOVERseq-only)
# ---------------------------------------------------------------------------

def _build_records(
    dfs: List[pd.DataFrame],
    dedupe_mode: str,
) -> Tuple[List[Record], List[Record]]:
    """
    Build Record objects from multiple gRNA-labelled DataFrames.

    Deduplication: sites sharing (chr, start, stop, strand) across files are
    collapsed into one Record; best Tier/LD, max Repro/UMI, all source tokens
    retained.

    Naming: every component stores raw_name = {gRNA}_Tier{N}_{X} (or {gRNA}_OnT).
    The Record's name_tokens are derived directly from component raw_names, joined
    with '__', capped at MAX_MERGED_NAME_LEN.
    """
    all_map: Dict = {}
    _ctr = [0]  # mutable counter for dedupe_mode == "none" unique keys

    for df in dfs:
        for _, row in df.iterrows():
            chr_   = str(row["chr"])
            start  = int(row["start"])
            stop   = int(row["stop"])
            strand = str(row["strand"])

            site_type_val = str(row.get("site_type", ""))
            is_on = bool(re.search(r"\bon\b", site_type_val.lower()))

            if dedupe_mode == "none":
                _ctr[0] += 1
                rec_key: tuple = (chr_, start, stop, strand, _ctr[0])
            else:
                rec_key = (chr_, start, stop, strand)

            rec = all_map.get(rec_key)
            if rec is None:
                rec = Record(
                    chr_, start, stop, strand,
                    name_tokens=[],
                    name_list=[],
                    sources_present=set(),
                    components=[],
                )
                all_map[rec_key] = rec

            rec.sources_present.add("UNCOVERseq")
            rec.is_on         = rec.is_on or is_on
            rec.uncov_present = True
            if is_on:
                rec.uncov_is_on = True

            # ── Build Component ──────────────────────────────────────────
            comp = Component(
                source="UNCOVERseq",
                site_type=site_type_val,
                raw_name=str(row.get("raw_name", "")),
            )
            comp.orig_chr    = chr_
            comp.orig_start  = start
            comp.orig_stop   = stop
            comp.orig_strand = strand

            # Assay-mapping pass-through fields
            comp.name     = str(row.get("site_name", ""))
            loc_raw = row.get("Location")
            comp.location = _str_or_none(loc_raw) or _build_location(chr_, start, stop, strand)
            ld_disp = row.get("Levenshtein Distance", row.get("UNCOV_ld"))
            comp.lev_dist = _str_or_none(ld_disp)
            comp.gene     = _str_or_none(row.get("Gene"))
            comp.annot    = _str_or_none(row.get("Annotation", row.get("Annot.")))
            comp.pool     = _str_or_none(row.get("Pool"))
            comp.mar_str  = _str_or_none(row.get("MAR"))
            comp.tier_str = _str_or_none(row.get("Tier", row.get("Tier_num")))
            comp.avg_pct_on_target_umi = _str_or_none(row.get("Avg. % On-Target UMI"))

            # ── Aggregate Record fields ─────────────────────────────────
            tier = _safe_int(row.get("Tier_num"))
            if tier is not None:
                rec.uncov_tier = (
                    min(rec.uncov_tier, tier) if rec.uncov_tier is not None else tier
                )
            comp.uncov_tier = tier

            ld = _safe_float(row.get("UNCOV_ld"))
            if ld is not None:
                rec.uncov_ld = (
                    min(rec.uncov_ld, ld) if rec.uncov_ld is not None else ld
                )
            comp.uncov_ld = ld

            nonon_idx = _safe_int(row.get("UNCOV_nonon_index"))
            if nonon_idx:
                rec.uncov_nonon_index = (
                    min(rec.uncov_nonon_index, nonon_idx) if rec.uncov_nonon_index is not None
                    else nonon_idx
                )
            comp.uncov_nonon_index = nonon_idx

            in_ord = _safe_int(row.get("__input_order__"))
            if in_ord is not None:
                rec.uncov_min_order = (
                    min(rec.uncov_min_order, in_ord) if rec.uncov_min_order is not None
                    else in_ord
                )
            comp.uncov_input_order = in_ord

            repro = _parse_repro(row.get("uncov_repro", row.get("Reproducibility")))
            if repro is not None:
                rec.uncov_repro = (
                    max(rec.uncov_repro, repro) if rec.uncov_repro is not None else repro
                )
            comp.uncov_repro = repro

            rec.components.append(comp)

    all_records = list(all_map.values())

    # ── Naming pass ──────────────────────────────────────────────────────
    # Use raw_name from each component ({gRNA}_Tier{N}_{X} or {gRNA}_OnT).
    # Multi-source sites join tokens with '__'.
    for r in all_records:
        parts = _ordered_raw_names_from_components(r.components)
        joined = "__".join(parts) if parts else ("OnT" if r.is_on else "unknown")
        if len(joined) > MAX_MERGED_NAME_LEN:
            joined = f"{r.chr}:{r.start}-{r.stop}"
        r.name_tokens = [joined]
        r.name_list   = parts if parts else [joined]

    on_records = [r for r in all_records if r.is_on]
    return all_records, on_records


# ---------------------------------------------------------------------------
# Prioritization — three-pass
# ---------------------------------------------------------------------------

def prioritize_ssll(
    all_recs: List[Record],
    max_sites: int,
) -> Tuple[List[Record], List[Record]]:
    """
    Pass 0: ON-target   (no LD restriction, always included)
    Pass 1: Tier ∈ {1, 2}                      | LD < 7  (filtered file order)
    Pass 2: Tier=3 & Repro>2 & Avg%UMI>0.5    | LD < 7  (filtered file order)
    Pass 3: Tier=3 & Repro>1 & Avg%UMI>0.5    | LD < 7  (filtered file order)
    Pass 4: all remaining LD < 7 sites         |         (filtered file order)
    """
    selected: List[Record] = []
    seen: Set[tuple] = set()

    def _key(r: Record) -> tuple:
        return (r.chr, r.start, r.stop, r.strand)

    def _input_order_key(r: Record) -> tuple:
        order = r.uncov_min_order if r.uncov_min_order is not None else 10 ** 9
        return (order, r.chr, r.start, r.stop, r.strand)

    def add_pass(pool: List[Record]) -> None:
        for r in sorted(pool, key=_input_order_key):
            if len(selected) >= max_sites:
                break
            k = _key(r)
            if k in seen:
                continue
            selected.append(r)
            seen.add(k)

    # Pass 0: ON-target — always include
    add_pass([r for r in all_recs if r.is_on])

    # Remaining passes require LD < 7
    ld_ok = [
        r for r in all_recs
        if not r.is_on and r.uncov_ld is not None and r.uncov_ld < 7
    ]

    # Pass 1: Tier 1 or 2
    if len(selected) < max_sites:
        add_pass([r for r in ld_ok if r.uncov_tier is not None and r.uncov_tier in (1, 2)])

    # Pass 2: Tier 3 & Repro > 2 & Avg% UMI > 0.5
    if len(selected) < max_sites:
        def _pass2(r: Record) -> bool:
            return (
                r.uncov_tier == 3
                and r.uncov_repro is not None and r.uncov_repro > 2
                and (_get_max_avg_umi(r) or 0.0) > 0.5
            )
        add_pass([r for r in ld_ok if _pass2(r)])

    # Pass 3: Tier 3 & Repro > 1 & Avg% UMI > 0.5
    if len(selected) < max_sites:
        def _pass3(r: Record) -> bool:
            return (
                r.uncov_tier == 3
                and r.uncov_repro is not None and r.uncov_repro > 1
                and (_get_max_avg_umi(r) or 0.0) > 0.5
            )
        add_pass([r for r in ld_ok if _pass3(r)])

    # Pass 4: all remaining LD < 7 sites (any tier / repro / UMI) — fill to max_sites
    if len(selected) < max_sites:
        add_pass(ld_ok)

    # Build excluded — stable LD/tier/row-order sort, deduped
    seen_keys: Set[tuple] = set(seen)
    all_sorted = sorted(
        all_recs,
        key=lambda r: (
            r.uncov_ld        if r.uncov_ld        is not None else 999.0,
            r.uncov_tier      if r.uncov_tier      is not None else 9,
            r.uncov_min_order if r.uncov_min_order is not None else 10 ** 9,
        ),
    )
    seen_exc: Set[tuple] = set(seen_keys)
    excluded: List[Record] = []
    for r in all_sorted:
        k = _key(r)
        if k in seen_exc:
            continue
        seen_exc.add(k)
        excluded.append(r)

    return selected, excluded


# ---------------------------------------------------------------------------
# Merge-mapping writer  (override: uses raw_name instead of display_token)
# ---------------------------------------------------------------------------

def _write_merge_mapping(
    path: str,
    all_records: List[Record],
    prioritized: List[Record],
) -> None:
    """
    TSV mapping of merged inputs → final output name.
    component_names column uses {gRNA}_Tier{N}_{X} tokens (raw_name).
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    rank: Dict[tuple, int] = {
        (r.chr, r.start, r.stop, r.strand): i
        for i, r in enumerate(prioritized, start=1)
    }
    cols = [
        "final_name", "chr", "start", "stop", "strand",
        "prioritized", "rank",
        "sources_present", "component_sources", "component_names",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in all_records:
            final_name = "_".join([t for t in r.name_tokens if t])
            k = (r.chr, r.start, r.stop, r.strand)
            row_d = {
                "final_name":       final_name,
                "chr":              r.chr,
                "start":            r.start,
                "stop":             r.stop,
                "strand":           r.strand or ".",
                "prioritized":      "Y" if k in rank else "N",
                "rank":             rank.get(k, ""),
                "sources_present":  ";".join(sorted(r.sources_present)),
                "component_sources": ";".join(c.source for c in r.components),
                "component_names":  ";".join(
                    (c.raw_name or c.source) for c in r.components
                ),
            }
            fh.write("\t".join(str(row_d[c]) for c in cols) + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description=(
            "Multi-gRNA UNCOVERseq prioritization (SSLL-055 / MGH). "
            "All input files carry equal weight; tiebreaker: LD ↑ → row-order ↑."
        )
    )
    p.add_argument(
        "--uncoverseq", required=True,
        help="Comma-separated list of UNCOVERseq CSV paths.",
    )
    p.add_argument(
        "--grna_names", required=True,
        help=(
            "Comma-separated gRNA labels, one per file in the same order as "
            "--uncoverseq.  Used as the prefix in {gRNA}_Tier{N}_{X} names."
        ),
    )
    p.add_argument(
        "--max_sites", type=int, default=None,
        help=f"Override per-run site cap (default {MAX_SITES}).",
    )
    p.add_argument(
        "--dedupe", choices=["by_coord", "none"], default="by_coord",
        help="Deduplication strategy (default: by_coord).",
    )
    p.add_argument(
        "--out_prefix", type=str, default=None,
        help="Path + filename prefix for all outputs.",
    )
    p.add_argument(
        "--outdir", type=str, default=None,
        help="Output directory (used only when --out_prefix is not provided).",
    )
    p.add_argument("--verbose", action="store_true")
    # Target-merge controls (mirrors base pipeline)
    p.add_argument(
        "--merge_targets",
        choices=["proximity", "none", "dedupe_coords"],
        default="proximity",
        help="Target merge strategy (default: proximity).",
    )
    p.add_argument("--target_distance", type=int, default=500)
    p.add_argument("--max_insert",      type=int, default=200)
    p.add_argument("--target_flank",    type=int, default=40)
    p.add_argument("--max_merge_size",  type=int, default=None)
    p.add_argument("--merge_name_delim", type=str, default="__")
    args = p.parse_args(argv)

    MAX_SITES_RUNTIME = MAX_SITES if args.max_sites is None else int(args.max_sites)
    if MAX_SITES_RUNTIME <= 0:
        raise SystemExit("--max_sites must be positive")

    # Parse file / label lists
    uncov_paths  = [x.strip() for x in args.uncoverseq.split(",")  if x.strip()]
    grna_labels  = [x.strip() for x in args.grna_names.split(",")  if x.strip()]
    if len(uncov_paths) != len(grna_labels):
        raise SystemExit(
            f"--uncoverseq has {len(uncov_paths)} entries but "
            f"--grna_names has {len(grna_labels)}. They must match."
        )
    for fp in uncov_paths:
        if not os.path.isfile(fp):
            raise SystemExit(f"File not found: {fp}")

    # Output path helpers
    if args.out_prefix:
        out_dir = os.path.dirname(args.out_prefix) or "."
        base    = os.path.basename(args.out_prefix)
    else:
        out_dir = args.outdir or "."
        base    = "SSLL-055_MGH"
    os.makedirs(out_dir, exist_ok=True)

    def _p(suffix: str) -> str:
        return os.path.join(out_dir, f"{base}{suffix}")

    paths = {
        "prioritized":     _p("_prioritized.bed"),
        "excluded":        _p("_excluded.bed"),
        "full":            _p("_all_candidates_full.bed"),
        "mapping":         _p("_merge_mapping.tsv"),
        "dedup_nonmerged": _p("_dedup_nonmerged_prioritized_sites.bed"),
        "am_prior":        _p("_prioritized_assay_mapping.csv"),
        "am_excl":         _p("_excluded_assay_mapping.csv"),
        "am_all":          _p("_all_candidates_assay_mapping.csv"),
    }

    # Load
    dfs = [_load_one_uncov(fp, gn) for fp, gn in zip(uncov_paths, grna_labels)]

    if args.verbose:
        sep = "-" * 68
        print(sep)
        print("SSLL-055 MGH  —  Multi-gRNA UNCOVERseq Prioritization")
        print(sep)
        print(f"  MAX_SITES (runtime): {MAX_SITES_RUNTIME}")
        print(f"  Dedupe mode:         {args.dedupe}")
        print(f"  Merge strategy:      {args.merge_targets}")
        print()
        for fp, gn, df in zip(uncov_paths, grna_labels, dfs):
            n_on  = int(df["Site Type"].str.contains(r"\bon\b", case=False, na=False).sum())
            n_off = len(df) - n_on
            print(f"  [{gn}]  {fp}")
            print(f"          {len(df)} rows  ({n_on} ON, {n_off} OFF)")
        print()
        print("  Priority passes (after LD < 7 filter):")
        print("    Pass 0: ON-target (no LD restriction)")
        print("    Pass 1: Tier ∈ {1, 2}")
        print("    Pass 2: Tier=3 & Repro>2 & Avg%UMI>0.5")
        print("    Pass 3: Tier=3 & Repro>1 & Avg%UMI>0.5  [if under max]")
        print("    Pass 4: all remaining LD<7 sites          [if under max]")
        print()
        print("  Outputs:")
        for label, path in paths.items():
            print(f"    {label:<20} {path}")
        print(sep)

    # Build records
    all_records, _ = _build_records(dfs, args.dedupe)

    # Proximity / coord merge
    if args.merge_targets == "proximity":
        all_records = merge_adjacent_records(
            all_records,
            target_distance=args.target_distance,
            max_insert=args.max_insert,
            target_flank=args.target_flank,
            max_merge_size=args.max_merge_size,
            name_delim=args.merge_name_delim,
        )
    elif args.merge_targets == "dedupe_coords":
        all_records = merge_duplicate_coords(
            all_records, name_delim=args.merge_name_delim
        )

    # Prioritize
    final, excluded = prioritize_ssll(all_records, MAX_SITES_RUNTIME)

    # Write outputs
    n_pri    = write_bed(paths["prioritized"], final)
    n_exc    = write_bed(paths["excluded"],    excluded)
    n_all    = write_bed(paths["full"],        all_records)
    _write_merge_mapping(paths["mapping"], all_records, final)
    rows_nm  = build_dedup_nonmerged_prioritized_rows(final)
    n_nm     = write_bed6(paths["dedup_nonmerged"], rows_nm) or 0
    n_am_p   = write_assay_mapping(paths["am_prior"], final)
    n_am_e   = write_assay_mapping(paths["am_excl"],  excluded)
    n_am_a   = write_assay_mapping(paths["am_all"],   all_records)

    print(f"Wrote {n_pri:>4} lines : {paths['prioritized']}")
    print(f"Wrote {n_exc:>4} lines : {paths['excluded']}")
    print(f"Wrote {n_all:>4} lines : {paths['full']}")
    print(f"Wrote mapping : {paths['mapping']}")
    print(f"Wrote {n_nm:>4} lines : {paths['dedup_nonmerged']}")
    print(f"Wrote {n_am_p:>4} rows  : {paths['am_prior']}")
    print(f"Wrote {n_am_e:>4} rows  : {paths['am_excl']}")
    print(f"Wrote {n_am_a:>4} rows  : {paths['am_all']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
