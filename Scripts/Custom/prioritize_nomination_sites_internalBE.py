# -*- coding: utf-8 -*-
"""
Prioritization pipeline for nomination sites -> rhAmpSeq confirmation.
Implements stepwise logic from "20260218_DraftPrioritizationLogic.docx".

Per requirements:
- CHANGEseq-BE and UNCOVERseq inputs MUST have headers: chr, start, stop, strand (strand required).
- Runtime cap: --max_sites of merged target(default 200).
- De-duplication key: (chr, start, stop, strand) unless --dedupe none.
- Output naming via --out_prefix (path + prefix) or --outdir (dir only).
- --verbose prints an echo banner with resolved paths, settings, and input counts.

Other behavior:
- ON-target naming: If one ONT → OnT_1. If multiple/incongruent ONTs, label by origin
    (OnT_CHANGEseq-BE_*, OnT_UNCOVERseq_*).
- Prioritization order:
    1) Top 60 CHANGEseq-BE by avg CHANGEseq-BE_rep*reads (desc), then Edit_distance (asc)
    2) Top 20 in-silico by existing insilico ranking
    3) Remaining slots from UNCOVERseq (LD < 7, original order)

Key outputs (tab-delimited, 0-based coordinates):
 1) prioritized_sites.bed → chr, start, stop, name, 0, strand
 2) excluded_sites.bed → any candidates not selected due to MaxSites
 3) all_candidates_full.bed → all candidates irrespective of MaxSites
 4) *_merge_mapping.tsv → map of merged components → final name

Additional outputs:
 5) *_prioritized_assay_mapping.csv
 6) *_excluded_assay_mapping.csv
 7) *_all_candidates_assay_mapping.csv

Assay mapping CSV (one row per FINAL target) columns:
 chr, start, stop, strand, Name, Location, Lev. Dist., Gene, Annot., Pool,
 MAR, Tier, Tier_Final, Site Type, Origin, Avg. % On-Target UMI,
 orig_chr, orig_start, orig_stop, orig_strand, orig_Name

Aggregation rules:
 - Values are aggregated across original components using ';'
 - If all values for a column are identical, collapse to that single value
 - Missing/blank values are written as 'NONE' before aggregation
 - Output 'Annot.' is derived from UNCOVERseq 'Annotation', insilico 'gam_region', and 'NONE' for CHANGEseq-BE
 - Component ordering for aggregation: CHANGEseq-BE → insilico → UNCOVERseq
 - Origins 'insilico_pop' and 'insilico_ref' are collapsed to 'insilico'
"""

from __future__ import annotations
import argparse
import os
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set
import pandas as pd
import re

# ------------------------------
# Constants & ranking helpers
# ------------------------------
MAX_SITES = 200
MAR_ORDER = {"High": 0, "Med.": 1, "Med": 1, "Medium": 1, "Low": 2}

# Name length cap for merged names (match bed_qc behavior)
MAX_MERGED_NAME_LEN = 50

# ------------------------------
# Robust file reading & header normalization
# ------------------------------
def _read_table_auto(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    df.columns = [c.encode('utf-8-sig').decode('utf-8-sig').strip() for c in df.columns]
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


def _find_first_col_by_token_prefix(df: pd.DataFrame, token_prefixes: List[str]) -> Optional[str]:
    """Return first column whose normalized token starts with any provided prefix."""
    for col in df.columns:
        token = "".join(str(col).lower().split())
        for prefix in token_prefixes:
            if token.startswith(prefix):
                return col
    return None

# ------------------------------
# Parsers (no coordinate conversions)
# ------------------------------
def parse_pos_field(pos: str) -> Tuple[str, int, int, str]:
    s = str(pos).strip()
    m_strand = re.search(r"\(([+-])\)\s*$", s)
    strand = m_strand.group(1) if m_strand else "."
    s_clean = re.sub(r"\(([+-])\)\s*$", "", s)
    if ":" not in s_clean:
        raise ValueError(f"Invalid pos (missing ':'): {pos}")
    chr_part, rest = s_clean.split(":", 1)
    rest = rest.replace(",", "-")
    if "-" not in rest:
        raise ValueError(f"Invalid pos (missing start-stop): {pos}")
    a, b = rest.split("-", 1)
    start = int(a.strip())
    stop = int(b.strip())
    return chr_part.strip(), start, stop, strand

# ------------------------------
# Safe type helpers & normalizers
# ------------------------------
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

def _norm_mar(v: Optional[str]) -> str:
    t = (v or "").strip()
    if t in MAR_ORDER:
        return "Med." if t in ("Med", "Medium") else t
    return "Low"

def _ensure_str_col(df: pd.DataFrame, col: str, default: str = "") -> None:
    if col in df.columns:
        df[col] = df[col].astype(str)
    else:
        df[col] = pd.Series([default] * len(df), index=df.index)

def _sanitize_name(v: Optional[str]) -> Optional[str]:
    t = (v or "").strip()
    if t == "":
        return None

    # allow double underscores originating from merge
    t = re.sub(r"[^A-Za-z0-9._\\-]+", "_", t)

    # collapse only sequences of 3 or more underscores
    t = re.sub(r"_{3,}", "_", t)

    return t or None

def _parse_repro(v) -> Optional[int]:
    """
    UNCOVERseq Reproducibility is always an integer 1–3.
    Anything else becomes None.
    """
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    try:
        n = int(str(v).strip())
        return n if 1 <= n <= 3 else None
    except:
        return None

def _build_location(chr_: str, start: int, stop: int, strand: str) -> str:
    s = (strand or ".")
    try:
        return f"{chr_}:{int(start)}-{int(stop)}({s})"
    except Exception:
        return f"{chr_}:{start}-{stop}({s})"

def _final_name_from_record(r: "Record") -> str:
    """Return the current final BED name for a record."""
    return "_".join([t for t in r.name_tokens if t])

# ------------------------------
# Loaders (no coordinate conversions)
# ------------------------------
def load_changeseq_be_csv(path: Optional[str]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=[
            "chr","start","stop","strand","site_type","raw_name","cust_name",
            "source","CUST_nonon_index","__cust_input_order__",
            "CSEQ_avg_reads","CSEQ_edit_distance"
        ])

    df = _read_table_auto(path)
    df = _canonicalize_columns(df, {
        "chromosome":"chr", "chr":"chr",
        "start":"start", "end":"stop", "stop":"stop",
        "strand":"strand",
        "sitetype":"Site Type", "name":"Name"
    })

    if "chr" not in df.columns and "Genomic_Coordinate" in df.columns:
        parsed = df["Genomic_Coordinate"].astype(str).str.extract(r"^([^:]+):(\d+)-(\d+)$")
        if parsed.shape[1] == 3:
            df["chr"] = parsed[0]
            df["start"] = parsed[1]
            df["stop"] = parsed[2]

    for col in ["chr","start","stop","strand"]:
        if col not in df.columns:
            raise ValueError(f"CHANGEseq-BE: missing required column '{col}' after normalization")

    df["chr"] = df["chr"].astype(str).str.strip()
    df["start"] = df["start"].apply(_safe_int)
    df["stop"] = df["stop"].apply(_safe_int)
    df["strand"]= df["strand"].astype(str).str.strip()
    df["__cust_input_order__"] = range(1, len(df) + 1)

    rep_cols = [c for c in df.columns if re.match(r"^change(seq)?-be_rep\d+_reads$", "".join(str(c).lower().split()))]
    if not rep_cols:
        rep_cols = [c for c in df.columns if re.match(r"^changeseq-be_rep\d+_reads$", "".join(str(c).lower().split()))]
    if rep_cols:
        rep_vals = df[rep_cols].applymap(_safe_float)
        df["CSEQ_avg_reads"] = rep_vals.mean(axis=1, skipna=True)
    else:
        avg_col = _find_first_col_by_token_prefix(df, ["avgchangeseq-be", "avgchangeseqbe", "avgchange-seq-be", "avgchange-seq"])
        df["CSEQ_avg_reads"] = df[avg_col].apply(_safe_float) if avg_col else pd.NA

    ed_col = _find_first_col_by_token_prefix(df, ["edit_distance", "editdistance"])
    df["CSEQ_edit_distance"] = df[ed_col].apply(_safe_float) if ed_col else pd.NA

    # Highest avg reads first, then lowest edit distance.
    df = df.sort_values(
        by=["CSEQ_avg_reads", "CSEQ_edit_distance", "__cust_input_order__"],
        ascending=[False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    df["__cust_input_order__"] = range(1, len(df) + 1)

    _ensure_str_col(df, "Site Type", "")
    is_on = df["Site Type"].str.contains(r"\bon\b", case=False, na=False)
    df["CUST_nonon_index"] = None
    mask = ~is_on
    df.loc[mask, "CUST_nonon_index"] = range(1, int(mask.sum()) + 1)
    df["site_type"] = df["Site Type"]

    if "Name" in df.columns:
        s = df["Name"].astype(str)
        s_clean = s.where(s.str.strip() != "", None)
        df["cust_name"] = s_clean.fillna("")
        df["raw_name"] = s_clean.fillna("CHANGEseq-BE")
    else:
        df["cust_name"] = ""
        df["raw_name"] = "CHANGEseq-BE"

    df["source"] = "customer"
    return df

def load_uncoverseq_csv(path: Optional[str]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=[
            "chr","start","stop","strand","Tier","Tier_num","Reproducibility","uncov_repro",
            "source","raw_name","__input_order__","UNCOV_nonon_index","site_type","UNCOV_ld"
        ])
    df = _read_table_auto(path)
    df = _canonicalize_columns(df, {
        "chr":"chr", "start":"start", "stop":"stop", "strand":"strand",
        "tier":"Tier", "reproducibility":"Reproducibility", "sitetype":"Site Type",
        "levenshteindistance":"Levenshtein Distance", "annotation":"Annotation"
    })
    for col in ["chr","start","stop","strand"]:
        if col not in df.columns:
            raise ValueError(f"UNCOVERseq: missing required column '{col}' after normalization")
    df["chr"] = df["chr"].astype(str).str.strip()
    df["start"] = df["start"].apply(_safe_int)
    df["stop"] = df["stop"].apply(_safe_int)
    df["strand"]= df["strand"].astype(str).str.strip()
    if "Tier" in df.columns:
        df["Tier_num"] = df["Tier"].apply(_safe_int)
    else:
        df["Tier_num"] = 9
    if "Reproducibility" not in df.columns:
        df["Reproducibility"] = pd.NA
    df["uncov_repro"] = df["Reproducibility"].apply(_parse_repro)
    _ensure_str_col(df, "Site Type", "")
    is_on = df["Site Type"].str.contains(r"\bon\b", case=False, na=False)
    # UNCOV LD — compute first so indexes can be assigned on the LD<7 subset
    if "Levenshtein Distance" in df.columns:
        df["UNCOV_ld"] = df["Levenshtein Distance"].apply(_safe_float)
    elif "levenshtein_distance" in df.columns:
        df["UNCOV_ld"] = df["levenshtein_distance"].apply(_safe_float)
    else:
        df["UNCOV_ld"] = pd.NA
    # Assign __input_order__ and UNCOV_nonon_index only on rows that pass LD < 7
    pass_ld = df["UNCOV_ld"].notna() & (df["UNCOV_ld"] < 7)
    df.attrs["uncov_ld_lt7_count"] = int(pass_ld.sum())
    df.attrs["uncov_ld_ge7_or_missing_count"] = int((~pass_ld).sum())
    df["__input_order__"] = None
    input_idx = 1
    for i in df.index:
        if pass_ld.at[i]:
            df.at[i, "__input_order__"] = input_idx
            input_idx += 1
    df["UNCOV_nonon_index"] = None
    nonon_idx = 1
    for i in df.index:
        if pass_ld.at[i] and not is_on.at[i]:
            df.at[i, "UNCOV_nonon_index"] = nonon_idx
            nonon_idx += 1
    df["site_type"] = df["Site Type"]
    df["source"] = "UNCOVERseq"
    df["raw_name"] = "UNCOVERseq"
    return df

def load_insilico_csv(path: Optional[str], source_label: str) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    df = _read_table_auto(path)
    df = _canonicalize_columns(df, {
        "pos":"pos", "levenshteindistance":"levenshtein_distance",
        "tier":"Tier", "mar":"MAR", "gam_region":"gam_region"
    })
    for col in ["pos","levenshtein_distance","Tier","MAR"]:
        if col not in df.columns:
            raise ValueError(f"In-silico ({source_label}): missing required column '{col}' after normalization")
    rows = [parse_pos_field(p) for p in df["pos"].astype(str)]
    df["chr"] = [r[0] for r in rows]
    df["start"] = [r[1] for r in rows]
    df["stop"] = [r[2] for r in rows]
    df["strand"] = [r[3] for r in rows]
    df["MAR"] = df["MAR"].map(_norm_mar)
    df["Tier_num"] = df["Tier"].apply(_safe_int)
    df["ld"] = df["levenshtein_distance"].apply(_safe_float)
    _ensure_str_col(df, "gam_region", "")
    df["source"] = source_label
    df["raw_name"] = source_label
    df["mar_rank"] = df["MAR"].map(lambda x: MAR_ORDER.get(x,2))
    # Sort to define per-list order
    df = df.sort_values(["ld","Tier_num","mar_rank"], ascending=[True,True,True], kind="mergesort")
    df["INS_index"] = range(1, len(df)+1)
    return df

# ------------------------------
# Record, Component & merging
# ------------------------------
@dataclass
class Component:
    """One source row that merged into a final record (by coord)."""
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
    # Original (pre-merge) coordinates (used in assay_mapping)
    orig_chr: Optional[str] = None
    orig_start: Optional[int] = None
    orig_stop: Optional[int] = None
    orig_strand: Optional[str] = None
    # --- Pass-through / synthesized fields for assay_mapping.csv ---
    name: Optional[str] = None
    location: Optional[str] = None
    lev_dist: Optional[str] = None
    gene: Optional[str] = None
    annot: Optional[str] = None
    pool: Optional[str] = None
    mar_str: Optional[str] = None
    tier_str: Optional[str] = None
    avg_pct_on_target_umi: Optional[str] = None
    cseq_avg_reads: Optional[str] = None
    cseq_edit_distance: Optional[str] = None

    def display_token(self) -> str:
        """Human-friendly component token per source for mapping file."""
        s = self.source
        if s == "customer":
            x = self.cust_nonon_index if self.cust_nonon_index is not None else 1
            return f"CHANGEseq-BE{x}"
        elif s == "UNCOVERseq":
            y = self.uncov_tier if self.uncov_tier is not None else 9
            x = self.uncov_nonon_index if self.uncov_nonon_index is not None else 1
            return f"Tier{y}_{x}"
        elif s == "insilico_pop":
            return f"insilico_pop{self.ins_index or 1}"
        elif s == "insilico_ref":
            return f"insilico_ref{self.ins_index or 1}"
        return s

@dataclass
class Record:
    chr: str
    start: int
    stop: int
    strand: str
    name_tokens: List[str]
    # Provenance TOKENS for 'orig_Name' (always ';'-joined tokens, even if no merges)
    # Example: ["Tier3_3","CHANGEseq-BE1","insilico_pop65","insilico_ref44"]
    name_list: List[str]
    sources_present: Set[str]
    components: List[Component]
    is_on: bool = False
    # UNCOV
    uncov_present: bool = False
    uncov_is_on: bool = False
    uncov_tier: Optional[int] = None
    uncov_nonon_index: Optional[int] = None
    uncov_min_order: Optional[int] = None
    uncov_ld: Optional[float] = None
    uncov_repro: Optional[int] = None # 1..3 only
    # Customer
    customer_present: bool = False
    customer_is_on: bool = False
    customer_nonon_index: Optional[int] = None
    customer_name: Optional[str] = None
    cust_min_order: Optional[int] = None
    cust_avg_reads_max: Optional[float] = None
    cust_edit_distance_min: Optional[float] = None
    # In-silico
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

def merge_sources(df_customer: pd.DataFrame,
                  df_uncov: pd.DataFrame,
                  df_pop: pd.DataFrame,
                  df_ref: pd.DataFrame,
                  dedupe_mode: str) -> Tuple[List[Record], List[Record]]:
    order = [df_customer, df_uncov, df_pop, df_ref]
    labels = ["customer","UNCOVERseq","insilico_pop","insilico_ref"]
    all_map: Dict[Tuple[str,int,int,str], Record] = {}

    def add_row(row: pd.Series, label: str):
        key = (str(row["chr"]), int(row["start"]), int(row["stop"]), str(row["strand"]) or ".")
        site_type_val = str(row.get("site_type",""))
        is_on = bool(re.search(r"\bon\b", site_type_val.lower()))
        tier = _safe_int(row.get("Tier", row.get("Tier_num", pd.NA)), None)

        if dedupe_mode == "none":
            rec = Record(key[0], key[1], key[2], key[3],
                         name_tokens=[],
                         name_list=[],
                         sources_present=set(),
                         components=[])
            all_map[(id(rec), len(all_map))] = rec
        else:
            rec = all_map.get(key)
            if rec is None:
                rec = Record(key[0], key[1], key[2], key[3],
                             name_tokens=[],
                             name_list=[],
                             sources_present=set(),
                             components=[])
                all_map[key] = rec

        rec.sources_present.add(label)
        rec.is_on = rec.is_on or is_on

        # Build component
        comp = Component(source=label, site_type=site_type_val, raw_name=str(row.get("raw_name","")))
        # Original coordinates
        comp.orig_chr = key[0]
        comp.orig_start = key[1]
        comp.orig_stop = key[2]
        comp.orig_strand = key[3]

        # ---- Pass-through / synthesized fields for assay_mapping.csv ----
        comp.name     = (None if pd.isna(row.get("Name", None)) else str(row.get("Name", "") or "") or None)
        # Prefer existing Location; otherwise synthesize/derive per source
        loc_in = row.get("Location", None)
        comp.location = (None if pd.isna(loc_in) else str(loc_in or "") or None)
        # Accept either short header "Lev. Dist." or canonical "Levenshtein Distance"
        ld_disp = row.get("Lev. Dist.", row.get("Levenshtein Distance", None))
        comp.lev_dist = (None if pd.isna(ld_disp) else str(ld_disp or "") or None)
        comp.gene     = (None if pd.isna(row.get("Gene", None)) else str(row.get("Gene", "") or "") or None)
        annot_in = row.get("Annot.", row.get("Annotation", None))
        if (annot_in is None or pd.isna(annot_in)) and label.startswith("insilico"):
            annot_in = row.get("gam_region", None)
        comp.annot    = (None if pd.isna(annot_in) else str(annot_in or "") or None)
        comp.pool     = (None if pd.isna(row.get("Pool", None)) else str(row.get("Pool", "") or "") or None)
        comp.mar_str  = (None if pd.isna(row.get("MAR", None)) else str(row.get("MAR", "") or "") or None)
        tier_disp = row.get("Tier", row.get("Tier_num", None))
        comp.tier_str = (None if pd.isna(tier_disp) else str(tier_disp or "") or None)
        comp.avg_pct_on_target_umi = (None if pd.isna(row.get("Avg. % On-Target UMI", None))
                                      else str(row.get("Avg. % On-Target UMI", "") or "") or None)

        # If Location still missing, derive it now
        if comp.location is None or str(comp.location).strip() == "":
            if label.startswith("insilico") and not pd.isna(row.get("pos", None)):
                comp.location = str(row.get("pos"))  # original pos string
            else:
                comp.location = _build_location(key[0], key[1], key[2], key[3])

        if label == "UNCOVERseq":
            rec.uncov_present = True
            rec.uncov_is_on = rec.uncov_is_on or is_on
            if tier is not None:
                rec.uncov_tier = min(rec.uncov_tier, tier) if rec.uncov_tier is not None else tier
            uncov_nonon_index = _safe_int(row.get("UNCOV_nonon_index", pd.NA), None)
            if uncov_nonon_index:
                rec.uncov_nonon_index = min(rec.uncov_nonon_index, uncov_nonon_index) if rec.uncov_nonon_index else uncov_nonon_index
            in_ord = _safe_int(row.get("__input__order__", row.get("__input_order__", pd.NA)), None)
            if in_ord is None:
                in_ord = _safe_int(row.get("__input_order__", pd.NA), None)
            if in_ord is not None:
                rec.uncov_min_order = min(rec.uncov_min_order, in_ord) if rec.uncov_min_order is not None else in_ord
            ld_uncov = _safe_float(row.get("UNCOV_ld", pd.NA), None)
            if ld_uncov is not None:
                rec.uncov_ld = min(rec.uncov_ld, ld_uncov) if rec.uncov_ld is not None else ld_uncov
            repro = _parse_repro(row.get("uncov_repro", row.get("Reproducibility", None)))
            if repro is not None:
                rec.uncov_repro = max(rec.uncov_repro, repro) if rec.uncov_repro is not None else repro
            comp.uncov_tier = _safe_int(row.get("Tier_num", row.get("Tier", None)), None)
            comp.uncov_nonon_index = uncov_nonon_index
            comp.uncov_input_order = in_ord
            comp.uncov_ld = ld_uncov
            comp.uncov_repro = repro

        elif label == "customer":
            rec.customer_present = True
            rec.customer_is_on = rec.customer_is_on or is_on
            cust_nonon_idx = _safe_int(row.get("CUST_nonon_index", pd.NA), None)
            if cust_nonon_idx:
                rec.customer_nonon_index = min(rec.customer_nonon_index, cust_nonon_idx) if rec.customer_nonon_index else cust_nonon_idx
            cin = _safe_int(row.get("__cust_input_order__", pd.NA), None)
            if cin is not None:
                rec.cust_min_order = min(rec.cust_min_order, cin) if rec.cust_min_order is not None else cin
            cavg = _safe_float(row.get("CSEQ_avg_reads", pd.NA), None)
            if cavg is not None:
                rec.cust_avg_reads_max = max(rec.cust_avg_reads_max, cavg) if rec.cust_avg_reads_max is not None else cavg
            ced = _safe_float(row.get("CSEQ_edit_distance", pd.NA), None)
            if ced is not None:
                rec.cust_edit_distance_min = min(rec.cust_edit_distance_min, ced) if rec.cust_edit_distance_min is not None else ced
            if cavg is not None:
                comp.cseq_avg_reads = str(cavg)
            if ced is not None:
                comp.cseq_edit_distance = str(ced)
            provided = row.get("cust_name", "")
            if isinstance(provided, str):
                nm = _sanitize_name(provided)
                if nm:
                    rec.customer_name = rec.customer_name or nm
            comp.cust_name = _sanitize_name(row.get("cust_name", ""))
            comp.cust_nonon_index = cust_nonon_idx
            comp.cust_input_order = cin

        elif label in ("insilico_pop","insilico_ref"):
            is_pop = (label == "insilico_pop")
            if is_pop:
                rec.ins_pop_present = True
                idx_val = _safe_int(row.get("INS_index", pd.NA), None)
                rec.ins_pop_index = min(rec.ins_pop_index, idx_val) if (rec.ins_pop_index and idx_val) else (rec.ins_pop_index or idx_val)
            else:
                rec.ins_ref_present = True
                idx_val = _safe_int(row.get("INS_index", pd.NA), None)
                rec.ins_ref_index = min(rec.ins_ref_index, idx_val) if (rec.ins_ref_index and idx_val) else (rec.ins_ref_index or idx_val)
            ld_ins = _safe_float(row.get("ld", pd.NA), None)
            if ld_ins is not None:
                rec.ins_ld_min = min(rec.ins_ld_min, ld_ins) if rec.ins_ld_min is not None else ld_ins
            if tier is not None:
                rec.ins_tier_min = min(rec.ins_tier_min, tier) if rec.ins_tier_min is not None else tier
            mar_rank = _safe_int(row.get("mar_rank", pd.NA), None)
            if mar_rank is not None:
                rec.ins_mar_rank_min = min(rec.ins_mar_rank_min, mar_rank) if rec.ins_mar_rank_min is not None else mar_rank
            gr = str(row.get("gam_region","") or "")
            if ld_ins == 3 and gr.strip().lower() == "exon":
                rec.ins_ld3_exon = True
            comp.ins_index = _safe_int(row.get("INS_index", pd.NA), None)
            comp.ins_ld = ld_ins
            comp.ins_tier = _safe_int(row.get("Tier_num", row.get("Tier", None)), None)
            comp.ins_mar_rank = mar_rank
            comp.gam_region = gr

        # Append component to record
        rec.components.append(comp)

    for df, label in zip(order, labels):
        if df is not None and not df.empty:
            for _, row in df.iterrows():
                add_row(row, label)

    all_records = list(all_map.values())

    # ---- Naming with source indices ----
    on_records = [r for r in all_records if r.is_on]
    multiple_on = len(on_records) > 1
    on_count_uncov = 0
    on_count_cust = 0
    for r in all_records:
        tokens: List[str] = []
        if r.is_on:
            if not multiple_on:
                tokens = ["OnT"]
            else:
                if r.uncov_is_on:
                    on_count_uncov += 1
                    tokens = [f"OnT_UNCOVERseq_{on_count_uncov}"]
                elif r.customer_is_on:
                    on_count_cust += 1
                    tokens = [f"OnT_CHANGEseq-BE_{on_count_cust}"]
                else:
                    tokens = ["OnT"]
        else:
            if r.uncov_present:
                tier_val = r.uncov_tier if r.uncov_tier is not None else 9
                x = r.uncov_nonon_index if r.uncov_nonon_index is not None else 1
                base = f"Tier{tier_val}_{x}"
                tokens = [base]
                if r.customer_present and not r.customer_is_on:
                    cx = r.customer_nonon_index if r.customer_nonon_index is not None else 1
                    tokens.append(f"CHANGEseq-BE{cx}")
                if r.ins_pop_present:
                    tokens.append(f"insilico_pop{r.ins_pop_index or 1}")
                if r.ins_ref_present:
                    tokens.append(f"insilico_ref{r.ins_ref_index or 1}")
            elif r.customer_present:
                if r.customer_is_on:
                    tokens = ["OnT"]
                else:
                    x = r.customer_nonon_index if r.customer_nonon_index is not None else 1
                    base = f"CHANGEseq-BE{x}"
                    tokens = [base]
                    if r.ins_pop_present:
                        tokens.append(f"insilico_pop{r.ins_pop_index or 1}")
                    if r.ins_ref_present:
                        tokens.append(f"insilico_ref{r.ins_ref_index or 1}")
            else:
                if r.ins_pop_present and r.ins_ref_present:
                    tokens = [f"insilico_pop{r.ins_pop_index or 1}", f"insilico_ref{r.ins_ref_index or 1}"]
                elif r.ins_pop_present:
                    tokens = [f"insilico_pop{r.ins_pop_index or 1}"]
                elif r.ins_ref_present:
                    tokens = [f"insilico_ref{r.ins_ref_index or 1}"]
                else:
                    tokens = ["unknown"]
        r.name_tokens = tokens
        # Initialize provenance token list for orig_Name with current tokens
        r.name_list = list(r.name_tokens)

    return all_records, on_records

# ------------------------------
# Proximity/coordinate merging (bed_qc-like)
# ------------------------------
def _recompute_aggregates_from_components(r: "Record") -> None:
    comps = r.components or []
    r.sources_present = set(c.source for c in comps)
    r.is_on = any(bool(re.search(r"\bon\b", (c.site_type or "").lower())) for c in comps)
    # UNCOVERseq
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
    # Customer
    r.customer_present = any(c.source == "customer" for c in comps)
    r.customer_is_on = any((c.source == "customer") and bool(re.search(r"\bon\b", (c.site_type or "").lower())) for c in comps)
    cust_nonon = [c.cust_nonon_index for c in comps if c.cust_nonon_index]
    r.customer_nonon_index = min(cust_nonon) if cust_nonon else None
    cust_orders = [c.cust_input_order for c in comps if c.cust_input_order is not None]
    r.cust_min_order = min(cust_orders) if cust_orders else None
    cust_names = [c.cust_name for c in comps if c.cust_name]
    r.customer_name = cust_names[0] if cust_names else None
    cust_avg_reads = [_safe_float(c.cseq_avg_reads, None) for c in comps if c.source == "customer"]
    cust_avg_reads = [v for v in cust_avg_reads if v is not None]
    r.cust_avg_reads_max = max(cust_avg_reads) if cust_avg_reads else None
    cust_edit = [_safe_float(c.cseq_edit_distance, None) for c in comps if c.source == "customer"]
    cust_edit = [v for v in cust_edit if v is not None]
    r.cust_edit_distance_min = min(cust_edit) if cust_edit else None
    # In-silico
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
    r.ins_ld3_exon = any((c.ins_ld == 3) and (str(getattr(c, "gam_region","")).strip().lower() == "exon") for c in comps)

def _merge_two_records(prev: "Record", curr: "Record", name_delim: str) -> "Record":
    chr_ = prev.chr
    start = min(prev.start, curr.start)
    has_customer = prev.customer_present or curr.customer_present
    has_insilico = (prev.ins_pop_present or prev.ins_ref_present or curr.ins_pop_present or curr.ins_ref_present)
    # For same-start customer+insilico pairs, prefer tighter endpoint.
    if prev.start == curr.start and has_customer and has_insilico:
        stop = min(prev.stop, curr.stop)
    else:
        # Emulate target_qc merge endpoint behavior: merged stop follows curr.stop.
        stop = curr.stop
    # Merge regardless of strand; use '.' if they differ
    strand = prev.strand if prev.strand == curr.strand else "."
    # Join final names; cap at 50 chars
    left  = _final_name_from_record(prev)
    right = _final_name_from_record(curr)
    merged_name = f"{left}{name_delim}{right}" if left and right else (left or right or "unknown")
    if len(merged_name) > MAX_MERGED_NAME_LEN:
        merged_name = f"{chr_}:{start}-{stop}"
    # Merge provenance token lists (preserve order)
    merged_name_list = (prev.name_list or []) + (curr.name_list or [])
    merged = Record(chr_, start, stop, strand,
                    name_tokens=[merged_name],
                    name_list=merged_name_list,
                    sources_present=set(),
                    components=(prev.components + curr.components))
    _recompute_aggregates_from_components(merged)
    return merged

def merge_adjacent_records(all_records: List["Record"],
                           target_distance: int,
                           max_insert: int = 200,
                           target_flank: int = 40,
                           max_merge_size: Optional[int] = None,
                           ##name_delim to match rhAmpSeq target_qc.py
                           name_delim: str = "__") -> List["Record"]:
    """
    Proximity-based merge (bed_qc-like):
      - same chromosome
      - merge if (gap = curr.start - prev.stop) < target_distance
        AND (merged_len = curr.stop - prev.start) <= max_merge_size
      - default max_merge_size: max_insert - target_flank - 15
    """
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
            else:
                out.append(prev)
                prev = curr
        else:
            out.append(prev)
            prev = curr
    if prev is not None:
        out.append(prev)
    return out

def merge_duplicate_coords(all_records: List["Record"], name_delim: str = ";") -> List["Record"]:
    """
    Merge records that share identical (chr,start,stop,strand) by concatenating components
    and joining names with delimiter (cap >50 as in bed_qc).
    Useful when dedupe_mode == 'none'.
    """
    if not all_records:
        return all_records
    from collections import defaultdict
    groups: Dict[Tuple[str,int,int,str], List[Record]] = defaultdict(list)
    for r in all_records:
        groups[(r.chr, r.start, r.stop, r.strand or ".")].append(r)
    out: List[Record] = []
    for (chr_, start, stop, strand), items in groups.items():
        if len(items) == 1:
            out.append(items[0])
            continue
        # Start with the first, then fold in others
        base = items[0]
        for nxt in items[1:]:
            base = _merge_two_records(base, nxt, name_delim)
            # Restore original span because coords are identical
            base.chr, base.start, base.stop, base.strand = chr_, start, stop, strand
        out.append(base)
    return out

# ------------------------------
# Prioritization (per doc)
# ------------------------------
def prioritize(all_recs: List[Record], max_sites: int):
    """
        Order:
            1) Top 60 CHANGEseq-BE sites by avg rep reads (desc), then Edit_distance (asc)
            2) Top 20 in-silico sites
            3) Remaining slots from UNCOVERseq (LD < 7, original order)
    """
    selected: List[Record] = []
    seen = set()

    def add_in_order(records: List[Record], key=None, limit=None):
        nonlocal selected, seen
        out = 0
        for r in sorted(records, key=key) if key else records:
            if len(selected) >= max_sites:
                break
            k = (r.chr, r.start, r.stop, r.strand)
            if k in seen:
                continue
            selected.append(r)
            seen.add(k)
            out += 1
            if limit is not None and out >= limit:
                break
        return out

    # 1) CHANGEseq-BE block (60)
    cust_block = [r for r in all_recs if r.customer_present]
    add_in_order(
        cust_block,
        key=lambda r: (
            -(r.cust_avg_reads_max if r.cust_avg_reads_max is not None else -1e9),
            (r.cust_edit_distance_min if r.cust_edit_distance_min is not None else 1e9),
            (r.cust_min_order if r.cust_min_order is not None else 10**12),
        ),
        limit=60,
    )
    if len(selected) >= max_sites:
        return selected, []

    # 2) In-silico block (20)
    ins_pool = [r for r in all_recs
                if (r.ins_pop_present or r.ins_ref_present)
                and (r.chr, r.start, r.stop, r.strand) not in seen]
    add_in_order(
        ins_pool,
        key=lambda r: (
            r.ins_ld_min if r.ins_ld_min is not None else 1e9,
            r.ins_tier_min if r.ins_tier_min is not None else 9,
            r.ins_mar_rank_min if r.ins_mar_rank_min is not None else 2,
            min([x for x in [r.ins_pop_index, r.ins_ref_index] if x is not None], default=10**9),
        ),
        limit=20,
    )
    if len(selected) >= max_sites:
        return selected, []

    # 3) UNCOVERseq remainder (apply LD < 7)
    uncov_pool_all = [r for r in all_recs
                      if r.uncov_present and r.uncov_ld is not None and r.uncov_ld < 7
                      and (r.chr, r.start, r.stop, r.strand) not in seen]
    add_in_order(
        uncov_pool_all,
        key=lambda r: (r.uncov_min_order if r.uncov_min_order is not None else 10**12),
    )

    # Build excluded (stable source-aware order)
    all_by_source_order = []
    all_by_source_order.extend(sorted([r for r in all_recs if r.customer_present],
                                      key=lambda r: (
                                          -(r.cust_avg_reads_max if r.cust_avg_reads_max is not None else -1e9),
                                          (r.cust_edit_distance_min if r.cust_edit_distance_min is not None else 1e9),
                                          (r.cust_min_order if r.cust_min_order is not None else 10**12),
                                      )))
    ins_tail = sorted([r for r in all_recs if (r.ins_pop_present or r.ins_ref_present)],
                      key=lambda r: (
                          r.ins_ld_min if r.ins_ld_min is not None else 1e9,
                          r.ins_tier_min if r.ins_tier_min is not None else 9,
                          r.ins_mar_rank_min if r.ins_mar_rank_min is not None else 2,
                          min([x for x in [r.ins_pop_index, r.ins_ref_index] if x is not None], default=10**9)
                      ))
    all_by_source_order.extend(ins_tail)
    all_by_source_order.extend(sorted([r for r in all_recs if r.uncov_present],
                                      key=lambda r: (r.uncov_min_order if r.uncov_min_order is not None else 10**12)))

    # Dedup for excluded
    seen2 = set((r.chr, r.start, r.stop, r.strand) for r in selected)
    excluded = []
    for r in all_by_source_order:
        key = (r.chr, r.start, r.stop, r.strand)
        if key in seen2:
            continue
        seen2.add(key)
        excluded.append(r)
    return selected, excluded

# ------------------------------
# I/O
# ------------------------------
def write_bed(path: str, recs: List[Record]) -> int:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    count = 0
    with open(path, "w") as f:
        for r in recs:
            f.write("\t".join(r.bed_row()) + "\n")
            count += 1
    return count
    
def build_dedup_nonmerged_prioritized_rows(prioritized: list) -> list:
    """
    Build a prioritized list of original (non-merged) component coordinates,
    deduplicated by (chr,start,stop,strand). Each row: [chr, start, stop, name, strand]
    'name' is the ORIGINAL component label, not the merged BED name:
      - 'OnT' if that component is ON-target (from its own Site Type)
      - UNCOVERseq: 'Tier{tier}_{nonon_index}'
            - customer:  'CHANGEseq-BE{index}'
      - insilico_pop/ref: 'Insilico_pop{index}' / 'Insilico_ref{index}'
    """
    rows = []
    seen = set()

    def comp_is_on(c):
        st = (c.site_type or '').lower()
        return bool(re.search(r"\bon\b", st))

    for r in prioritized:
        for c in (r.components or []):
            chr_   = c.orig_chr   or r.chr
            start  = c.orig_start if c.orig_start is not None else r.start
            stop   = c.orig_stop  if c.orig_stop  is not None else r.stop
            strand = (c.orig_strand or r.strand or ".")

            key = (chr_, int(start), int(stop), strand)
            if key in seen:
                continue
            seen.add(key)

            # --- original per-component 'Name' (no merged tokens) ---
            if comp_is_on(c):
                orig_name = "OnT"
            elif c.source == "UNCOVERseq":
                t = c.uncov_tier if c.uncov_tier is not None else 9
                x = c.uncov_nonon_index if c.uncov_nonon_index is not None else 1
                orig_name = f"Tier{t}_{x}"
            elif c.source == "customer":
                x = c.cust_nonon_index if c.cust_nonon_index is not None else 1
                orig_name = f"CHANGEseq-BE{x}"
            elif c.source == "insilico_pop":
                orig_name = f"insilico_pop{c.ins_index or 1}"
            elif c.source == "insilico_ref":
                orig_name = f"insilico_ref{c.ins_index or 1}"
            else:
                # Fallback: display token or raw name
                try:
                    orig_name = c.display_token()
                except Exception:
                    orig_name = c.raw_name or "orig"

            rows.append([chr_, str(int(start)), str(int(stop)), orig_name, "0", strand])

    return rows

def write_bed6(path: str, rows: list) -> int:
    """
    Write a 6-column BED-like file:
    chr, start, stop, name, 0, strand

    Coordinates are assumed 0-based;
    'strand' is '.' or '+' or '-'.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    n = 0
    with open(path, "w") as f:
        for r in rows:
            # Ensure exactly 6 columns are written
            line = "\t".join([str(r[0]), str(r[1]), str(r[2]), str(r[3]), str(r[4]), str(r[5])])
            f.write(line + "\n")
            n += 1
    return n

def _semi_join(values: List[Optional[str]]) -> str:
    """
    Join with ';' after normalizing blanks to 'NONE'.
    If all normalized values are identical, collapse to that single value.
    Preserve input order; do NOT de-duplicate unless all values match.
    """
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
    # Collapse insilico_* to 'insilico'
    if source.startswith("insilico"):
        return "insilico"
    if source == "customer":
        return "CHANGEseq-BE"
    return source

def _component_sort_key(c: Component) -> int:
    """
    Ordering priority for aggregation:
    CHANGEseq-BE (0) < insilico (1) < UNCOVERseq (2)
    Stable sort keeps original per-source order.
    """
    s = c.source or ""
    if s == "customer":
        return 0
    if s.startswith("insilico"):
        return 1
    if s == "UNCOVERseq":
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
    """
    Emit assay_mapping.csv (one row per FINAL target).
    Columns:
      chr, start, stop, strand, Name, Location, Lev. Dist., Gene, Annot., Pool,
      MAR, Tier, Tier_Final, Site Type, Origin, Avg. % On-Target UMI,
      orig_chr, orig_start, orig_stop, orig_strand, orig_Name

    Values are ';'-joined across merged originals in order:
            CHANGEseq-BE -> insilico -> UNCOVERseq.
        Output 'Annot.' is sourced from UNCOVERseq 'Annotation' and insilico 'gam_region'; CHANGEseq-BE contributes 'NONE'.
    Blanks become 'NONE'. If all values in a column are identical, collapse to a single value.
    """
    import csv
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    headers = [
        "chr","start","stop","strand",
        "Name","Location","Lev. Dist.","Gene","Annot.","Pool",
        "MAR","Tier","Tier_Final","Site Type","Origin","Avg. % On-Target UMI",
        "CHANGEseq-BE avg reads","CHANGEseq-BE Edit_distance",
        "orig_chr","orig_start","orig_stop","orig_strand","orig_Name",
    ]
    n = 0
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in recs:
            comps = sorted(r.components, key=_component_sort_key)

            # Final bed name (4th BED column, possibly capped to <=50 chars or replaced by span)
            final_name = "_".join([t for t in r.name_tokens if t])
            # Provenance tokens for orig_Name (semicolon-joined)
            name_uncapped = ";".join([s for s in (r.name_list or []) if s is not None and str(s).strip() != ""])

            # Collect per-field lists in the specified order
            locs       = [c.location for c in comps]
            levd       = [c.lev_dist for c in comps]
            genes      = [c.gene     for c in comps]
            annots     = [c.annot    for c in comps]
            pools      = [c.pool     for c in comps]
            mars       = [c.mar_str  for c in comps]
            tiers      = [c.tier_str for c in comps]
            stypes     = [c.site_type for c in comps]
            origins    = [_origin_display(c.source or "") for c in comps]
            avg_pct    = [c.avg_pct_on_target_umi for c in comps]
            cseq_avg   = [c.cseq_avg_reads for c in comps if c.source == "customer"]
            cseq_ed    = [c.cseq_edit_distance for c in comps if c.source == "customer"]
            # Original coords
            o_chr      = [c.orig_chr   for c in comps]
            o_start    = [str(c.orig_start) if c.orig_start is not None else None for c in comps]
            o_stop     = [str(c.orig_stop)  if c.orig_stop  is not None else None for c in comps]
            o_strand   = [c.orig_strand for c in comps]

            row = [
                r.chr, r.start, r.stop, (r.strand or "."),
                final_name,               # Name = actual BED name (capped/coord as needed)
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
                _semi_join(cseq_avg if cseq_avg else [None]),
                _semi_join(cseq_ed if cseq_ed else [None]),
                _semi_join(o_chr),
                _semi_join(o_start),
                _semi_join(o_stop),
                _semi_join(o_strand),
                name_uncapped if name_uncapped else final_name,  # orig_Name = tokens joined with ';'
            ]
            w.writerow(row)
            n += 1
    return n

def write_merge_mapping(path: str,
                        all_records: List[Record],
                        prioritized: List[Record]) -> int:
    """
    Writes a TSV mapping of merged inputs -> final output name.
    Columns:
      final_name, chr, start, stop, strand, prioritized, rank,
      sources_present, component_sources, component_names
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    # Build rank map for prioritized
    rank = {}
    for i, r in enumerate(prioritized, start=1):
        rank[(r.chr, r.start, r.stop, r.strand)] = i
    cols = [
        "final_name","chr","start","stop","strand","prioritized","rank",
        "sources_present","component_sources","component_names"
    ]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in all_records:
            final_name = "_".join([t for t in r.name_tokens if t])
            key = (r.chr, r.start, r.stop, r.strand)
            is_sel = key in rank
            rnk = rank.get(key, "")
            comp_sources = [c.source for c in r.components]
            comp_names = [c.display_token() for c in r.components]
            row = {
                "final_name": final_name,
                "chr": r.chr,
                "start": r.start,
                "stop": r.stop,
                "strand": r.strand or ".",
                "prioritized": "Y" if is_sel else "N",
                "rank": rnk,
                "sources_present": ";".join(sorted(list(r.sources_present))),
                "component_sources": ";".join(comp_sources),
                "component_names": ";".join(comp_names),
            }
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
    return 0

# ------------------------------
# CLI
# ------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description="Nomination site prioritization → rhAmpSeq")
    p.add_argument("--changeseq_be", type=str, help="CHANGEseq-BE CSV/TSV with headers including Chromosome/Start/End/Strand and CHANGEseq-BE_rep*_reads")
    p.add_argument("--uncoverseq", type=str, help="UNCOVERseq CSV/TSV with headers: chr,start,stop,strand (+ optional Tier/Reproducibility/Site Type)")
    p.add_argument("--insilico_pop", type=str, help="In-silico POP CSV/TSV with headers: pos,levenshtein_distance,Tier,MAR")
    p.add_argument("--insilico_ref", type=str, help="In-silico REF CSV/TSV with headers: pos,levenshtein_distance,Tier,MAR")
    p.add_argument("--max_sites", type=int, default=None, help="Override global cap (default 200)")
    p.add_argument("--dedupe", choices=["by_coord","none"], default="by_coord", help="Default: by_coord")
    p.add_argument("--out_prefix", type=str, default=None, help="Path + filename prefix for outputs")
    p.add_argument("--outdir", type=str, default=None, help="Output directory (used only if --out_prefix not provided)")
    p.add_argument("--verbose", action="store_true", help="Echo banner with resolved paths, settings, and input counts")
    # --- Target merging (bed_qc-like) ---
    p.add_argument("--merge_targets",
                   choices=["proximity","none","dedupe_coords"],
                   default="proximity",
                   help="Target merge strategy (default: proximity — match rhAmpSeq design bed qc). "
                        "Use 'none' to disable or 'dedupe_coords' to merge only exact-duplicate coordinates.")
    p.add_argument("--target_distance", type=int, default=500,
                   help="Gap between adjacent targets (bp) for proximity merge: merge if (next.start - prev.stop) < target_distance [default: 500 to match rhAmpSeq design bed qc]")
    p.add_argument("--max_insert", type=int, default=200,
                   help="Max insert used to derive max_merge_size when not provided (default: 200)")
    p.add_argument("--target_flank", type=int, default=40,
                   help="Flank length subtracted when deriving max_merge_size (default: 20)")
    p.add_argument("--max_merge_size", type=int, default=None,
                   help="Explicit cap on merged span length; overrides derived value if set")
    p.add_argument("--merge_name_delim", type=str, default="__",
                   help="Delimiter used when joining names during merge. default matches rhAmpSeq (default: '__')")
    args = p.parse_args(argv)

    MAX_SITES_RUNTIME = MAX_SITES if args.max_sites is None else int(args.max_sites)
    if MAX_SITES_RUNTIME <= 0:
        raise SystemExit("--max_sites must be positive")

    if args.out_prefix:
        out_dir = os.path.dirname(args.out_prefix) or "."
        os.makedirs(out_dir, exist_ok=True)
        base = os.path.basename(args.out_prefix)
        prioritized_path = os.path.join(out_dir, f"{base}_prioritized.bed")
        excluded_path = os.path.join(out_dir, f"{base}_excluded.bed")
        full_path = os.path.join(out_dir, f"{base}_all_candidates_full.bed")
        mapping_path = os.path.join(out_dir, f"{base}_merge_mapping.tsv")
        dedup_nonmerged_prioritized_path = os.path.join(out_dir, f"{base}_dedup_nonmerged_prioritized_sites.bed")
    else:
        out_dir = args.outdir or "."
        os.makedirs(out_dir, exist_ok=True)
        prioritized_path = os.path.join(out_dir, "prioritized_sites.bed")
        excluded_path = os.path.join(out_dir, "excluded_sites.bed")
        full_path = os.path.join(out_dir, "all_candidates_full.bed")
        mapping_path = os.path.join(out_dir, "merge_mapping.tsv")
        dedup_nonmerged_prioritized_path = os.path.join(out_dir, "dedup_nonmerged_prioritized_sites.bed")

    # Assay mapping CSV paths
    if args.out_prefix:
        base = os.path.basename(args.out_prefix)
        am_prior = os.path.join(out_dir, f"{base}_prioritized_assay_mapping.csv")
        am_excl  = os.path.join(out_dir, f"{base}_excluded_assay_mapping.csv")
        am_all   = os.path.join(out_dir, f"{base}_all_candidates_assay_mapping.csv")
    else:
        am_prior = os.path.join(out_dir, "prioritized_assay_mapping.csv")
        am_excl  = os.path.join(out_dir, "excluded_assay_mapping.csv")
        am_all   = os.path.join(out_dir, "all_candidates_assay_mapping.csv")

    # Load sources
    df_customer = load_changeseq_be_csv(args.changeseq_be) if args.changeseq_be else pd.DataFrame()
    df_uncov = load_uncoverseq_csv(args.uncoverseq) if args.uncoverseq else pd.DataFrame()
    df_pop = load_insilico_csv(args.insilico_pop, "insilico_pop") if args.insilico_pop else pd.DataFrame()
    df_ref = load_insilico_csv(args.insilico_ref, "insilico_ref") if args.insilico_ref else pd.DataFrame()

    if args.verbose:
        sep = "-" * 60
        print(sep)
        print("VERBOSE RUN SETTINGS")
        print(sep)
        print(f"MAX_SITES (runtime): {MAX_SITES_RUNTIME}")
        print(f"Deduplication mode: {args.dedupe}")
        print(f"Merge strategy: {args.merge_targets}")
        if args.merge_targets == "proximity":
            derived_mms = (args.max_insert - args.target_flank - 15) if args.max_merge_size is None else args.max_merge_size
            print(f" target_distance: {args.target_distance}")
            print(f" max_insert: {args.max_insert}")
            print(f" target_flank: {args.target_flank}")
            print(f" max_merge_size: {derived_mms}")
            print(f" merge_name_delim: {args.merge_name_delim}")
        print()
        print("Output paths:")
        print(f" prioritized: {prioritized_path}")
        print(f" excluded: {excluded_path}")
        print(f" full list: {full_path}")
        print(f" mapping: {mapping_path}")
        print(f" assay map (prioritized): {am_prior}")
        print(f" assay map (excluded):    {am_excl}")
        print(f" assay map (all):         {am_all}")
        print(f" dedup (non-merged, prioritized): {dedup_nonmerged_prioritized_path}")
        print()
        print("Input summary:")
        print(f" CHANGEseq-BE rows: {len(df_customer)}")
        print(f" UNCOVERseq rows: {len(df_uncov)}")
        print(f" In-silico POP: {len(df_pop)}")
        print(f" In-silico REF: {len(df_ref)}")
        print(sep)
        print()

    # Merge & (optionally) proximity/coord-merge, then prioritize
    all_records, _ = merge_sources(df_customer, df_uncov, df_pop, df_ref, args.dedupe)

    if args.merge_targets == "proximity":
        all_records = merge_adjacent_records(
            all_records,
            target_distance=args.target_distance,
            max_insert=args.max_insert,
            target_flank=args.target_flank,
            max_merge_size=args.max_merge_size,
            name_delim=args.merge_name_delim
        )
    elif args.merge_targets == "dedupe_coords":
        all_records = merge_duplicate_coords(all_records, name_delim=args.merge_name_delim)
    # else: 'none' => leave as-is

    final, excluded = prioritize(all_records, MAX_SITES_RUNTIME)
    uncov_ld_lt7_excluded = sum(
        1
        for r in excluded
        if r.uncov_present and r.uncov_ld is not None and r.uncov_ld < 7
    )

    # Write outputs & report counts
    n_prioritized = write_bed(prioritized_path, final)
    n_excluded = write_bed(excluded_path, excluded)
    n_full = write_bed(full_path, all_records)
    write_merge_mapping(mapping_path, all_records, final)
    # Write the non-merged, deduplicated prioritized 6-col BED
    rows_nm = build_dedup_nonmerged_prioritized_rows(final)
    n_nm = write_bed6(dedup_nonmerged_prioritized_path, rows_nm)
    n_nm = 0 if n_nm is None else int(n_nm)

    # Assay mapping CSVs
    n_am_prior = write_assay_mapping(am_prior, final)
    n_am_excl  = write_assay_mapping(am_excl,  excluded)
    n_am_all   = write_assay_mapping(am_all,   all_records)

    print(f"Wrote {n_prioritized} lines: {prioritized_path}")
    print(f"Wrote {n_excluded} lines: {excluded_path}")
    print(f"UNCOVERseq LD<7 sites excluded (in excluded_sites list): {uncov_ld_lt7_excluded}")
    print(f"Wrote {n_full} lines: {full_path}")
    print(f"Wrote mapping: {mapping_path}")
    print(f"Wrote {n_nm:>4} lines: {dedup_nonmerged_prioritized_path}")
    print(f"Wrote {n_am_prior:>4} rows: {am_prior}")
    print(f"Wrote {n_am_excl:>4} rows: {am_excl}")
    print(f"Wrote {n_am_all:>4} rows: {am_all}")

    if args.verbose:
        print("Completed.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
