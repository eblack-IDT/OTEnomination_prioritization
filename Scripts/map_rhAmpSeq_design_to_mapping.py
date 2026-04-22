#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


POOL_BY_PREFIX = {
    "p0": "Primary pool",
    "p2": "Secondary pool",
    "sc": "Singleton collection",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Populate the Pool column in assay_mapping.csv using rhAmpSeq ADF design files "
            "(p0.<project>.adf, p2.<project>.adf, sc.<project>.adf)."
        )
    )
    parser.add_argument(
        "--adf-dir",
        type=Path,
        required=True,
        help="Directory containing p0/p2/sc ADF files",
    )
    parser.add_argument(
        "--assay-mapping",
        type=Path,
        required=True,
        help="Path to existing assay_mapping.csv",
    )
    parser.add_argument(
        "--projectname",
        type=str,
        required=True,
        help="Project name used in ADF file names (e.g., p0.<projectname>.adf)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output CSV path (default: overwrite assay_mapping.csv in place)."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting output file when --output is provided",
    )
    return parser.parse_args()


def read_table_auto(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    df.columns = [col.encode("utf-8-sig").decode("utf-8-sig").strip() for col in df.columns]
    return df


def load_tids_by_pool(adf_dir: Path, projectname: str) -> dict[str, set[str]]:
    tids_by_pool: dict[str, set[str]] = {}

    for prefix, pool_label in POOL_BY_PREFIX.items():
        adf_path = adf_dir / f"{prefix}.{projectname}.adf"
        if not adf_path.exists():
            tids_by_pool[pool_label] = set()
            continue

        df = read_table_auto(adf_path)
        if "TID" not in df.columns:
            raise KeyError(f"Missing required column 'TID' in {adf_path}")

        tids = {
            str(value).strip()
            for value in df["TID"].tolist()
            if str(value).strip() != ""
        }
        tids_by_pool[pool_label] = tids

    return tids_by_pool


def build_name_to_pools(tids_by_pool: dict[str, set[str]]) -> dict[str, list[str]]:
    all_names: set[str] = set()
    for tids in tids_by_pool.values():
        all_names.update(tids)

    ordered_pool_labels = [
        POOL_BY_PREFIX["p0"],
        POOL_BY_PREFIX["p2"],
        POOL_BY_PREFIX["sc"],
    ]

    name_to_pools: dict[str, list[str]] = {}
    for name in all_names:
        labels = [label for label in ordered_pool_labels if name in tids_by_pool.get(label, set())]
        name_to_pools[name] = labels
    return name_to_pools


def resolve_output_path(mapping_path: Path, output_arg: Path | None) -> Path:
    default_name = f"{mapping_path.stem}.updated{mapping_path.suffix}"

    if output_arg is None:
        return mapping_path.with_name(default_name)

    resolved_output = output_arg.resolve()
    if resolved_output.exists() and resolved_output.is_dir():
        return resolved_output / default_name
    if output_arg.suffix == "":
        return resolved_output / default_name
    return resolved_output


def build_singleton_pcr2_pool(mapping_df: pd.DataFrame) -> pd.Series:
    if "Singleton_PCR2_pool" in mapping_df.columns:
        singleton_values = mapping_df["Singleton_PCR2_pool"].astype(str)
    else:
        singleton_values = pd.Series("", index=mapping_df.index, dtype=str)

    return singleton_values.where(
        mapping_df["Pool"] == POOL_BY_PREFIX["sc"],
        "NA",
    )


def main() -> None:
    args = parse_args()

    adf_dir = args.adf_dir.resolve()
    mapping_path = args.assay_mapping.resolve()

    if not adf_dir.exists() or not adf_dir.is_dir():
        raise NotADirectoryError(f"ADF directory not found: {adf_dir}")
    if not mapping_path.exists():
        raise FileNotFoundError(f"assay_mapping.csv not found: {mapping_path}")

    out_path = resolve_output_path(mapping_path=mapping_path, output_arg=args.output)
    if args.output is not None and out_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists: {out_path}. Use --overwrite to replace it.")

    mapping_df = read_table_auto(mapping_path)
    if "Name" not in mapping_df.columns:
        raise KeyError(f"Missing required column 'Name' in {mapping_path}")

    if "Pool" not in mapping_df.columns:
        mapping_df["Pool"] = ""

    tids_by_pool = load_tids_by_pool(adf_dir=adf_dir, projectname=args.projectname)
    name_to_pools = build_name_to_pools(tids_by_pool)

    def map_pool(name_value: str) -> str:
        name = str(name_value).strip()
        if name == "":
            return "No designs"

        labels = name_to_pools.get(name, [])
        if not labels:
            return "No designs"
        if len(labels) == 1:
            return labels[0]
        return "; ".join(labels)

    mapping_df["Pool"] = mapping_df["Name"].apply(map_pool)
    mapping_df["Singleton_PCR2_pool"] = build_singleton_pcr2_pool(mapping_df)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_df.to_csv(out_path, index=False, encoding="utf-8")

    n_total = len(mapping_df)
    n_no_design = int((mapping_df["Pool"] == "No designs").sum())
    n_p0 = int((mapping_df["Pool"] == POOL_BY_PREFIX["p0"]).sum())
    n_p2 = int((mapping_df["Pool"] == POOL_BY_PREFIX["p2"]).sum())
    n_sc = int((mapping_df["Pool"] == POOL_BY_PREFIX["sc"]).sum())
    n_multi = int(
        mapping_df["Pool"].astype(str).str.contains(";", regex=False).sum()
    )

    print(f"ADF directory: {adf_dir}")
    print(f"assay_mapping input: {mapping_path}")
    print(f"Output: {out_path}")
    print(f"Rows processed: {n_total}")
    print(f"Primary pool: {n_p0}")
    print(f"Secondary pool: {n_p2}")
    print(f"Singleton collection: {n_sc}")
    print(f"Multiple pool matches: {n_multi}")
    print(f"No designs: {n_no_design}")


if __name__ == "__main__":
    main()
