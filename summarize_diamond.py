#!/usr/bin/env python3
"""Summarize DIAMOND BLASTx tabular output with full NCBI taxonomy lineages."""
import argparse
import os
import sys

import pandas as pd


# ---------------------------------------------------------------------------------
# TAXONOMY PARSING UTILS
# ---------------------------------------------------------------------------------
def load_taxonomy_names(names_dmp_path):
    """
    Load names.dmp into a dict: taxid -> scientific_name.
    names.dmp format: tax_id | name_txt | unique name | name_class |
    """
    print(f"   Loading Taxonomy Names from {names_dmp_path}...")
    names_map = {}
    try:
        with open(names_dmp_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.split("\t|\t")
                if len(parts) < 4:
                    continue
                taxid = int(parts[0])
                name = parts[1].strip()
                name_class = parts[3].replace("\t|", "").strip()
                if name_class == "scientific name":
                    names_map[taxid] = name
    except Exception as e:
        print(f"Warning: error loading names.dmp: {e}", file=sys.stderr)
    return names_map


def load_merged_map(merged_dmp_path):
    """
    Load merged.dmp into a dict: old_taxid -> new_taxid.
    merged.dmp format: old_tax_id | new_tax_id |
    """
    print(f"   Loading Merged TaxIDs from {merged_dmp_path}...")
    merged_map = {}
    try:
        if not merged_dmp_path or not os.path.exists(merged_dmp_path):
            print(
                f"Warning: merged.dmp not found at {merged_dmp_path}. Skipping old taxid resolution.",
                file=sys.stderr,
            )
            return {}

        with open(merged_dmp_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.split("\t|\t")
                if len(parts) < 2:
                    continue
                old_id = int(parts[0])
                new_id = int(parts[1].replace("\t|", "").strip())
                merged_map[old_id] = new_id
    except Exception as e:
        print(f"Warning: error loading merged.dmp: {e}", file=sys.stderr)
    return merged_map


def load_taxonomy_nodes(nodes_dmp_path):
    """
    Load nodes.dmp into a dict: taxid -> {'parent': parent_taxid, 'rank': rank}.
    nodes.dmp format: tax_id | parent_tax_id | rank | ...
    """
    print(f"   Loading Taxonomy Nodes from {nodes_dmp_path}...")
    nodes_map = {}
    try:
        with open(nodes_dmp_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.split("\t|\t")
                if len(parts) < 3:
                    continue
                taxid = int(parts[0])
                parent = int(parts[1])
                rank = parts[2].strip()
                nodes_map[taxid] = {"parent": parent, "rank": rank}
    except Exception as e:
        print(f"Warning: error loading nodes.dmp: {e}", file=sys.stderr)
    return nodes_map


def parse_first_taxid(raw_value):
    """Parse the first integer TaxID from a DIAMOND staxids field (e.g. '123;456')."""
    if pd.isna(raw_value):
        return 0
    text = str(raw_value).strip()
    if not text or text == "0":
        return 0
    first = text.split(";")[0].strip()
    try:
        return int(first)
    except ValueError:
        return 0


def get_lineage(taxid, nodes_map, names_map):
    """
    Traverse the taxonomy tree upwards to populate standard taxonomic ranks.
    Returns a dict with keys: superkingdom, domain, kingdom, phylum, class,
    order, family, genus, species.
    """
    lineage = {
        "superkingdom": "Unknown",
        "domain": "Unknown",
        "kingdom": "Unknown",
        "phylum": "Unknown",
        "class": "Unknown",
        "order": "Unknown",
        "family": "Unknown",
        "genus": "Unknown",
        "species": "Unknown",
    }

    current_id = taxid

    # Explicit ancestor mapping (TaxID -> label) for consistent high-level grouping
    known_roots = {
        10239: "Viruses",
        2: "Bacteria",
        2157: "Archaea",
        2759: "Eukaryota",
    }

    for _ in range(100):
        if current_id not in nodes_map:
            break

        if current_id in known_roots:
            lineage["superkingdom"] = known_roots[current_id]

        node = nodes_map[current_id]
        rank = node["rank"]
        parent_id = node["parent"]

        if rank in lineage and lineage[rank] == "Unknown":
            lineage[rank] = names_map.get(current_id, f"tid_{current_id}")

        if current_id == 1 or current_id == parent_id:
            break

        current_id = parent_id

    # Consolidate Domain/Superkingdom (NCBI renamed superkingdom -> domain)
    if lineage["domain"] != "Unknown":
        lineage["superkingdom"] = lineage["domain"]

    return lineage


# ---------------------------------------------------------------------------------
# DATAFRAME PROCESSING
# ---------------------------------------------------------------------------------
def process_diamond_report(
    input_file, nodes_map, names_map, merged_map, output_prefix=None
):
    """Summarize a DIAMOND TSV report into taxonomy summary and ranked hit CSVs."""
    cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "staxids",
        "sscinames",
        "sskingdoms",
        "skingdoms",
    ]

    print(f"\nProcessing {input_file}...")
    try:
        df = pd.read_csv(input_file, sep="\t", names=cols, converters={"staxids": str})
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return None

    if df.empty:
        print("File is empty.")
        return None

    df["main_taxid"] = df["staxids"].apply(parse_first_taxid)

    print("   Resolving Lineages (Domain/Kingdom, Genus, Species)...")
    unique_taxids = df["main_taxid"].unique()
    tax_lookup = {}

    for tid in unique_taxids:
        if tid == 0:
            tax_lookup[tid] = {
                k: "Unclassified" for k in ("superkingdom", "genus", "species")
            }
        else:
            final_tid = tid
            if final_tid not in nodes_map and final_tid in merged_map:
                final_tid = merged_map[final_tid]
            tax_lookup[tid] = get_lineage(final_tid, nodes_map, names_map)

    df["Kingdom"] = df["main_taxid"].map(
        lambda x: tax_lookup[x].get("superkingdom", "Unknown")
    )
    df["Genus"] = df["main_taxid"].map(lambda x: tax_lookup[x].get("genus", "Unknown"))
    df["Species"] = df["main_taxid"].map(
        lambda x: tax_lookup[x].get("species", "Unknown")
    )

    print("   Generating Summary...")
    group_cols = ["Kingdom", "Genus", "Species"]
    summary = (
        df.groupby(group_cols)
        .agg(
            {
                "qseqid": "nunique",
                "pident": "mean",
                "length": "mean",
                "evalue": "min",
            }
        )
        .rename(columns={"qseqid": "read_count", "pident": "avg_similarity"})
        .sort_values(by=["read_count"], ascending=False)
    )

    print("   Sorting Details...")
    ranked_df = df.sort_values(
        by=["evalue", "bitscore", "length"], ascending=[True, False, False]
    )
    out_cols = [
        "Kingdom",
        "Genus",
        "Species",
        "qseqid",
        "evalue",
        "pident",
        "length",
        "sseqid",
        "sscinames",
    ]
    ranked_df = ranked_df[out_cols]

    if not output_prefix:
        output_prefix = os.path.splitext(input_file)[0]

    out_summary = f"{output_prefix}_summary_taxonomy.csv"
    out_ranked = f"{output_prefix}_ranked_taxonomy.csv"

    summary.to_csv(out_summary, float_format="%.6g")
    ranked_df.to_csv(out_ranked, index=False, float_format="%.6g")

    print(f"   Created: {out_summary}")
    print(f"   Created: {out_ranked}")
    return out_summary, out_ranked


def main(argv=None):
    taxdmp_dir = os.environ.get("TAXDMP_DIR", "taxdmp")
    parser = argparse.ArgumentParser(
        description="Summarize DIAMOND TSV reports with full NCBI taxonomy lineages."
    )
    parser.add_argument("input_files", nargs="+", help="Input DIAMOND .tsv files")
    parser.add_argument(
        "--nodes",
        default=os.path.join(taxdmp_dir, "nodes.dmp"),
        help="Path to nodes.dmp (default: $TAXDMP_DIR/nodes.dmp or taxdmp/nodes.dmp)",
    )
    parser.add_argument(
        "--names",
        default=os.path.join(taxdmp_dir, "names.dmp"),
        help="Path to names.dmp (default: $TAXDMP_DIR/names.dmp or taxdmp/names.dmp)",
    )
    parser.add_argument(
        "--merged",
        default=os.path.join(taxdmp_dir, "merged.dmp"),
        help="Path to merged.dmp (default: $TAXDMP_DIR/merged.dmp or taxdmp/merged.dmp)",
    )
    args = parser.parse_args(argv)

    nodes_map = load_taxonomy_nodes(args.nodes)
    names_map = load_taxonomy_names(args.names)
    merged_map = load_merged_map(args.merged)

    for input_file in args.input_files:
        process_diamond_report(input_file, nodes_map, names_map, merged_map)


if __name__ == "__main__":
    main()
