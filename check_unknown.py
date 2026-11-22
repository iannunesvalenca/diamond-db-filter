#!/usr/bin/env python3
"""Inspect TaxIDs in keep_taxids.txt that do not map to target roots."""
import argparse
import os
import sys

from audit_taxids import ROOTS, find_root, load_parent_map


def load_scientific_names(names_path):
    """Load NCBI names.dmp scientific names into a dict mapping taxid -> name."""
    names = {}
    with open(names_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[int(parts[0])] = parts[1]
    return names


def find_unknown_taxids(nodes_path, keep_taxids_path, roots=ROOTS):
    """Return TaxIDs from keep_taxids_path that do not descend from any target root."""
    parent = load_parent_map(nodes_path)
    unknown = []
    with open(keep_taxids_path, "r", encoding="utf-8") as f:
        for line in f:
            raw = line.strip()
            if not raw:
                continue
            taxid = int(raw)
            if find_root(taxid, parent, roots) is None:
                unknown.append(taxid)
    return unknown


def main(argv=None):
    taxdmp_dir = os.environ.get("TAXDMP_DIR", "taxdmp")
    parser = argparse.ArgumentParser(
        description="Inspect TaxIDs that do not map to the target taxonomic roots."
    )
    parser.add_argument(
        "--nodes",
        default=os.path.join(taxdmp_dir, "nodes.dmp"),
        help="Path to NCBI nodes.dmp (default: $TAXDMP_DIR/nodes.dmp or taxdmp/nodes.dmp)",
    )
    parser.add_argument(
        "--names",
        default=os.path.join(taxdmp_dir, "names.dmp"),
        help="Path to NCBI names.dmp (default: $TAXDMP_DIR/names.dmp or taxdmp/names.dmp)",
    )
    parser.add_argument(
        "--keep-taxids",
        default="keep_taxids.txt",
        help="Path to file containing one TaxID per line (default: keep_taxids.txt)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Maximum number of unknown TaxIDs to display (default: 20)",
    )
    args = parser.parse_args(argv)

    print(f"Loading {args.names}...", file=sys.stderr)
    names = load_scientific_names(args.names)
    print(f"Checking {args.keep_taxids} against {args.nodes}...", file=sys.stderr)
    unknown = find_unknown_taxids(args.nodes, args.keep_taxids)

    print(f"\nFound {len(unknown)} unknown taxids:")
    for t in unknown[: args.limit]:
        name = names.get(t, "NO NAME")
        print(f"  {t}: {name}")
    if len(unknown) > args.limit:
        print(f"  ... and {len(unknown) - args.limit} more")


if __name__ == "__main__":
    main()
