#!/usr/bin/env python3
"""Audit TaxID distribution across major target taxonomic groups."""
import argparse
import os
import sys
from collections import defaultdict

ROOTS = {
    9606: "Human",
    2: "Bacteria",
    2157: "Archaea",
    4751: "Fungi",
    10239: "Viruses",
}


def load_parent_map(nodes_path):
    """Load NCBI nodes.dmp into a dict mapping taxid -> parent_taxid."""
    parent = {}
    with open(nodes_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]
            parent[int(parts[0])] = int(parts[1])
    return parent


def find_root(taxid, parent, roots=ROOTS):
    """Walk up the taxonomy tree and return the matching root TaxID, or None."""
    t = taxid
    visited = set()
    while t in parent and t not in visited:
        visited.add(t)
        if t in roots:
            return t
        nxt = parent[t]
        if nxt == t:
            break
        t = nxt
    return None


def audit_taxids(nodes_path, keep_taxids_path, roots=ROOTS):
    """Count TaxIDs in keep_taxids_path grouped by major taxonomic root."""
    parent = load_parent_map(nodes_path)
    counts = defaultdict(int)
    seen = set()

    with open(keep_taxids_path, "r", encoding="utf-8") as f:
        for line in f:
            raw = line.strip()
            if not raw:
                continue
            taxid = int(raw)
            if taxid in seen:
                continue
            seen.add(taxid)
            root = find_root(taxid, parent, roots)
            if root is not None:
                counts[roots[root]] += 1
            else:
                counts["Unknown/Other"] += 1

    return dict(counts)


def main(argv=None):
    taxdmp_dir = os.environ.get("TAXDMP_DIR", "taxdmp")
    parser = argparse.ArgumentParser(
        description="Audit TaxID distribution by taxonomic group."
    )
    parser.add_argument(
        "--nodes",
        default=os.path.join(taxdmp_dir, "nodes.dmp"),
        help="Path to NCBI nodes.dmp (default: $TAXDMP_DIR/nodes.dmp or taxdmp/nodes.dmp)",
    )
    parser.add_argument(
        "--keep-taxids",
        default="keep_taxids.txt",
        help="Path to file containing one TaxID per line (default: keep_taxids.txt)",
    )
    args = parser.parse_args(argv)

    print(f"Loading {args.nodes} and auditing {args.keep_taxids}...", file=sys.stderr)
    counts = audit_taxids(args.nodes, args.keep_taxids)

    print("\n=== Taxids per group ===")
    for name, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        print(f"{name}: {count:,} taxids")
    print(f"Total: {sum(counts.values()):,} taxids")


if __name__ == "__main__":
    main()
