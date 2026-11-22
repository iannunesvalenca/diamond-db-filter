#!/usr/bin/env python3
"""
Estimate accession counts per taxonomic group by sampling a prot.accession2taxid map.
Since full NCBI mapping files can exceed 150 GB, this script samples every Nth line.
"""
import argparse
import os
import sys
from collections import defaultdict

from audit_taxids import ROOTS, load_parent_map


def estimate_distribution(
    nodes_path,
    keep_taxids_path,
    prot_map_path,
    sample_rate=1000,
    roots=ROOTS,
):
    """Sample prot_map_path and return (counts, total_sampled, kept_sampled)."""
    if sample_rate < 1:
        raise ValueError("sample_rate must be >= 1")

    parent = load_parent_map(nodes_path)
    taxid_cache = {}

    def get_root(taxid):
        if taxid in taxid_cache:
            return taxid_cache[taxid]
        t = taxid
        visited = []
        visited_set = set()
        found_root = None
        while t in parent and t not in taxid_cache and t not in visited_set:
            visited.append(t)
            visited_set.add(t)
            if t in roots:
                found_root = t
                break
            nxt = parent[t]
            if nxt == t:
                break
            t = nxt
        if t in taxid_cache:
            found_root = taxid_cache[t]
        for v in visited:
            taxid_cache[v] = found_root
        return found_root

    keep_taxids = set()
    with open(keep_taxids_path, "r", encoding="utf-8") as f:
        for line in f:
            raw = line.strip()
            if raw:
                keep_taxids.add(int(raw))

    counts = defaultdict(int)
    total_sampled = 0
    kept_sampled = 0

    with open(prot_map_path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i % sample_rate != 0:
                continue
            if i > 0 and i % 10_000_000 == 0:
                print(f"  Processed {i:,} lines...", file=sys.stderr)

            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            try:
                taxid = int(parts[1])
            except ValueError:
                continue  # Skip header or malformed rows

            total_sampled += 1
            if taxid in keep_taxids:
                kept_sampled += 1
                root = get_root(taxid)
                if root is not None:
                    counts[roots[root]] += 1
                else:
                    counts["Unknown/Other"] += 1
            else:
                counts["Excluded"] += 1

    return dict(counts), total_sampled, kept_sampled


def main(argv=None):
    taxdmp_dir = os.environ.get("TAXDMP_DIR", "taxdmp")
    parser = argparse.ArgumentParser(
        description="Estimate kept vs. excluded accession distribution by sampling prot.map."
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
    parser.add_argument(
        "--prot-map",
        default="prot.map",
        help="Path to two-column accession<TAB>taxid mapping file (default: prot.map)",
    )
    parser.add_argument(
        "--sample-rate",
        type=int,
        default=1000,
        help="Sample every Nth line (default: 1000)",
    )
    args = parser.parse_args(argv)

    print(
        f"Sampling {args.prot_map} (1 in {args.sample_rate})...",
        file=sys.stderr,
    )
    counts, total_sampled, kept_sampled = estimate_distribution(
        args.nodes, args.keep_taxids, args.prot_map, sample_rate=args.sample_rate
    )

    print(f"\n=== Sampled {total_sampled:,} accessions (1 in {args.sample_rate}) ===")
    if total_sampled == 0:
        print("No valid accession rows sampled.")
        return

    print(f"Kept: {kept_sampled:,} ({100 * kept_sampled / total_sampled:.1f}%)")
    print("\nEstimated distribution (extrapolated by sample rate):")
    for name, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        if name == "Excluded":
            continue
        estimated = count * args.sample_rate
        pct = 100 * count / kept_sampled if kept_sampled > 0 else 0.0
        print(f"  {name}: ~{estimated:,} accessions ({pct:.1f}%)")

    excluded = counts.get("Excluded", 0)
    excluded_pct = 100 * excluded / total_sampled
    print(f"\nExcluded: ~{excluded * args.sample_rate:,} accessions ({excluded_pct:.1f}%)")


if __name__ == "__main__":
    main()
