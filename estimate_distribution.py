#!/usr/bin/env python3
"""
Count accessions per taxonomic group by sampling prot.map.
Since prot.map is 151GB, we sample every Nth line.
"""
import sys
from collections import defaultdict

SAMPLE_RATE = 1000  # Sample every 1000th line

# Define roots
roots = {
    9606: "Human",
    2: "Bacteria", 
    2157: "Archaea",
    4751: "Fungi",
    10239: "Viruses"
}

# Load taxonomy tree
print("Loading nodes.dmp...", file=sys.stderr)
parent = {}
with open("/home/ian/taxdmp/nodes.dmp") as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]
        taxid, par = int(parts[0]), int(parts[1])
        parent[taxid] = par

# Cache: taxid -> root taxid
print("Building taxid->root cache...", file=sys.stderr)
taxid_cache = {}
def get_root(taxid):
    if taxid in taxid_cache:
        return taxid_cache[taxid]
    t = taxid
    visited = []
    found_root = None
    while t in parent and t not in taxid_cache:
        visited.append(t)
        if t in roots:
            found_root = t
            break
        t = parent[t]
    if t in taxid_cache:
        found_root = taxid_cache[t]
    for v in visited:
        taxid_cache[v] = found_root
    return found_root

# Load keep_taxids as a set
print("Loading keep_taxids.txt...", file=sys.stderr)
keep_taxids = set()
with open("/mnt/f/diamond_db_nr/keep_taxids.txt") as f:
    for line in f:
        keep_taxids.add(int(line.strip()))

# Sample prot.map and count
print(f"Sampling prot.map (1 in {SAMPLE_RATE})...", file=sys.stderr)
counts = defaultdict(int)
total_sampled = 0
kept_sampled = 0

with open("/mnt/f/diamond_db_nr/prot.map") as f:
    for i, line in enumerate(f):
        if i % SAMPLE_RATE != 0:
            continue
        if i % 10_000_000 == 0:
            print(f"  Processed {i:,} lines...", file=sys.stderr)
        
        total_sampled += 1
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            try:
                taxid = int(parts[1])
            except ValueError:
                continue  # Skip header or malformed
            
            if taxid in keep_taxids:
                kept_sampled += 1
                root = get_root(taxid)
                if root:
                    counts[roots[root]] += 1
                else:
                    counts["Unknown/Other"] += 1
            else:
                counts["Excluded"] += 1

print(f"\n=== Sampled {total_sampled:,} accessions (1 in {SAMPLE_RATE}) ===")
print(f"Kept: {kept_sampled:,} ({100*kept_sampled/total_sampled:.1f}%)")
print(f"\nEstimated distribution (extrapolated to full):")

# Estimate total from keep_seqids count
total_kept = 1_301_685_513  # From wc -l earlier
for name, count in sorted(counts.items(), key=lambda x: -x[1]):
    if name == "Excluded":
        continue
    estimated = int(count / kept_sampled * total_kept) if kept_sampled > 0 else 0
    pct = 100 * count / kept_sampled if kept_sampled > 0 else 0
    print(f"  {name}: ~{estimated:,} accessions ({pct:.1f}%)")

excluded = counts.get("Excluded", 0)
excluded_pct = 100 * excluded / total_sampled if total_sampled > 0 else 0
print(f"\nExcluded: ~{int(excluded/total_sampled * 2_600_000_000):,} accessions ({excluded_pct:.1f}%)")
