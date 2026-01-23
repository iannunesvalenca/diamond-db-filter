#!/usr/bin/env python3
"""Audit taxid distribution by taxonomic group."""
import sys
from collections import defaultdict

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

# Build taxid->root mapping for all taxids in keep_taxids.txt
print("Mapping taxids to roots...", file=sys.stderr)
taxid_to_root = {}
with open("/mnt/f/diamond_db_nr/keep_taxids.txt") as f:
    for line in f:
        taxid = int(line.strip())
        # Walk up to find root
        t = taxid
        found_root = None
        visited = set()
        while t in parent and t not in visited:
            visited.add(t)
            if t in roots:
                found_root = t
                break
            t = parent[t]
        taxid_to_root[taxid] = found_root

# Count per root
print("Counting taxids per root...", file=sys.stderr)
counts = defaultdict(int)
for taxid, root in taxid_to_root.items():
    if root:
        counts[roots[root]] += 1
    else:
        counts["Unknown/Other"] += 1

print("\n=== Taxids per group ===")
for name, count in sorted(counts.items(), key=lambda x: -x[1]):
    print(f"{name}: {count:,} taxids")
print(f"Total: {sum(counts.values()):,} taxids")
