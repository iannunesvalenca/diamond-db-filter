#!/usr/bin/env python3
"""Check what the Unknown/Other taxids are."""
import sys

# Define roots
roots = {
    9606: "Human",
    2: "Bacteria", 
    2157: "Archaea",
    4751: "Fungi",
    10239: "Viruses"
}

# Load taxonomy tree and names
print("Loading nodes.dmp...", file=sys.stderr)
parent = {}
with open("/home/ian/taxdmp/nodes.dmp") as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]
        taxid, par = int(parts[0]), int(parts[1])
        parent[taxid] = par

print("Loading names.dmp...", file=sys.stderr)
names = {}
with open("/home/ian/taxdmp/names.dmp") as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]
        if parts[3] == "scientific name":
            names[int(parts[0])] = parts[1]

# Find Unknown/Other taxids
print("\n=== Unknown/Other taxids ===", file=sys.stderr)
unknown = []
with open("/mnt/f/diamond_db_nr/keep_taxids.txt") as f:
    for line in f:
        taxid = int(line.strip())
        t = taxid
        found_root = None
        visited = set()
        while t in parent and t not in visited:
            visited.add(t)
            if t in roots:
                found_root = t
                break
            t = parent[t]
        if not found_root:
            unknown.append(taxid)

print(f"\nFound {len(unknown)} unknown taxids:")
for t in unknown[:20]:
    name = names.get(t, "NO NAME")
    print(f"  {t}: {name}")
if len(unknown) > 20:
    print(f"  ... and {len(unknown)-20} more")
