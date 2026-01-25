import pandas as pd
import argparse
import sys
import os
import csv

# ---------------------------------------------------------------------------------
# TAXONOMY PARSING UTILS
# ---------------------------------------------------------------------------------
def load_taxonomy_names(names_dmp_path):
    """
    Load names.dmp into a dict: taxid -> scientific_name
    names.dmp format: tax_id | name_txt | unique name | name_class |
    We only care about name_class == 'scientific name'
    """
    print(f"   Loading Taxonomy Names from {names_dmp_path}...")
    names_map = {}
    try:
        with open(names_dmp_path, 'r') as f:
            for line in f:
                parts = line.split('\t|\t')
                taxid = int(parts[0])
                name = parts[1].strip()
                # Fix parsing of the last column which has trailing delimiters
                name_class = parts[3].replace('\t|', '').strip()
                
                if name_class == 'scientific name':
                    names_map[taxid] = name
    except Exception as e:
        print(f"⚠️  Error loading names.dmp: {e}")
    return names_map

def load_merged_map(merged_dmp_path):
    """
    Load merged.dmp into a dict: old_taxid -> new_taxid
    merged.dmp format: old_tax_id | new_tax_id |
    """
    print(f"   Loading Merged TaxIDs from {merged_dmp_path}...")
    merged_map = {}
    try:
        if not os.path.exists(merged_dmp_path):
            print(f"⚠️  merged.dmp not found at {merged_dmp_path}. Skipping old taxid resolution.")
            return {}

        with open(merged_dmp_path, 'r') as f:
            for line in f:
                parts = line.split('\t|\t')
                old_id = int(parts[0])
                # second part might have trailing delimiters like '\t|\n'
                new_id = int(parts[1].replace('\t|', '').strip())
                merged_map[old_id] = new_id
    except Exception as e:
        print(f"⚠️  Error loading merged.dmp: {e}")
    return merged_map


def load_taxonomy_nodes(nodes_dmp_path):
    """
    Load nodes.dmp into a dict: taxid -> (parent_taxid, rank)
    nodes.dmp format: tax_id | parent_tax_id | rank | ...
    """
    print(f"   Loading Taxonomy Nodes from {nodes_dmp_path}...")
    nodes_map = {}
    # Key ranks we want to extract
    target_ranks = {'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'kingdom', 'domain'}
    
    try:
        with open(nodes_dmp_path, 'r') as f:
            for line in f:
                parts = line.split('\t|\t')
                taxid = int(parts[0])
                parent = int(parts[1])
                rank = parts[2].strip()
                
                nodes_map[taxid] = {'parent': parent, 'rank': rank}
    except Exception as e:
        print(f"⚠️  Error loading nodes.dmp: {e}")
    return nodes_map

def get_lineage(taxid, nodes_map, names_map):
    """
    Traverse the tree upwards to find ranks.
    Returns a dict with keys: kingdom, genus, species, etc.
    """
    lineage = {
        'superkingdom': 'Unknown',
        'domain': 'Unknown',
        'kingdom': 'Unknown',
        'phylum': 'Unknown',
        'class': 'Unknown',
        'order': 'Unknown',
        'family': 'Unknown',
        'genus': 'Unknown',
        'species': 'Unknown'
    }
    
    current_id = taxid
    
    # Explicit Ancestor Mapping (ID -> Name) to force consistent High-Level Grouping
    # 10239 = Viruses (superkingdom/root)
    # 2 = Bacteria (superkingdom)
    # 2157 = Archaea (superkingdom)
    # 2759 = Eukaryota (superkingdom/domain)
    known_roots = {
        10239: 'Viruses',
        2: 'Bacteria', 
        2157: 'Archaea',
        2759: 'Eukaryota'
    }

    # Safety: Max depth 100 to prevent infinite loops (though tree shouldn't have them)
    for _ in range(100):
        if current_id not in nodes_map:
            break
            
        # Check if we hit a known root
        if current_id in known_roots:
            lineage['superkingdom'] = known_roots[current_id]

        node = nodes_map[current_id]
        rank = node['rank']
        parent_id = node['parent']
        
        # Capture rank if it's one we assume
        if rank in lineage and lineage[rank] == 'Unknown':
            name = names_map.get(current_id, f"tid_{current_id}")
            lineage[rank] = name
            
        # Stop at root (usually 1, parent of 1 is 1)
        if current_id == 1 or current_id == parent_id:
            break
            
        current_id = parent_id
    
    # Consolidate Domain/Superkingdom
    if lineage['domain'] != 'Unknown':
        lineage['superkingdom'] = lineage['domain']
        
    return lineage

# ---------------------------------------------------------------------------------
# DATAFRAME PROCESSING
# ---------------------------------------------------------------------------------
def process_diamond_report(input_file, nodes_map, names_map, merged_map, output_prefix=None):
    cols = [
        "staxids", "sscinames", "qseqid", "sseqid", "pident", "length", 
        "mismatch", "gapopen", "qstart", "qend", "sstart", "send", 
        "evalue", "bitscore", "full_qseq"
    ]

    print(f"\nProcessing {input_file}...")
    try:
        # Some rows might have multiple taxids separated by ';'. taking the first one.
        df = pd.read_csv(input_file, sep="\t", names=cols, converters={'staxids': str})
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if df.empty:
        print("File is empty.")
        return

    # Process TaxIDs (handle '123;456' case by taking first)
    df['main_taxid'] = df['staxids'].apply(lambda x: int(x.split(';')[0]) if pd.notna(x) and x != '0' else 0)

    # ---------------------------------------------------------
    # 3. ADD TAXONOMY COLUMNS
    # ---------------------------------------------------------
    print("   Resolving Lineages (Domin/Kingdom, Genus, Species)...")
    
    # Cache lookups to speed up massive dataframes
    unique_taxids = df['main_taxid'].unique()
    tax_lookup = {}
    
    for tid in unique_taxids:
        if tid == 0:
            tax_lookup[tid] = {k: 'Unclassified' for k in ['superkingdom', 'genus', 'species']}
        else:
            final_tid = tid
            # Check if we need to merge
            if final_tid not in nodes_map and final_tid in merged_map:
                 final_tid = merged_map[final_tid]
                 
            tax_lookup[tid] = get_lineage(final_tid, nodes_map, names_map)
            
    # Map back to Dataframe
    # Use Superkingdom (which now includes Domain) as the "Kingdom" column for high-level grouping
    df['Kingdom'] = df['main_taxid'].map(lambda x: tax_lookup[x].get('superkingdom', 'Unknown'))
    df['Genus']   = df['main_taxid'].map(lambda x: tax_lookup[x].get('genus', 'Unknown'))
    df['Species'] = df['main_taxid'].map(lambda x: tax_lookup[x].get('species', 'Unknown'))


    # ---------------------------------------------------------
    # 4. CREATE SUMMARY (Grouped by Kingdom > Genus > Species)
    # ---------------------------------------------------------
    print("   Generating Summary...")
    group_cols = ['Kingdom', 'Genus', 'Species']
    summary = df.groupby(group_cols).agg({
        'qseqid': 'nunique',
        'pident': 'mean',
        'length': 'mean',
        'evalue': 'min'
    }).rename(columns={'qseqid': 'read_count', 'pident': 'avg_similarity'})
    
    summary = summary.sort_values(by=['read_count'], ascending=False)
    
    
    # ---------------------------------------------------------
    # 5. CREATE DETAILED RANKED OUTPUT
    # ---------------------------------------------------------
    print("   Sorting Details...")
    ranked_df = df.sort_values(by=['evalue', 'bitscore', 'length'], ascending=[True, False, False])
    
    out_cols = [
        "Kingdom", "Genus", "Species", "qseqid", "evalue", 
        "pident", "length", "sseqid", "full_qseq", "sscinames"
    ]
    ranked_df = ranked_df[out_cols]
    
    # ---------------------------------------------------------
    # 6. OUTPUT FILES
    # ---------------------------------------------------------
    if not output_prefix:
        output_prefix = os.path.splitext(input_file)[0]

    out_summary = f"{output_prefix}_summary_taxonomy.csv"
    out_ranked = f"{output_prefix}_ranked_taxonomy.csv"

    summary.to_csv(out_summary, float_format='%.2f')
    ranked_df.to_csv(out_ranked, index=False, float_format='%.2f')
    
    print(f"   ✅ Created: {out_summary}")
    print(f"   ✅ Created: {out_ranked}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize DIAMOND TSV with Taxonomy")
    parser.add_argument("input_files", nargs='+', help="Input .tsv files")
    parser.add_argument("--nodes", default="/home/ian/taxdmp/nodes.dmp", help="Path to nodes.dmp")
    parser.add_argument("--names", default="/home/ian/taxdmp/names.dmp", help="Path to names.dmp")
    parser.add_argument("--merged", default="/home/ian/taxdmp/merged.dmp", help="Path to merged.dmp")
    
    args = parser.parse_args()

    # Load Tax Dumps once
    nodes_map = load_taxonomy_nodes(args.nodes)
    names_map = load_taxonomy_names(args.names)
    merged_map = load_merged_map(args.merged)
    
    for f in args.input_files:
        process_diamond_report(f, nodes_map, names_map, merged_map)
