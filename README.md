# Efficient DIAMOND Database Construction for Metagenomics
**Targeting: Viruses, Bacteria, Archaea, Fungi, and Human**

This repository contains the scripts and methodology used to filter the massive NCBI NR database (~186 GB compressed) into a smaller, specialized DIAMOND protein database (~12 GB) suitable for metagenomic analysis on consumer hardware.

## The Problem
Building a custom DIAMOND database from the full NR dataset typically requires significant RAM (often >64GB) when using tools like `seqkit` to filter by a list of ~1.3 billion accession IDs. On a standard machine (e.g., 48GB RAM), this leads to **Out-Of-Memory (OOM)** crashes.

Additionally, the standard `prot.accession2taxid.FULL` mapping file from NCBI (~151 GB) often lacks the specific header required by current versions of `diamond makedb`, causing "Invalid mapping file" errors.

## The Solution
This workflow implements two key strategies to overcome these limitations:

1.  **"Split-and-Filter" Strategy**: Instead of loading 1.3 billion accessions into RAM at once, we list them all but split the *list* into smaller chunks (e.g., 3 parts). We then run multiple low-memory passes of `seqkit grep` and merge the results.
    
2.  **Named Pipe for Taxonomy Mapping**: To fix the header issue without creating a duplicate 151 GB file, we use a named pipe (FIFO) to prepend the correct header (`accession.version\ttaxid`) to the mapping file stream on-the-fly.

## Prerequisites
- **DIAMOND** (v2.x)
- **SeqKit**
- **Python 3**
- Standard Linux utilities (`split`, `cat`, `mkfifo`, `wget`)
- ~300GB Disk Space (for intermediate files, though much is reclaimable)

## Workflow Steps

### 1. Generate Target Accession List
*Note: This step assumes you have `keep_taxids.txt` (a list of TaxIDs you want to keep).*

First, get the full accession-to-taxid mapping file from NCBI (`prot.accession2taxid.FULL.gz`). Then, filter it to get your "keep list" of accessions.

### 2. The "Split-and-Filter" Process
Use the `run_filter.sh` script. This script:
1.  Takes your huge `keep_seqids.txt`.
2.  Splits it into manageable chunks (adjust lines per chunk in the script if needed).
3.  Runs `seqkit grep` for each chunk against the massive `nr.gz`.
4.  Concatenates the outputs into `nr_filtered.faa.gz`.

```bash
bash run_filter.sh
```

### 3. Build the DIAMOND Database
We use a named pipe to feed the mapping file with the correct header.

1.  **Create the Header File:**
    ```bash
    printf "accession.version\ttaxid\n" > header.txt
    ```

2.  **Create the Pipe & Start the Build:**
    ```bash
    mkfifo prot_map.pipe
    
    # Start the stream in the background
    cat header.txt /path/to/prot.map > prot_map.pipe &
    
    # Run DIAMOND interactively reading from the pipe
    diamond makedb \
      --in nr_filtered.faa.gz \
      -d nr_filtered \
      --taxonmap prot_map.pipe \
      --taxonnodes nodes.dmp \
      --taxonnames names.dmp
    ```

## Included Scripts
- **`run_filter.sh`**: Automates the multi-pass `seqkit` filtering.
- **`audit_taxids.py`**: Checks distribution of TaxIDs across major kingdoms.
- **`check_unknown.py`**: Investigates TaxIDs that don't map to standard roots.
- **`estimate_distribution.py`**: fast statistical sampling of the map file.

## Expected Results
- **Input:** ~186 GB `nr.gz`
- **Output:** ~12 GB `nr_filtered.dmnd`
- **Content:** Comprehensive coverage of microbial and human sequences, with plants and invertebrates filtered out.
