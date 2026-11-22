# Efficient DIAMOND Database Construction for Metagenomics
**Targeting: Viruses, Bacteria, Archaea, Fungi, and Human**

This repository contains the scripts and methodology used to filter the massive NCBI NR database (~186 GB compressed) into a smaller, specialized DIAMOND protein database (~12 GB) suitable for metagenomic analysis on consumer hardware.

## The Problem
Building a custom DIAMOND database from the full NR dataset typically requires significant RAM (often >64 GB) when using tools like `seqkit` to filter by a list of ~1.3 billion accession IDs. On a standard workstation (e.g., 48 GB RAM), this leads to **Out-Of-Memory (OOM)** crashes.

Additionally, the standard `prot.accession2taxid.FULL` mapping file from NCBI (~151 GB) often lacks the specific header required by current versions of `diamond makedb`, causing "Invalid mapping file" errors.

## The Solution
This workflow implements two key strategies to overcome these limitations:

1. **"Split-and-Filter" Strategy**: Instead of loading 1.3 billion accessions into RAM at once, split the accession list into smaller chunks (e.g., 3 parts), run multiple low-memory passes of `seqkit grep`, and concatenate the gzipped outputs.
2. **Named Pipe for Taxonomy Mapping**: To fix the header issue without duplicating a 151 GB file, use a named pipe (FIFO) to prepend the required header (`accession.version\ttaxid`) on-the-fly.

## Prerequisites
- **DIAMOND** (v2.x)
- **SeqKit**
- **pigz**
- **Python 3** with dependencies from `requirements.txt` (`pip install -r requirements.txt`)
- Standard Linux utilities (`split`, `cat`, `mkfifo`, `wget`)
- ~300 GB disk space (for intermediate files, most of which are reclaimable)

## Workflow Steps

### 1. Generate Target Accession List
*Assumes you have `keep_taxids.txt` (one NCBI TaxID per line) and NCBI taxonomy dump files (`nodes.dmp`, `names.dmp`, `merged.dmp`).*

Download `prot.accession2taxid.FULL.gz` from NCBI and extract the matching protein accessions into `keep_seqids.txt`, then split it into chunks:

```bash
split -n l/3 -d keep_seqids.txt keep_seqids_part_
```

You can inspect your TaxID list and estimate the retained accession distribution beforehand:

```bash
python3 audit_taxids.py --nodes /path/to/taxdmp/nodes.dmp --keep-taxids keep_taxids.txt
python3 check_unknown.py --nodes /path/to/taxdmp/nodes.dmp --names /path/to/taxdmp/names.dmp --keep-taxids keep_taxids.txt
python3 estimate_distribution.py --nodes /path/to/taxdmp/nodes.dmp --keep-taxids keep_taxids.txt --prot-map /path/to/prot.map
```

### 2. The "Split-and-Filter" Process
Configure input/output paths via environment variables and run `run_filter.sh`:

```bash
NR_GZ=/path/to/nr.gz \
SEQIDS_PREFIX=keep_seqids_part_ \
OUT_DIR=filtered_nr \
THREADS=16 \
bash run_filter.sh
```

Set `CLEANUP_PARTS=1` to delete per-chunk `.faa.gz` files after merging into `${OUT_DIR}/nr_filtered.faa.gz`.

### 3. Build the DIAMOND Database
Use a named pipe to feed the mapping file with the required header:

```bash
printf "accession.version\ttaxid\n" > header.txt
mkfifo prot_map.pipe

# Start the stream in the background
cat header.txt /path/to/prot.map > prot_map.pipe &

# Build the database reading taxonomy mapping from the pipe
diamond makedb \
  --in filtered_nr/nr_filtered.faa.gz \
  -d nr_filtered \
  --taxonmap prot_map.pipe \
  --taxonnodes /path/to/taxdmp/nodes.dmp \
  --taxonnames /path/to/taxdmp/names.dmp
```

### 4. Batch BLASTx & Taxonomy Summarization
Run `batch_blast_barcodes.sh` on a directory of `.fastq` or `.fastq.gz` files:

```bash
TAXDMP_DIR=/path/to/taxdmp \
THREADS=16 \
bash batch_blast_barcodes.sh /path/to/barcodes nr_filtered.dmnd diamond_reports_viral
```

Or run `summarize_diamond.py` directly on existing DIAMOND TSV reports:

```bash
python3 summarize_diamond.py report.tsv \
  --nodes /path/to/taxdmp/nodes.dmp \
  --names /path/to/taxdmp/names.dmp \
  --merged /path/to/taxdmp/merged.dmp
```

## Included Scripts
- **`run_filter.sh`**: Multi-pass `seqkit grep` filtering over split accession chunks.
- **`batch_blast_barcodes.sh`**: Runs `diamond blastx` over barcoded FASTQ files and invokes `summarize_diamond.py`.
- **`summarize_diamond.py`**: Resolves full NCBI taxonomy lineages (including `merged.dmp` and domain/superkingdom mapping) and outputs per-taxon summaries and ranked hits.
- **`audit_taxids.py`**: Counts TaxIDs across major target taxonomic groups.
- **`check_unknown.py`**: Lists TaxIDs in `keep_taxids.txt` that do not descend from the target roots.
- **`estimate_distribution.py`**: Estimates kept vs. excluded accession counts by sampling `prot.map`.

## Testing

```bash
pip install -r requirements.txt pytest
python -m pytest tests -v
```

## Expected Results
- **Input:** ~186 GB `nr.gz`
- **Output:** ~12 GB `nr_filtered.dmnd`
- **Content:** Comprehensive coverage of microbial and human sequences, with plants and invertebrates filtered out.

## License
MIT — see [`LICENSE`](LICENSE).
