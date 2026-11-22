# Efficient DIAMOND Database Construction for Metagenomics

**Target groups:** Viruses, Bacteria, Archaea, Fungi, and Human

Scripts and notes for filtering the full NCBI NR database (~186 GB compressed) into a smaller DIAMOND protein database (~12 GB) that runs comfortably on a regular workstation.

## Why this exists

Filtering full NR with `seqkit` against ~1.3 billion accession IDs usually takes more than 64 GB of RAM and crashes with OOM errors on a 48 GB machine.

On top of that, NCBI's `prot.accession2taxid.FULL` file (~151 GB) often does not have the exact header newer versions of `diamond makedb` expect, which causes "Invalid mapping file" errors.

## How it works

Two simple steps get around both limits:

1. **Split and filter:** Split the accession list into smaller chunks (for example, 3 parts). Run `seqkit grep` in separate low-memory passes, then concatenate the gzipped outputs.
2. **Named pipe for the taxonomy map:** Use a FIFO named pipe to prepend `accession.version\ttaxid` on the fly instead of duplicating a 151 GB file on disk.

## Requirements

- **DIAMOND** (v2.x)
- **SeqKit**
- **pigz**
- **Python 3** with `pandas` (`pip install -r requirements.txt`)
- Standard Linux tools (`split`, `cat`, `mkfifo`, `wget`)
- Around 300 GB of temporary disk space for intermediate files

## Workflow

### 1. Build and split your accession list

Start with `keep_taxids.txt` (one NCBI TaxID per line) and NCBI taxonomy dump files (`nodes.dmp`, `names.dmp`, `merged.dmp`).

Download `prot.accession2taxid.FULL.gz` from NCBI, extract the matching protein accessions into `keep_seqids.txt`, and split the list into chunks:

```bash
split -n l/3 -d keep_seqids.txt keep_seqids_part_
```

You can check your TaxID list and sample the map file beforehand:

```bash
python3 audit_taxids.py --nodes /path/to/taxdmp/nodes.dmp --keep-taxids keep_taxids.txt
python3 check_unknown.py --nodes /path/to/taxdmp/nodes.dmp --names /path/to/taxdmp/names.dmp --keep-taxids keep_taxids.txt
python3 estimate_distribution.py --nodes /path/to/taxdmp/nodes.dmp --keep-taxids keep_taxids.txt --prot-map /path/to/prot.map
```

### 2. Filter NR in chunks

Pass your paths via environment variables and run `run_filter.sh`:

```bash
NR_GZ=/path/to/nr.gz \
SEQIDS_PREFIX=keep_seqids_part_ \
OUT_DIR=filtered_nr \
THREADS=16 \
bash run_filter.sh
```

Set `CLEANUP_PARTS=1` to delete the per-chunk `.faa.gz` files after they are merged into `${OUT_DIR}/nr_filtered.faa.gz`.

### 3. Build the DIAMOND database

Feed the mapping file through a named pipe with the header prepended:

```bash
printf "accession.version\ttaxid\n" > header.txt
mkfifo prot_map.pipe

# Stream the header and mapping file in the background
cat header.txt /path/to/prot.map > prot_map.pipe &

# Build the DIAMOND database from the pipe
diamond makedb \
  --in filtered_nr/nr_filtered.faa.gz \
  -d nr_filtered \
  --taxonmap prot_map.pipe \
  --taxonnodes /path/to/taxdmp/nodes.dmp \
  --taxonnames /path/to/taxdmp/names.dmp
```

### 4. Run batch BLASTx and summarize taxonomy

Run `batch_blast_barcodes.sh` on a folder of `.fastq` or `.fastq.gz` files:

```bash
TAXDMP_DIR=/path/to/taxdmp \
THREADS=16 \
bash batch_blast_barcodes.sh /path/to/barcodes nr_filtered.dmnd diamond_reports_viral
```

Or summarize existing DIAMOND TSV files directly:

```bash
python3 summarize_diamond.py report.tsv \
  --nodes /path/to/taxdmp/nodes.dmp \
  --names /path/to/taxdmp/names.dmp \
  --merged /path/to/taxdmp/merged.dmp
```

## Scripts

- **`run_filter.sh`:** Runs multi-pass `seqkit grep` over split accession chunks and merges the result.
- **`batch_blast_barcodes.sh`:** Runs `diamond blastx` over barcoded FASTQ files and calls `summarize_diamond.py`.
- **`summarize_diamond.py`:** Resolves full NCBI lineages (including `merged.dmp` and domain/superkingdom ranks) and writes summary and ranked CSV tables.
- **`audit_taxids.py`:** Counts TaxIDs per target group.
- **`check_unknown.py`:** Lists TaxIDs in `keep_taxids.txt` that do not fall under the target roots.
- **`estimate_distribution.py`:** Samples `prot.map` to estimate kept and excluded accession counts.

## Testing

```bash
pip install -r requirements.txt pytest
python -m pytest tests -v
```

## Typical Output Size

- **Input:** ~186 GB `nr.gz`
- **Output:** ~12 GB `nr_filtered.dmnd`
- **Coverage:** Keeps microbial and human sequences while dropping plants and invertebrates.

## License

MIT. See [`LICENSE`](LICENSE).
