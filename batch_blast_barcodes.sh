#!/bin/bash

# Batch DIAMOND BLASTx for Barcoded FASTQ Files
# OPTIMIZED VERSION: Viral Discovery Mode + Readable Taxonomy
#
# Changes:
# - Added `staxids sscinames sskingdoms skingdoms` to output format.
# - You can now read the .tsv file and see "Homo sapiens", "Metazoa", etc. directly.

INPUT_DIR="/mnt/e/Ian_Drive_Backup/Ian/Data/library_14_Ian_20240229/library_14_Ian_20240229/concatenated_barcodes"
OUTPUT_DIR="diamond_reports_viral"
DB="/home/ian/diamond_db/nr_filtered.dmnd"

# Check for WSL Mount accessibility
if [ ! -d "$INPUT_DIR" ]; then
    echo "❌ Error: Cannot find input directory in WSL."
    echo "   Path checked: $INPUT_DIR"
    echo "   Please run: sudo mount -t drvfs E: /mnt/e"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "🚀 Starting Viral-Optimized Batch BLASTx..."
echo "   Database: $DB"

count=0
for fastq in "$INPUT_DIR"/*.fastq; do
    [ -e "$fastq" ] || continue
    
    filename=$(basename "$fastq" .fastq)
    output_file="${OUTPUT_DIR}/${filename}_viral_report.tsv"
    
    echo "---------------------------------------------------"
    echo "Processing: $filename"
    
    # Run DIAMOND
    # Uses user-specified params: sensitive, max-target-seqs 2, evalue 1e-5
    # Added --threads 16 to maximize usage
    # Outputting staxids sscinames sskingdoms skingdoms for quick inspection
    # OPTIMIZED: Top 1 hit only, No full sequence output (save space)
    diamond blastx \
        -q "$fastq" \
        -d "$DB" \
        -o "$output_file" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms \
        --evalue 0.001 \
        --max-target-seqs 1 \
        --very-sensitive \
        --frameshift 15 \
        --threads 16
        
    if [ $? -eq 0 ]; then
        echo "✅ BLAST Done."
        
        # AUTOMATICALLY GENERATE CLEAN SUMMARY
        # This fixes the "Orthornavirae" -> "Unknown" issue by looking up the full lineage
        # and creating separate CSVs with "Kingdom" = "Viruses" etc.
        echo "   Running Taxonomy Summarizer..."
        python3 /home/ian/diamond_work/summarize_diamond.py \
            "$output_file" \
            --nodes /home/ian/taxdmp/nodes.dmp \
            --names /home/ian/taxdmp/names.dmp \
            --merged /home/ian/taxdmp/merged.dmp
            
    else
        echo "❌ BLAST failed for $filename"
    fi
    ((count++))
done

if [ "$count" -eq 0 ]; then
    echo "⚠️  No .fastq files found."
else
    echo "🎉 Processing complete! ($count files)"
fi
