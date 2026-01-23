#!/bin/bash
set -e

# Log start
echo "Process started at $(date)" > /mnt/f/diamond_db_nr/filter_log.txt

# Pass 1
echo "Starting Pass 1 (chunk 00) at $(date)..." >> /mnt/f/diamond_db_nr/filter_log.txt
pigz -dc /home/ian/diamond_db/nr.gz | \
  /home/ian/miniconda3/bin/seqkit grep -f /home/ian/diamond_work/keep_seqids_part_00 --id-regexp '^([^\s]+)' -j 16 | \
  pigz -p 8 > /mnt/f/diamond_db_nr/nr_filtered_part00.faa.gz

# Pass 2
echo "Starting Pass 2 (chunk 01) at $(date)..." >> /mnt/f/diamond_db_nr/filter_log.txt
pigz -dc /home/ian/diamond_db/nr.gz | \
  /home/ian/miniconda3/bin/seqkit grep -f /home/ian/diamond_work/keep_seqids_part_01 --id-regexp '^([^\s]+)' -j 16 | \
  pigz -p 8 > /mnt/f/diamond_db_nr/nr_filtered_part01.faa.gz

# Pass 3
echo "Starting Pass 3 (chunk 02) at $(date)..." >> /mnt/f/diamond_db_nr/filter_log.txt
pigz -dc /home/ian/diamond_db/nr.gz | \
  /home/ian/miniconda3/bin/seqkit grep -f /home/ian/diamond_work/keep_seqids_part_02 --id-regexp '^([^\s]+)' -j 16 | \
  pigz -p 8 > /mnt/f/diamond_db_nr/nr_filtered_part02.faa.gz

# Merge
echo "Merging results at $(date)..." >> /mnt/f/diamond_db_nr/filter_log.txt
cat /mnt/f/diamond_db_nr/nr_filtered_part00.faa.gz \
    /mnt/f/diamond_db_nr/nr_filtered_part01.faa.gz \
    /mnt/f/diamond_db_nr/nr_filtered_part02.faa.gz > /mnt/f/diamond_db_nr/nr_filtered.faa.gz

# Clean up parts (optional, keeping for safety for now)
# rm /mnt/f/diamond_db_nr/nr_filtered_part*.faa.gz

echo "Completed successfully at $(date)!" >> /mnt/f/diamond_db_nr/filter_log.txt
