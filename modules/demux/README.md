# Demux 
Demux the IsoSeq reads based on the barcode and the original file each read was derived from
- Uses output of collapse to determine read id and pbid
- Uses the output of refine to map barcode reads to collapsed reads
- Uses the individual reports of lima to build file-reads map


Tables produced
---

collapsed_reads_with_mappings.tsv
- collapsed read names, along with the barcodes and files each read maps to

demuxed_reads_by_barcode_and_file.tsv
- Table with following columns: pbid, barcode, filename, fl_count 

demuxed_pivot_by_barcode_and_file.tsv
- Multi-indexed column that pivots demuxed_reads_by_barcode_and_file.tsv by reads and barcode


demuxed_pivot_by_barcode.tsv
- demuxed reads pivotted by barcode
    - index : pbid
    - columns : barcode 
    - values : fl_count

demuxed_pivot_by_file.tsv
- demuxed reads pivotted by file
    - index : pbid
    - columns : filename
    - values : fl_count

Arguments
---
|Argument | Description |
|---------|-------------|
| -idir / --nextflow_results_directory | Directory used to save files from isoseq run. same as used in isoseq.nf
| -odir / --output_directory | Directory to save files
| --name | name of experiment. same as used in isoseq.nf |


Run 
---
```
python demux.py \
-idir ./results \
-odir ./results/demux \
--name example
```