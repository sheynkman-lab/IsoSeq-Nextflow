# IsoSeq-Nextflow
## Nextflow implementation of IsoSeq3
---
## Steps
The IsoSeq pipeline is implemented in the following steps

1. lima
    * demultiplexing and primer removal 
2. merge
    * combine ccs reads based on barcode pair
3. refine
    * polyA tail trimming
    * concatemer removal
4. cluster
    * hierarchical, n*log(n) clusttering, alignment of shorter to longer sequences
    * iterative cluster merging
    * generate consensus for each read cluster using QV guided PoA
5. align
    * align reads to genome using pbmm2
6. collapse
    * collapse reads 


## Running the Pipeline
> nextflow run isoseq3.nf
--ccs_reads <ccs_reads> 
--barcodes <barcode> 
--genome_fasta <genome_fasta> 
--name example

## Input
1. --ccs_reads
    * ccs.bam file or directory containing cc.bam files
2. --barcodes
    * barcode primers to use when demultiplexing
    * 3' barcodes must end in _3p
    * 5' barcodes must end in _5p
3. --genome_fasta 
    * genome fasta file to use in alignment
4. --name
    *  name of experiment