# IsoSeq-Nextflow
## Nextflow implementation of IsoSeq3
---
## Steps
The IsoSeq pipeline is implemented in the following steps

1. lima
    * demultiplexing and primer removal 
2. merge
    * combine ccs reads based on barcode pair
    * only valid barcodes pairs are merged - all combinations of 3p__5p barcode pairs as defined in barcodes.fasta
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


## Running the Pipeline with Docker
> nextflow run isoseq3.nf
--ccs_reads ccs/read/dir
--barcodes barcodes.fasta
--genome_fasta genome.fasta
--name example

## Running the Pipeline with Singularity
> nextflow run isoseq3.nf
--ccs_reads ccs/read/dir 
--barcodes barcode.fasta
--genome_fasta genome.fasta
--name example
-with-singularity
-without-docker
 
## Input
1. --ccs_reads
    * ccs.bam file or directory containing ccs.bam files
2. --barcodes
    * barcode primers to use when demultiplexing
    * 3' barcodes must end in _3p
    * 5' barcodes must end in _5p
3. --genome_fasta 
    * genome fasta file to use in alignment
4. --name
    *  name of experiment

![IsoSeq Workflow](https://github.com/sheynkman-lab/IsoSeq-Nextflow/blob/main/media/isoseq-workflow.svg)
