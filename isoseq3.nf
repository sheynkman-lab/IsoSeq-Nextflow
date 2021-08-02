log.info "Demultiplexing - N F  ~  version 0.1"
log.info "====================================="
// Header log info
log.info "\nPARAMETERS SUMMARY"
log.info "ccs_reads                 : ${params.ccs_reads}"
log.info "barcodes                  : ${params.barcodes}"

if (!params.ccs_reads) exit 1, "Cannot find file for parameter --ccs_reads: ${params.ccs_reads}"
if (params.ccs_reads.endsWith(".ccs.bam")){
    ch_ccs_reads = Channel.value(file(params.ccs_reads))
}
else{
    ch_ccs_reads = Channel.fromPath("${params.ccs_reads}/*.ccs.bam")
}

if (!params.barcodes) exit 1, "Cannot find file for parameter --barcodes: ${params.barcodes}"
ch_barcodes = Channel.value(file(params.barcodes))


_ch_all_barcodes = Channel
     .fromPath("${params.barcodes}")
     .splitFasta( record: [id: true, seqString: false ])
_ch_all_barcodes.into{
    _ch_all_barcodes_3p
    _ch_all_barcodes_5p
}
_ch_3prime_barcodes = _ch_all_barcodes_3p
    .filter { record -> record.id =~ /_3p$/ }
_ch_5prime_barcodes = _ch_all_barcodes_5p
    .filter { record -> record.id =~ /_5p$/ }

ch_barcode_pairs_list = _ch_3prime_barcodes.combine(_ch_5prime_barcodes)
/*
Generates the string representation of all valid barcodes
Valid parcodes are defined as a  3'--5' 
*/
process find_barcode_pairs{
    input:
        val(barcode) from ch_barcode_pairs_list
    output:
        val(barcode_string) into ch_barcode_pairs
    exec:
        barcode_string = barcode[0].id + "--" + barcode[1].id

}

ch_barcode_pairs.into{
    ch_barcode_pairs_view
    ch_barcode_pairs_use
}

if (params.genome_fasta.endsWith('.gz')){
   ch_genome_fasta = Channel.value(file(params.genome_fasta))
} else {
   ch_genome_fasta_uncompressed = Channel.value(file(params.genome_fasta))
}
if (params.genome_fasta.endsWith('.gz')) {
   process gunzip_gencome_fasta {
   tag "decompress gzipped genome fasta"
   cpus 1

   input:
   file(genome_fasta) from ch_genome_fasta

   output:
   file("*.{fa,fasta}") into ch_genome_fasta_uncompressed

   script:
   """
   gunzip -f ${genome_fasta}
   """
   }
}

/*
Demultiplexes and removes barcode and primers from reads
*/
process lima{
    tag "${ccs_read}"
    label "isoseq3"
    publishDir "${params.outdir}/isoseq3/lima/${ccs_read.name.split('\\.')[0]}", mode: "copy"
    input:
        file(ccs_read) from ch_ccs_reads.flatten()
        file(barcodes) from ch_barcodes

    output:
        file("*.bam") into ch_individual_lima_bam 
        file("*")

    script:
        sample_name = ccs_read.name.split('\\.')[0]
        """
        lima $ccs_read $barcodes ${sample_name}.bam --split-bam-named --isoseq --peek-guess
        """
}

/*
Gathers demultiplexed barcodes by their barcode pair.
*/
ch_individual_lima_bam
    .flatten()
    .map { file ->
        def key = file.name.toString().tokenize('.')[-2]
        return tuple(key, file)
     }
    .groupTuple()
    .set{ ch_lima_grouped_by_barcode }
/*
Merge the ccs reads together based on their barcode pairs
Filters out barcode pairs that are not valid
*/
process merge {
    tag "$barcode"
    publishDir "${params.outdir}/isoseq3/merge/${barcode}", mode: "copy"
    label "isoseq3"

    input:
        set barcode, file(bam_files) from ch_lima_grouped_by_barcode
        val barcode_pairs from ch_barcode_pairs_use.collect()
    when:
        barcode_pairs.contains(barcode)
    output:
        tuple val(barcode), file("*.bam") into ch_merged_reads
    script:
        """
        samtools merge ${barcode}.bam $bam_files
        """
}

/*
polyA tail trimming and concatemer removal
*/
process refine{
    tag "$barcode"
    publishDir "${params.outdir}/isoseq3/refine/${barcode}", mode: "copy"
    label "isoseq3"
    input:
        tuple val(barcode), file(merged_bam) from ch_merged_reads
        file(barcode_fasta) from ch_barcodes
    output:
        file("*.flnc.bam") into ch_refined_reads
        file("*")
    script:
        """
        isoseq3 refine --require-polya $merged_bam $barcode_fasta ${barcode}.flnc.bam
        """
}

/*
Hierarchical clustering of all reads
*/
process cluster {
    publishDir "${params.outdir}/isoseq3/cluster", mode: "copy"
    label "isoseq3"
    cpus params.max_cpus
    input: 
        file(refined_reads) from ch_refined_reads.collect()
    output:
        file("clustered.hq.bam") into ch_clustered_reads
        file("*")

    script:
        """
        ls $refined_reads > flnc.fofn
        isoseq3 cluster flnc.fofn clustered.bam --verbose --use-qvs
        """
}

/*
Align reads to genome
*/
process align {
    publishDir "${params.outdir}/isoseq3/align", mode: "copy"
    label "isoseq3"
    cpus params.max_cpus
    input: 
        file(clustered_reads) from ch_clustered_reads
        file(genome_fasta) from ch_genome_fasta_uncompressed
    output:
        file("*")
        file("${params.name}.aligned.bam") into ch_aligned_reads

    script:
        """
        pbmm2 align $genome_fasta $clustered_reads ${params.name}.aligned.bam --preset ISOSEQ --sort -j ${task.cpus} --log-level INFO
        """
}

/*
Collapse reads
*/
process collapse {
    publishDir "${params.outdir}/isoseq3/collapse", mode: "copy"
    label "isoseq3"
    cpus params.max_cpus
    input: 
        file(aligned_reads) from ch_aligned_reads
    output:
        file("*")

    script:
        """
        isoseq3 collapse $aligned_reads ${params.name}.collapsed.gff
        """
}

log.info "====================================="
log.info("Barcode pairs found")
ch_barcode_pairs_view.view()
log.info "====================================="