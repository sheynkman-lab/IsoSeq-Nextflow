import argparse

import pandas as pd 
import os
import logging

"""
Demux the IsoSeq reads based on the barcode and the original file each read was derived from
Uses the output of refine to map barcode reads to collapsed reads
Uses the individual reports of lima to build file-reads map

Tables produced
---------------

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
"""

def read_collapsed_reads(nextflow_results_dir, name):
    
    collapsed_ifile = os.path.join(nextflow_results_dir, f'isoseq3/collapse/{name}.collapsed.read_stat.txt')
    logging.info(f'Reading collapsed reads...\n{collapsed_ifile}')
    collapsed_reads = pd.read_table(collapsed_ifile)
    logging.info('Finished reading collapsed reads')
    return collapsed_reads

def build_reads_file_map(nextflow_results_dir):
    logging.info('Building reads file map from lima reports')

    lima_dir = os.path.join(nextflow_results_dir, 'isoseq3/lima')
    read_file_data = []

    for dir in os.listdir(lima_dir):
        logging.info(f'Reading lima report for {dir}')
        lima_df = pd.read_table(os.path.join(lima_dir, dir, f'{dir}.lima.report'), usecols=['ZMW'])
        lima_df['id'] = lima_df['ZMW'].apply(lambda x: x + '/ccs')
        lima_df['file'] = dir
        read_file_data.append(lima_df)
    read_file_data = pd.concat(read_file_data)
    reads_file_map = dict(zip(read_file_data.id, read_file_data.file))
    logging.info('Finished building reads file map')
    return reads_file_map


def build_reads_barcode_map(nextflow_results_dir):

    refined_full_report_ifile = os.path.join(nextflow_results_dir, 'isoseq3/refine/full.flnc.report.csv')
    logging.info(f'Building reads barcode map\n{refined_full_report_ifile}')

    refined_data = pd.read_csv(refined_full_report_ifile, usecols=['id', 'primer'])
    reads_barcode_map = dict(zip(refined_data.id, refined_data.primer))
    logging.info('Finished building reads barcode map')
    return reads_barcode_map

def add_barcode_and_file(collapsed_reads, reads_barcode_map, reads_file_map):
    logging.info('Mapping reads to associated barcodes and files')
    collapsed_reads['barcode'] = collapsed_reads['id'].map(reads_barcode_map)
    collapsed_reads['file'] = collapsed_reads['id'].map(reads_file_map)
    logging.info('Finished mapping reads to associated barcodes and files')
    return collapsed_reads

def demux_reads(collapsed_reads):
    logging.info('Demuxing reads by barcode and file')
    demuxed_reads = (
        collapsed_reads
            .groupby(['pbid', 'barcode', 'file'])
            .size()
            .reset_index(name='fl_count')
    )
    return demuxed_reads

def pivot_demuxed_reads(demuxed_reads):
    logging.info('Pivot of demuxed reads by barcode and file')
    demuxed_pivot = (
        demuxed_reads
        .pivot(index='pbid', columns=['barcode', 'file'], values='fl_count')
        .fillna(0)
        .astype(int)
    )
    return demuxed_pivot

def pivot_reads_by_barcode(collapsed_reads):
    logging.info('Pivot of reads by barcode')
    demuxed_reads_barcode = (
        collapsed_reads
            .groupby(['pbid', 'barcode'])
            .size()
            .reset_index(name='fl_count')
    )

    demuxed_pivot_barcode = (
        demuxed_reads_barcode
        .pivot(index='pbid', columns='barcode', values='fl_count')
        .fillna(0)
        .astype(int)
    )
    return demuxed_pivot_barcode

def pivot_reads_by_file(collapsed_reads):
    logging.info('Pivot of reads by file')
    demuxed_reads_file = (
        collapsed_reads
            .groupby(['pbid', 'file'])
            .size()
            .reset_index(name='fl_count')
    )

    demuxed_pivot_file = (
        demuxed_reads_file
        .pivot(index='pbid', columns='file', values='fl_count')
        .fillna(0)
        .astype(int)
    )
    return demuxed_pivot_file


def save(collapsed_reads,demuxed_reads,demuxed_pivot, demuxed_pivot_barcode, demuxed_pivot_file, name, odir):
    logging.info('Saving files...')
    collapsed_reads.to_csv(os.path.join(odir, f'{name}.collapsed_reads_with_mapping.tsv'), sep='\t', index=False)
    demuxed_reads.to_csv(os.path.join(odir, f'{name}.demuxed_reads_by_barcode_and_file.tsv'), sep='\t', index=False)
    demuxed_pivot.to_csv(os.path.join(odir, f'{name}.demuxed_pivot_by_barcode_and_file.tsv'), sep='\t')
    demuxed_pivot_barcode.to_csv(os.path.join(odir, f'{name}.demuxed_pivot_by_barcode.tsv'), sep='\t')
    demuxed_pivot_file.to_csv(os.path.join(odir, f'{name}.demuxed_pivot_by_file.tsv'), sep='\t')

    

    




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-idir', '--nextflow_results_directory', help='output directory of Isoseq-Nextflow')
    parser.add_argument('-odir', '--output_directory')
    parser.add_argument('--name', help='Name used during Isoseq-Nextflow')
    args = parser.parse_args()
    nextflow_results_dir = args.nextflow_results_directory
    output_directory = args.output_directory
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    logging.basicConfig(filename=os.path.join(output_directory,'demux.log'),  level=logging.INFO)

    collapsed_reads = read_collapsed_reads(nextflow_results_dir, args.name)
    reads_file_map = build_reads_file_map(nextflow_results_dir)
    reads_barcode_map = build_reads_barcode_map(nextflow_results_dir)
    collapsed_reads = add_barcode_and_file(collapsed_reads, reads_barcode_map, reads_file_map)
    demuxed_reads = demux_reads(collapsed_reads)
    demuxed_pivot = pivot_demuxed_reads(demuxed_reads)
    demuxed_pivot_barcode = pivot_reads_by_barcode(collapsed_reads)
    demuxed_pivot_file = pivot_reads_by_file(collapsed_reads)
    save(collapsed_reads, demuxed_reads, demuxed_pivot, demuxed_pivot_barcode, demuxed_pivot_file, args.name, output_directory)



if __name__ == '__main__':
    main()
