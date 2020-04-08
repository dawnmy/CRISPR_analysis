'''
File: CRT_parser.py
Created Date: January 27th 2020
Author: ZL Deng <dawnmsg(at)gmail.com>
---------------------------------------
Last Modified: 27th January 2020 9:52:22 pm
'''

import re
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@click.command()
@click.argument('crt', type=click.Path(exists=True))
@click.option('--inseq', type=click.Path(exists=True), default=None, help='Input fasta file with CRISPR regions [Required when outseq is enabled]')
@click.option('-o', '--output', type=str, required=True, help='The output file name')
@click.option('--outseq', type=str, default=None, help='Output fasta file of CRISPR regions and their down and up-stream sequences.')
def crispr_parser(crt, inseq, output, outseq):
    """convert the CRT output of CRISPR loci into tab-delimited file (similar to GFF/GTF)

    Arguments:
        crt {file path} -- The CRT output file
        inseq {file path} -- Input sequence in fasta format
        output {file path} -- The output tab file
        outseq {file path} -- Output fasta file with CRISPR regions and their flanking regions
    """
    #id_pattern = re.compile(r'^CRISPR\s[0-9]+')
    if outseq is not None: 
        if inseq is None:
            raise Exception('Please specify the fasta file for genome sequence if you need to output the sequences of CRISPR and flanking regions')
        else:
            output_fasta = True
    else:
        output_fasta = False

    id_range_pattern = re.compile(r'^CRISPR\s([0-9]+)\s+Range:\s([0-9]+)\s-\s([0-9]+)')
    found = False
    used_cols = [0, 4]
    crispr_regions = {}
    out_fh = open(output, 'w')
    with open(crt, 'r') as crt_fh:
        genome = crt_fh.readline().split()[1]
        for line in crt_fh:
            line = line.strip()
            if line.startswith('CRISPR '):

                #crispr_id = id_pattern.findall(line)[0].split()[1]
                crispr_id, range_start, range_stop = id_range_pattern.findall(line)[0]
                crispr_regions['CRISPR_' + crispr_id] = (int(range_start), int(range_stop))
                repeat_id = 0
                found = True
                size = None

            else:
                if found and line.split()[0].isdigit():
                    if size is None:
                        start, size = extract_by_idx(
                            line.split(), used_cols)
                        size = int(size.strip(","))
                    else:
                        start = line.split()[0]
                    start = int(start)

                    repeat_id += 1
                    out_line = '\t'.join([genome, "CRISPR", str(start), str(
                        start + size - 1), 'crispr_id={};repeat_id={}'.format(crispr_id, repeat_id)])
                    out_fh.write(out_line + "\n")
                elif found and next(crt_fh).startswith('Repeats:'):
                    found = False
                else:
                    continue
        out_fh.close()
    if output_fasta:
        #with open open(outseq, 'w') as seq_fh:
        extract_crispr_flanking(crispr_regions, inseq, outseq)
        


def extract_by_idx(l, idx):
    """extract the items of a list based on given indice

    Arguments:
        l {list} -- a list where you want to extract items
        idx {list} -- a list given the indice

    Returns:
        list -- a list of extracted items
    """
    return list(map(l.__getitem__, idx))


def extract_crispr_flanking(seq_region_dict, seq, outseq_file):
    """extract the CRISPR regions and their flanking regions from input fasta file
    
    Arguments:
        seq_region_dict {dict} -- The regions with start and end postion to extract
        seq {file path} -- The input sequence in fasta format
        outseq_file {file path} -- The output fasta file with extracted regions

    """
    seq_record = SeqIO.read(seq, 'fasta')
    gene_and_flanking = []
    for seq_name, region in seq_region_dict.items():
        start, end = region
        start = start - 1
        gene_region = seq_record.seq[start:end]
        upstream = seq_record.seq[start-1000:start]
        downstream = seq_record.seq[end:end+1000]
        gene_record = SeqRecord(gene_region, seq_name, description='CRISPR region')
        upstream_record = SeqRecord(upstream, '{}_upstream_1kb'.format(seq_name), description='CRISPR upstream 1kb')
        downstream_record = SeqRecord(downstream, '{}_downstream_1kb'.format(seq_name), description='CRISPR downstream 1kb')
        gene_and_flanking.extend([gene_record, upstream_record, downstream_record])

    SeqIO.write(gene_and_flanking, outseq_file, "fasta")


if __name__ == "__main__":
    crispr_parser()
