#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import logging


def main():
    parser = argparse.ArgumentParser(description='Add locus tags and gene names to gff file')
    parser.add_argument('-f', '--fna', help='Input FNA file to build snpEff config portion.')
    parser.add_argument('-n', '--name', help='Genome name.')
    parser.add_argument('-o', '--output', help='Output file.')
    args = parser.parse_args()

    fna_file_path = args.fna
    genome_name = args.name
    genome_id = genome_name.replace(' ', '_')
    output_file_path = args.output
    contig_ids = []

    # Read contigs IDs
    logging.info('Read contig IDs file.')
    with open(fna_file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            contig_ids.append(record.id)

    # Write snpEff config part for the genome
    logging.info('Write snpEff config part for the genome.')
    with open(output_file_path, 'w') as output_file:
        output_file.write(f'# {genome_name}\n')
        output_file.write(f'{genome_id}.genome: {genome_name}\n')
        output_file.write(f'    {genome_id}.chromosomes: {", ".join(contig_ids)}\n')
        for contig in contig_ids:
            output_file.write(f'    {genome_id}.{contig}.codonTable: Bacterial_and_Plant_Plastid\n')


if __name__ == "__main__":
    main()
