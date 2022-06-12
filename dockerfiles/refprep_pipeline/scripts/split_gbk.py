#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser(description='Split GBK file with multiple contigs into separate GBK files with individual contigs.')
parser.add_argument('-g', '--gbk', help='GBK file that needs to be splitted')
parser.add_argument('-o', '--out', help='Output folder')
args = parser.parse_args()

# create output directory
if not os.path.exists(args.out):
    try:
        os.mkdir(args.out)
    except OSError:
        print(f'was not able to create {args.out}')
    else:
        print(f'{args.out} was successfully created')


# split GBK and output only contigs with genes
with open(args.gbk, 'r') as gbk_file:
    for record in SeqIO.parse(gbk_file, 'gb'):
        cds = [1 for feature in record.features if feature.type == 'CDS']
        if len(cds) > 0:
            print(record.id)
            # BioPython follows NCBI standard guidelines for GBK and has a very strict limit on the
            # Locus ID length which is violated by most assembler programs.
            original_id = record.id
            truncated_id = original_id[:10] + 't'
            record.id = truncated_id
            record.name = truncated_id
            with open(os.path.join(args.out, truncated_id + '_truncated.gbk'), 'w') as gbk_out_truncared:
                SeqIO.write(record, gbk_out_truncared, 'gb')

            # Fix truncated ID
            with open(os.path.join(args.out, truncated_id + '_truncated.gbk'), 'r') as gbk_out_truncated:
                gbk_content = gbk_out_truncated.read()
                gbk_content = gbk_content.replace(truncated_id, original_id)
                with open(os.path.join(args.out, original_id + '.gbk'), 'w') as gbk_out:
                    gbk_out.write(gbk_content)

            try:
                os.remove(os.path.join(args.out, truncated_id + '_truncated.gbk'))
            except OSError:
                print("Could not remove file " + os.path.join(args.out, truncated_id + '_truncated.gbk\nSomething wrong'))

