#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Adds locus_tag field to the GBK file and sort features by position.')
parser.add_argument('-g', '--gbk', help='Input GBK file.')
parser.add_argument('-l', '--lt_file', default='', help='File with locus tag information. Required coulmns:'
                                                        'accession, start, end, strand, refseq_locus_tag, gene.')
parser.add_argument('--sort_only', action='store_true', help='If selectd only sorting will be done.')
parser.add_argument('-o', '--out_gbk', help='Output GBK file.')
args = parser.parse_args()


def sort_coordinates(accession, start, end):
    if int(end) >= int(start):
        return accession + '_' + start + '_' + end
    else:
        return accession + '_' + end + '_' + start


# Collect new annotations in a dictionary
if not args.sort_only:
    annotations = pd.read_csv(args.lt_file, sep='\t')
    annotations = annotations.fillna('-')

    feature_ids = [sort_coordinates(accession, start, end) for accession, start, end in
                   zip(
                       annotations['accession'].astype(str).to_list(),
                       annotations['start'].astype(str).to_list(),
                       annotations['end'].astype(str).to_list()
                   )]
                   
    # sort coordinates as in BioPython record coordinates are in ascending order

    # In features save locus_tag and gene.
    features_new = [[lt, gene] for lt, gene in
                    zip(
                        annotations['refseq_locus_tag'].to_list(),
                        annotations['gene'].to_list()
                    )]

    new_featues_dict = dict(zip(feature_ids, features_new))


new_gbk = []
with open(args.gbk, 'r') as gb_file:
    for record in SeqIO.parse(gb_file, 'gb'):
        # Initialize dictionary to collect modified features
        # start_of_the_feature -> feature
        mod_features = dict()
        # Initialize list to collect
        for feature in record.features:
            # Get the start of the feature
            start = int(feature.location.start)
            # Add locus_tag and gene qualifiers if they are exist
            if not args.sort_only:
                feature_id = record.annotations.get('accessions')[0] + \
                     '_' + \
                     str(feature.location.start + 1) + \
                     '_' + \
                     str(feature.location.end)
                if new_featues_dict.get(feature_id, '') != '' and new_featues_dict.get(feature_id)[0] != '-':
                    feature.qualifiers['locus_tag'] = new_featues_dict.get(feature_id)[0]
                if new_featues_dict.get(feature_id, '') != '' and new_featues_dict.get(feature_id)[1] != '-':
                    feature.qualifiers['gene'] = new_featues_dict.get(feature_id)[1]

            # If start alrady exists increase it by 1 to not override previous feature
            while start in mod_features:
                start = start + 1
            mod_features[start] = feature

        # rebuild features
        record.features = [mod_features[start] for start in sorted(mod_features.keys())]
        new_gbk.append(record)

with open(args.out_gbk, 'w') as out_file:
    SeqIO.write(new_gbk, out_file, 'gb')

