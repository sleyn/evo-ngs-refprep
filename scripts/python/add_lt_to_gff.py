#!/usr/bin/env python3

import argparse
import re
import pandas as pd
from Bio import SeqIO
import logging


class GffParser:
    def __init__(self):
        self.gff_repr = {'header': [], 'sequence-region': [], 'content': [], 'fasta': []}

    # Read gff file
    def read_gff(self, gff_file_path):
        logging.info('Read GFF file')
        # If file has as ##FASTA section switch behaviour 
        # to read all subsequent lines to the 'fasta' key
        fasta_switch = False
        
        with open(gff_file_path, 'r') as gff_file:
            for line in gff_file:
                # Remove new line symbol from the line
                line = line.strip()
                
                # Skip empty lines
                if len(line) == 0:
                    continue
                
                # Check if we got to the fasta section
                if line == '##FASTA':
                    fasta_switch = True
                
                # If we are in the fasta section append everything to the 'fasta' key
                if fasta_switch:
                    self.gff_repr['fasta'].append(line)
                    continue
                
                # Check if line is a header line
                if line[0] == '#':
                    # Skip line with sequence region information as it will be added based on the FNA
                    if line.split('\t')[0] == '##sequence-region':
                        continue
                    else:
                        # Add all other headers
                        self.gff_repr['header'].append(line)
                else:
                    items = line.split('\t')
                    # Separate attributes from other fields
                    attributes = re.findall(r'(\w+)=([^;]+)', items[8])
                    items = items[:-1]

                    items_dict = dict(
                        zip(
                            ['contig', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'],
                            items
                        )
                    )

                    # Add attributes
                    attr_names = [a[0] for a in attributes]
                    attr_values = [a[1] for a in attributes]
                    # Change "Name" to "product" for functional attribute
                    if 'Name' in attr_names:
                        attr_names[attr_names.index('Name')] = 'product'
                    items_dict['attributes'] = dict(zip(attr_names, attr_values))
                    self.gff_repr['content'].append(items_dict)

    # Add "sequence-region" header fields.
    def add_sequence_region(self, fna_file_path):
        logging.info('Add "sequence-region" header fields')
        with open(fna_file_path, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                self.gff_repr['sequence-region'].append(
                    f'##sequence-region\t{record.id}\t1\t{len(record)}'
                )

    @staticmethod
    # Sort coordinates if the coordinates are not in ascending order
    def _sort_coordinates(accession, start, end):
        if int(end) >= int(start):
            return accession + '_' + start + '_' + end
        else:
            return accession + '_' + end + '_' + start

    # Modify gff fields
    # Requires file with columns:
    # - accession
    # - start
    # - end
    # - strand
    # - refseq_locus_tag
    # - gene
    def add_attributes(self, feature_file):
        logging.info('Add "locus_tag" and "gene" fields')
        annotations = pd.read_csv(feature_file, sep='\t')
        annotations = annotations.fillna('-')

        # Sort coordinates and make ids for each feature
        feature_ids = [self._sort_coordinates(accession, start, end) for accession, start, end in
                       zip(
                           annotations['accession'].astype(str).to_list(),
                           annotations['start'].astype(str).to_list(),
                           annotations['end'].astype(str).to_list()
                       )]

        # In features save locus_tag and gene.
        features_new = [[locus_tag, gene] for locus_tag, gene in
                        zip(
                            annotations['refseq_locus_tag'].to_list(),
                            annotations['gene'].to_list()
                        )]

        new_featues_dict = dict(zip(feature_ids, features_new))

        for i in range(len(self.gff_repr['content'])):
            item = self.gff_repr['content'][i]
            item_id = item['contig'] + '_' + item['start'] + '_' + item['end']
            if item_id in new_featues_dict:
                if new_featues_dict[item_id][0] != '-':
                    self.gff_repr['content'][i]['attributes']['locus_tag'] = new_featues_dict[item_id][0]
                if new_featues_dict[item_id][1] != '-':
                    self.gff_repr['content'][i]['attributes']['gene'] = new_featues_dict[item_id][1]

    # Write GFF file:
    def write_gff(self, out_gff_path):
        logging.info('Write output GFF file')
        with open(out_gff_path, 'w') as out_gff:
            # Write headers
            for header_line in self.gff_repr['header']:
                out_gff.write(header_line + '\n')

            # Write sequence-region lines
            for sequence_region in self.gff_repr['sequence-region']:
                out_gff.write(sequence_region + '\n')

            # Write genome items
            for item in self.gff_repr['content']:
                out_gff.write(
                    '\t'.join([
                        item['contig'],
                        item['source'],
                        item['type'],
                        item['start'],
                        item['end'],
                        item['score'],
                        item['strand'],
                        item['phase']
                    ])
                )

                # Write attributes
                attr_field = []
                for attr in item['attributes']:
                    attr_field.append(f"{attr}={item['attributes'][attr]}")
                out_gff.write('\t' + ';'.join(attr_field) + '\n')
            
            # Write fasta section
            out_gff.write('\n'.join(self.gff_repr['fasta']))


def main():
    parser = argparse.ArgumentParser(description='Add locus tags and gene names to gff file')
    parser.add_argument('-g', '--gff', help='Input RAST GFF file')
    parser.add_argument('-n', '--fna', help='Input FNA file to build headers')
    parser.add_argument('-f', '--features', help='Features table. Required coulmns:'
                                                 'accession, start, end, strand, refseq_locus_tag, gene.')
    parser.add_argument('-o', '--out_gff', help='Output file.')
    args = parser.parse_args()
    
    gff = GffParser()
    gff.read_gff(args.gff)
    gff.add_sequence_region(args.fna)
    gff.add_attributes(args.features)
    gff.write_gff(args.out_gff)


if __name__ == "__main__":
    main()
