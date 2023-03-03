#!/usr/bin/env python3
"""Change chromosome/scaffold names from ENSEMBL / GENCODE to UCSC."""
import argparse
import os


MAPPINGS_DIR = '../../../chr_mappings'
MAPPINGS_FILES = {
    'Human_hg38': [
        'hg38_UCSC2Gencode.txt',
    ],
    'mouse': [
        'mm10_UCSC2Gencode.txt',
    ],
    'fly': [
        'dm3_UCSC2Ensembl.txt',
    ],
    'yeast': [
        'R64-1-1_UCSC2ensembl.txt',
    ],
}


def parse_mapping_file(filename):
    """Parse mapping file."""
    mappings = {}
    with open(filename) as handle:
        for line in handle:
            line = line.strip().split()
            if len(line) < 2:
                # if line does not conatin names in both sources just skip it
                continue
            mappings[line[0]] = line[1]
    return mappings


def parse_mappings(species):
    """Parse file with chromosome mappings."""
    mappings = dict()
    if species not in MAPPINGS_FILES:
        raise ValueError('Species "{}" not supported.'.format(species))

    if MAPPINGS_FILES[species] == ['NO_MAPPING_FILE']:
        return mappings

    for basename in MAPPINGS_FILES[species]:
        filename = os.path.join(MAPPINGS_DIR, basename)
        mappings.update(parse_mapping_file(filename))
    return mappings


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="change ucsc to gencode/ensembl")
    parser.add_argument('--infile', help="Input bed file.")
    parser.add_argument('--outfile', help="Output bed file.")
    parser.add_argument('--species', help="Species.")
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    mappings = parse_mappings(args.species.strip("'").capitalize())

    with open(args.infile) as infile, open(args.outfile, 'wt') as outfile:
        for line in infile:
            line = line.strip().split('\t')
            if mappings and line[0] not in mappings:
                continue
            outfile.write('\t'.join([mappings.get(line[0], line[0])] + line[1:]) + '\n')


if __name__ == "__main__":
    main()
