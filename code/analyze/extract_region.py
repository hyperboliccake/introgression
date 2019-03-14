#!/usr/bin/env python3
import argparse
import os
import pickle
import gzip
import sys


def main():
    args = parse_args()
    args = validate_args(args)
    index = pickle.load(open(args['pickle'], 'rb'))
    locations = decode_regions(args['regions'], index, args['list_sort'])
    with gzip.open(args['filename'], 'rt') as reader:
        write_regions(reader, locations, args['suppress_header'])


def parse_args(args=None):
    '''
    Read in input arguments or the supplied list of strings
    Returns a dictionary of options
    '''
    parser = argparse.ArgumentParser(
        description='retrieve regions from indexed file.')

    parser.add_argument('regions',
                        nargs='+',
                        help='one or more region ids to retrieve')
    parser.add_argument('--filename',
                        required=True,
                        help='fa.gz file to look for regions')
    parser.add_argument('--list_sort',
                        action='store_true',
                        help='sort regions by the input order. Defualt sort by'
                        ' disk location')
    parser.add_argument('--suppress_header',
                        action='store_true',
                        help='suppress printing of header line in stdout')

    return vars(parser.parse_args(args))


def validate_args(args):
    '''
    Performs checks and conversions of input, raises ValueErrors if invalid
    '''
    if not os.path.exists(args['filename']):
        raise ValueError(f'{args["filename"]} not found')

    if args['filename'][-6:] != '.fa.gz':
        raise ValueError(f'{args["filename"]} expected to be .fa.gz')

    args['pickle'] = args['filename'][:-6] + '.pkl'
    if not os.path.exists(args['pickle']):
        raise ValueError(f'{args["pickle"]} not found with region file')

    parsed_regions = []
    for region in args['regions']:
        r = region
        if r[0] == 'r':
            r = r[1:]
        if not r.isdigit():
            raise ValueError(f'{region} could not be parsed')
        parsed_regions.append(int(r))
    args['regions'] = parsed_regions

    return args


def decode_regions(regions, index, retain_sort):
    '''
    Converts list of regions to file locations based on index dictionary
    Retain_sort controls if the output list order is determined by the
    region order or the disk location (i.e. values of index dict)
    '''

    try:
        result = [index[r] for r in regions]
    except KeyError as e:
        raise KeyError(f'r{e} not found in index')

    if retain_sort:
        return result
    else:
        return sorted(result)


def write_regions(reader, locations, suppress_header, num_lines=15):
    '''
    Writes the regions specified by index to stdout
    If print_header is false, ignore first line after location
    '''
    if suppress_header is True:
        num_lines -= 1

    for location in locations:
        reader.seek(location)
        if suppress_header is True:
            reader.readline()
        for i in range(num_lines):
            line = reader.readline()
            if line == '':
                print(f'{location} outside of file', file=sys.stderr)
                break
            else:
                print(line, end='')


if __name__ == '__main__':
    main()
