#!/usr/bin/env python
from __future__ import print_function
import contextlib
import logging
logging.basicConfig(format='%(levelname)s:%(asctime)s:%(message)s',
                    level=logging.DEBUG)


@contextlib.contextmanager
def smart_filereader(filename=None):
    """
    https://stackoverflow.com/a/17603000/1875814
    Use either filename or default to sys.stdout in a
    with... context manager
    """
    if filename and filename != '-':
        fh = open(filename)
    else:
        fh = sys.stdin

    try:
        yield fh
    finally:
        if fh is not sys.stdin:
            fh.close()


def parse_args():
    import argparse
    import textwrap
    the_description = '''\
    Scan SAM reads (either in a file passed via -r switch, or piped from
    e.g. samtools view) for sequences matching PCR primer pairs (from a
    file)'''
    parser = argparse.ArgumentParser(textwrap.dedent(the_description))
    parser.add_argument('primers', help='File containing primer sequences')
    parser.add_argument('-r', '--reads',
                        help='File containing reads (optional - can be piped)',
                        default='')
    return parser.parse_args()


def revcomp(s):
    complements = dict(
        A='T', C='G', G='C', T='A',
        M='K', R='Y', W='W', S='S',
        Y='R', K='M', V='B', H='D',
        D='H', B='V', N='N'
    )
    return ''.join(complements[char] for char in s.upper()[::-1])


def match_position(l):
    """
    Find first position in l that is not None, and return index
    """
    for pos, ele in enumerate(l):
        if ele is not None:
            return pos


def read_primers(filename):
    primers = []
    with open(filename) as fl:
        for line in fl:
            p1, p2 = line.rstrip().split()
            primers.append((p1, p2))
    return primers


def build_regex(primers_list):
    import re
    rgxs = [re.compile(r'(.*){}(.+){}(.*)'.format(primer[1], primer[0]))
            for primer in primers_list]
    rgxs.extend(
            re.compile(r'(.*){}(.+){}(.*)'.format(revcomp(primer[0]),
                       revcomp(primer[1])))
            for primer in primers_list)
    multirgx = re.compile('|'.join('(?:{0})'.format(x.pattern) for x in rgxs))
    return multirgx


if __name__ == '__main__':
    import os
    import sys
    args = parse_args()

    # Convert primers to regular expression
    if not os.path.exists(args.primers):
        logging.error('Primers file not found: {}'.format(args.primers))
        sys.exit()
    else:
        primers = read_primers(args.primers)
        logging.info('Read primers from {}'.format(args.primers))
        multirgx = build_regex(primers)
        logging.info('Built regular expression')

    # Work out if reads are coming from standard input or a file
    reads_filename = None
    if sys.stdin.isatty():  # not a pipe
        reads_filename = (args.reads)
        if not os.path.exists(reads_filename):
            logging.error('Reads file not found: {}'.format(reads_filename))
            sys.exit()

    print('READ\tCHROM\tPOS\tPRIMERID\tPREEXTENSION\tEXTENSION\t'
          'SEQUENCE\tLIGATION\tPOSTLIGATION\tUMI\tPRIMER.ORIENTATION\tREAD.ORIENTATION')
    try:
        with smart_filereader(reads_filename) as readsfile:
            lines_processed = 0
            for line in readsfile:
                lines_processed += 1
                fields = line.rstrip().split('\t')
                name = fields[0]
                flag = int(fields[1])
                chrom = fields[2]
                pos = fields[3]
                read_is_reverse = flag & 16 == 16
                seq = fields[9]
                search = multirgx.search(seq)
                if search:
                    rgx_pos = match_position(search.groups())
                    primer_pos = (rgx_pos // 3)  # Each regex has 3 capture groups -> original primer list is 3 times shorter
                    primer_is_reverse = primer_pos >= len(primers)
                    primer_pos = primer_pos % len(primers)  # Shift reverse-match position back to forward primer sequence
                    extension = primers[primer_pos][0]
                    ligation = primers[primer_pos][1]
                    pre_extension, captured, post_ligation = search.groups()[rgx_pos:rgx_pos+3]

                    if primer_is_reverse:
                        umi = post_ligation[:5]

                    else:
                        umi = pre_extension[-5:]

                    # Print the read name, index of matching primer (0-based), the matched sequences, and their orientation.
                    # Orientation note: output is displayed in same orientation as primer sequence,
                    # but if orientation is marked '-' the read sequence in the BAM is from the opposite strand
                    print(
                        ('{read}\t{chrom}\t{pos}\t{primerid}\t{preextension}\t{extension}\t'
                         '{sequence}\t{ligation}\t{postligation}\t{umi}\t'
                         '{primerorientation}\t{readorientation}')
                        .format(read=name, chrom=chrom, pos=pos, primerid=primer_pos,
                                preextension=pre_extension, extension=extension,
                                sequence=captured, ligation=ligation, umi=umi, postligation=post_ligation,
                                primerorientation=('-' if primer_is_reverse else '+'),
                                readorientation=('-' if read_is_reverse else '+')))
                if lines_processed % 100 == 0:
                    logging.debug('Processed {} lines'.format(lines_processed))
    except (BrokenPipeError, KeyboardInterrupt) as err:
        logging.debug("Broken pipe")
        sys.stderr.close()
        sys.exit()
