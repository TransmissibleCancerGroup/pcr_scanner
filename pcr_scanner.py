#!/usr/bin/env python
from __future__ import print_function
import contextlib
import logging
import fuzzysearch
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
    return primers, [(revcomp(p1), revcomp(p2)) for (p1, p2) in primers]


def fuzzy_search(string, primer_pairs):
    """
    Search `string` for both primer sequences
    """
    for (i, primer_pair) in enumerate(primer_pairs):
        search1 = fuzzysearch.find_near_matches(primer_pair[0], string, max_l_dist=2)
        if search1:
            search2 = fuzzysearch.find_near_matches(primer_pair[1], string, max_l_dist=2)
            if search2:
                return i, search1[0], search2[0]


if __name__ == '__main__':
    import os
    import sys
    args = parse_args()

    # Convert primers to regular expression
    if not os.path.exists(args.primers):
        logging.error('Primers file not found: {}'.format(args.primers))
        sys.exit()
    else:
        primers, rev_primers = read_primers(args.primers)
        logging.info('Read primers from {}'.format(args.primers))
        

    # Work out if reads are coming from standard input or a file
    reads_filename = None
    if sys.stdin.isatty():  # not a pipe
        reads_filename = (args.reads)
        if not os.path.exists(reads_filename):
            logging.error('Reads file not found: {}'.format(reads_filename))
            sys.exit()

    print('READ\tCHROM\tPOS\tPRIMERID\tPRELEFT\tLEFTPRIMER\t'
          'SEQUENCE\tRIGHTPRIMER\tPOSTRIGHT\tUMI\tPRIMER.ORIENTATION\tREAD.ORIENTATION\tSEQ')
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
                search = fuzzy_search(seq, primers)
                if search:
                    primer_is_reverse = False
                else:
                    search = fuzzy_search(seq, rev_primers)
                    if search:
                        primer_is_reverse = True
                if search:
                    primer_pos, match1, match2 = search
                    if (match1.start, match1.end) < (match2.start, match2.end):
                        l_s, l_e, r_s, r_e = match1.start, match1.end, match2.start, match2.end

                    else:
                        l_s, l_e, r_s, r_e = match2.start, match2.end, match1.start, match1.end
                    
                    pre_left = seq[:l_s]
                    left_primer = seq[l_s:l_e]
                    captured = seq[l_e:r_s]
                    right_primer = seq[r_s:r_e]
                    post_right = seq[r_e:]
                    
                    if primer_is_reverse:
                        umi = post_right[:5]

                    else:
                        umi = pre_left[-5:]

                    # Print the read name, index of matching primer (0-based), the matched sequences, and their orientation.
                    # Orientation note: output is displayed in same orientation as primer sequence,
                    # but if orientation is marked '-' the read sequence in the BAM is from the opposite strand
                    print(
                        ('{read}\t{chrom}\t{pos}\t{primerid}\t{preleft}\t{left}\t'
                         '{sequence}\t{right}\t{postright}\t{umi}\t'
                         '{primerorientation}\t{readorientation}\t{seq}')
                        .format(read=name, chrom=chrom, pos=pos, primerid=primer_pos,
                                preleft=pre_left, left=left_primer,
                                sequence=captured, right=right_primer, umi=umi, postright=post_right,
                                primerorientation=('-' if primer_is_reverse else '+'),
                                readorientation=('-' if read_is_reverse else '+'),
                                seq=seq))
                if lines_processed % 1000 == 0:
                    logging.debug('Processed {} lines'.format(lines_processed))
    except (BrokenPipeError, KeyboardInterrupt) as err:
        logging.debug("Broken pipe")
        sys.stderr.close()
        sys.exit()
