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

    print('READ\tPRIMERPOS\tFIVEPRIME\tPRIMER1\t'
          'SEQUENCE\tPRIMER2\tTHREEPRIME\tORIENTATION')
    try:
        with smart_filereader(reads_filename) as readsfile:
            lines_processed = 0
            for line in readsfile:
                lines_processed += 1
                fields = line.rstrip().split('\t')
                name = fields[0]
                seq = fields[9]
                search = multirgx.search(seq)
                if search:
                    rgx_pos = match_position(search.groups())
                    primer_pos = (rgx_pos // 3)  # Each regex has 3 capture groups -> original primer list is 3 times shorter
                    is_reverse = primer_pos >= len(primers)
                    primer_pos = primer_pos % len(primers)  # Shift reverse-match position back to forward primer sequence
                    primer1 = primers[primer_pos][1]
                    primer2 = primers[primer_pos][0]
                    if is_reverse:
                        three_prime, captured, five_prime = [revcomp(match) for match in search.groups()[rgx_pos:rgx_pos+3]]
                    else:
                        five_prime, captured, three_prime = search.groups()[rgx_pos:rgx_pos+3]

                    # Print the read name, index of matching primer (0-based), the matched sequences, and their orientation.
                    # Orientation note: output is displayed in same orientation as primer sequence,
                    # but if orientation is marked '-' the read sequence in the BAM is from the opposite strand
                    print(
                        '{read}\t{primerpos}\t{fiveprime}\t{primer1}\t{sequence}\t{primer2}\t{threeprime}\t{orientation}'
                        .format(read=name, primerpos=primer_pos, fiveprime=five_prime, primer1=primer1,
                                sequence=captured, primer2=primer2, threeprime=three_prime,
                                orientation=('-' if is_reverse else '+')))
                if lines_processed % 100 == 0:
                    logging.debug('Processed {} lines'.format(lines_processed))
    except (BrokenPipeError, KeyboardInterrupt) as err:
        logging.debug("Broken pipe")
        sys.stderr.close()
        sys.exit()
