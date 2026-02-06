"""=====================================================================================================================
codon_reference.py

Calculate reference codon usage tables based on file of CDS. Assumes that the input sequences are CDS (no intron or UTR)
and that the first base is the beginning of the first codon in the coding frame, and the sequence ends at a stop codon.
Sequence length should be multiple of three.

Writes three files codon_rf0.txt, codon_rf1.txt, codon_rf2.txt

TODO add option for counts or frequencies
TODO add option of codon preference

2026-02-03 gribskov
====================================================================================================================="""
import argparse
import datetime
import textwrap as _textwrap

from codon import Codon
from include.sequence.fasta import Fasta


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """-----------------------------------------------------------------------------------------------------------------
    Custom formatter for argparse.  Less ugly breaking of lines.
    _textwrap.wrap(text, 90) sets wrap at 90 characters
    -----------------------------------------------------------------------------------------------------------------"""

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        self._max_help_position = 30
        return _textwrap.wrap(text, 90)


def getopts():
    """-----------------------------------------------------------------------------------------------------------------
    Get command line options with argparse

    :return: namespace with options
    -----------------------------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(description='Generate reference codon usage tables',
                                 formatter_class=CustomFormatter)
    # cl.add_argument('--atype', type=str, default='excel',
    #                 help='annotation file type (excel)')

    cl.add_argument('-c', '--count', action='store_true',
                    help='Write raw counts (default: False')
    cl.add_argument('-p', '--preference', action='store_true',
                    help='Calculate codon preference instead of codon usage (default: False')

    # files are strings, not opened automatically by argparse
    cl.add_argument('cds_file', type=str, help='reference genome coding sequences')
    cl.add_argument('codon', type=str, help='base name for codon table names')

    return cl.parse_args()


# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    current = datetime.datetime.now().strftime('%X %G.%b.%d ')
    print(f'Codon_reference - {current}')
    opt = getopts()

    print(f'\tcoding sequences: {opt.cds_file}')
    print(f'\toutput file prefix: {opt.codon}')
    if opt.count:
        print(f'\toutput type: counts')
    else:
        print(f'\toutput type: frequencies')

    print('\nreading coding sequences from {opt.cds_file}')

    cds = Fasta(filename=opt.cds_file)

    # rf holds the codon information for the three reading frames
    rf = [Codon(), Codon(), Codon()]
    seq_n = 0
    for seq in cds:
        # sum codon usage counts over all CDS sequences for all three reading frames
        seq_n += 1
        print(f'\t{seq_n}\t{cds.id}')
        for frame in range(3):
            # codon counts for all three frames, includes stop codons
            rf[frame].add_from_dna(cds.seq, frame)

    print(f'\nsequences read from {opt.cds_file}: {seq_n}')
    for i in range(3):
        print(f'\tcodons in frame {i}: {rf[i].n}')

    if opt.count:
        # print raw counts
        for i in range(3):
            outfile = f'{opt.codon}_rf{i}.count.txt'
            rf[i].to_file(outfile, fieldwidth=6, decimal=0)

    else:
        # print frequencies
        for i in range(3):
            outfile = f'{opt.codon}_rf{i}.frequency.txt'
            rf[i] /= rf[i].n
            rf[i].to_file(outfile, fieldwidth=6, decimal=4)

    exit(0)

    # codon preference not done yet
    familycount = []
    for frame in range(3):
        # calculate frequency of codon families in all three frames
        familycount.append(rf[frame].family())

    preference = []
    for frame in range(3):
        # codon preference for all three frames
        preference.append(rf[frame] / familycount[frame])

    exit(0)
