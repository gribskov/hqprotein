"""=====================================================================================================================
codon_reference.py

Calculate reference codon usage tables based on file of CDS. Assumes that the input sequences are CDS (no intron or UTR)
and that the first base is the beginning of the first codon in the coding frame.

Writes three files codon_rf0.txt, codon_rf1.txt, codon_rf2.txt

TODO add option for counts or frequencies
TODO add option of codon preference

2026-02-03 gribskov
====================================================================================================================="""
import argparse
import datetime
import textwrap as _textwrap
from include.fasta import Fasta
from codon import Codon


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
    cds = Fasta(filename=opt.cds_file)

    # rf holds the codon information for the three reading frames
    rf = [Codon(), Codon(), Codon()]

    seq_n = 0
    while cds.next():
        # sum codon usage counts over all CDS sequences for all three reading frames
        seq_n += 1
        print(f'{seq_n}\t{cds.id}')
        for frame in range(3):
            # codon counts for all three frames, includes stop codons
            rf[frame].add_from_dna(cds.seq, frame)


    exit(0)
