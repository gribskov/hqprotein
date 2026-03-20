"""=====================================================================================================================
pssm.py
class for storing and manipulating pssm (position specific scoring matrices)

2026-03-19 gribskov
====================================================================================================================="""
import datetime
import sys
import uuid
from copy import deepcopy

import numpy as np


# import yaml

class PSSM:
    """=================================================================================================================
    Position Specific Scoring matrices are a matrix where the rows correspond to characters in a sequence alphabet,
    columns correspond to positions in a sequence.

    Values can be integer counts, or frequencies (i.e., P(character|position). It is up to you to know what kind of
    information is stored. One way to discover is to look st PSSM.n which shows the maximum total counts in a column.
    This will normally be 1.0 for frequencies and a much larger number for counts.

    UUIDs are automatically created and get updated when the count/frequency content changes but this is still a little
    incomplete. Changes are tracked in the comments
    TODO add comments for adding counts

    Synopsis:
    # create a PSSM
    pssm = PSSM(title='donor', rows='ACGT', offset=pre)

    # add counts from a sequence
    pssm.add_counts_sequence(site)

    # convert counts to frequencies
    donor_frequencies = pssm.frequency()

    # calculate positional Shannon entropy in bits (range: 0 - 2 for DNA)
    I = donor_frequencies.information()
    ================================================================================================================="""
    # print formats are shared by all pssms via class variables
    fieldwidth=0
    precision=0

    def __init__(self, title='pssm', rows='ACGT', columns=15, offset=0):
        """-------------------------------------------------------------------------------------------------------------
        uid: str            uid for this pssm; first block of uuid4
        comment: str        any descriptive information, default='pssm'
        creation: str       creation date: YYYY-MM-DD HH:MM:SS -> see timestamp to globally change format
        rows: str           alphabet for pssm, usually 'ACGT' or 'ACDEFGHIKLMNPQRSTVWY'; default='ACGT'
        columns: int        number of sequence positions in pssm; default = 15
        offset: int         offset for position to mark for occurrence
        n: int              number of observations; expect 1.0 for a frequency matrix
        matrix: np array    the actual data, rows are the letters of the sequence alphabet, columns are sequence positions
                            can be either counts or frequencies
        a2i: dict           for conversion of sequence characters to integer index (used in add_counts_sequence())
        -------------------------------------------------------------------------------------------------------------"""
        self.uid = PSSM.uid()
        self.title = title
        self.comment = ''
        self.creation = PSSM.timestamp()
        self.rows = rows
        self.cols = columns
        self.offset = offset
        self.n = 0
        self.matrix = np.zeros((len(rows), columns), dtype=np.float64)
        self.a2i = {rows[i]: i for i in range(len(rows))}

        # __str__ does not take arguments so define fieldwidth and precision within the object. These are now global
        # self.fieldwidth = fwidth
        # self.precision = fprecision

    def __str__(self):
        """-------------------------------------------------------------------------------------------------------------
        formatted version of site. since __str__() does not accept parameters, the formatting parameters are  included
        in the object itself as global variables (PSSM.fieldwidth, PSSM.precision)
        TODO maybe divider should be global so it can be easily changed

        example formatted text (frequency):
        -----
        title: acceptor
        uid: e0c45f62
        creation: 2026-03-20 121651
        n: 1.0
        comment
            (2026-03-20 122020) copied from 64882c19 (acceptor)
            (2026-03-20 122020) converted to frequency
        pssm 4,15
            A  0.247  0.240  0.217  0.205  0.180  0.189  0.341  0.072  0.999  0.000      |  0.275  0.213  0.249  0.230  0.232
            C  0.326  0.290  0.306  0.306  0.290  0.332  0.250  0.678  0.000  0.001      |  0.222  0.306  0.333  0.313  0.311
            G  0.157  0.187  0.173  0.174  0.173  0.125  0.250  0.006  0.000  0.999      |  0.340  0.182  0.190  0.215  0.226
            T  0.270  0.283  0.304  0.314  0.357  0.354  0.159  0.244  0.000  0.000      |  0.163  0.299  0.228  0.243  0.230
        -----
        :return: str    formatted string pssm metadata and values
        -------------------------------------------------------------------------------------------------------------"""
        fmt = f'{PSSM.fieldwidth}.{PSSM.precision}f'
        divider = f'{'|':>{PSSM.fieldwidth}s}'
        dividerpos = self.offset

        self.set_n()
        outstr = f'title: {self.title}\nuid: {self.uid}\ncreation: {self.creation}\nn: {self.n}\n'
        all_comments = self.comment.split('\n')
        if all_comments:
            outstr += f'comment\n'
            for line in all_comments:
                if line:
                    outstr += f'\t{line}\n'

        outstr += f'pssm {len(self.rows)},{self.cols}\n'
        matrix = self.matrix
        for i, row in enumerate(matrix):
            outstr += f'\t{self.rows[i]}'
            # for pos in range(len(row)):
            for pos, value in enumerate(row):
                if pos == dividerpos:
                    outstr += f'{divider:{PSSM.fieldwidth}}'

                outstr += f'{value:{fmt}}'
            outstr += '\n'

        return outstr

    @staticmethod
    def timestamp():
        """-------------------------------------------------------------------------------------------------------------
        Formatted timestamp for logging: 2026-03-20 121651
        Warning: adding colons will probably break reading the formatted file which uses colons as a marker in parsing

        :return: str
        -------------------------------------------------------------------------------------------------------------"""
        ftime = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
        return ftime

    @staticmethod
    def uid():
        """-------------------------------------------------------------------------------------------------------------
        generate uuid using uuid4. only the first 8 character group is used. example: e0c45f6

        :return:str     uuid string
        -------------------------------------------------------------------------------------------------------------"""
        return str(uuid.uuid4())[:8]

    def read(self, filename):
        """-------------------------------------------------------------------------------------------------------------
        read a pssm written to a file by __str__ as in print(pssm, file=fh). Format:
        donor:
        example formatted text (frequency):
        -----
        title: acceptor
        uid: e0c45f62
        creation: 2026-03-20 121651
        n: 1.0
        comment
            (2026-03-20 122020) copied from 64882c19 (acceptor)
            (2026-03-20 122020) converted to frequency
        pssm 4,15
            A  0.247  0.240  0.217  0.205  0.180  0.189  0.341  0.072  0.999  0.000      |  0.275  0.213  0.249  0.230  0.232
            C  0.326  0.290  0.306  0.306  0.290  0.332  0.250  0.678  0.000  0.001      |  0.222  0.306  0.333  0.313  0.311
            G  0.157  0.187  0.173  0.174  0.173  0.125  0.250  0.006  0.000  0.999      |  0.340  0.182  0.190  0.215  0.226
            T  0.270  0.283  0.304  0.314  0.357  0.354  0.159  0.244  0.000  0.000      |  0.163  0.299  0.228  0.243  0.230
        -----
        :param filename: str    path to a readable file
        :return: float          number of donor sequences, 1.0 indicates frequencies
        -------------------------------------------------------------------------------------------------------------"""
        try:
            fh = open(filename, 'r')
        except (IOError, OSError):
            sys.stderr.write(f'SpliceSite::read - unable to open ({filename}) for reading')
            exit(1)

        # pssm = self
        for line in fh:
            line = line.rstrip()
            if not line: continue

            if line.startswith('comment'):
                # comment lines begin with tab, the first line without a tab is the next section (should be pssm)
                for line in fh:
                    if line.startswith('\t'):
                        self.comment += f'{line[1:]}\n'
                    else:
                        break

            if line.startswith('pssm'):
                tag, shape = line.split(' ', maxsplit=1)
                rows, cols = shape.strip().split(',')
                rows = int(rows.strip())
                cols = int(cols.strip())
                # matrix = np.array((rows, cols), dtype=np.float64)
                matrix = np.zeros((rows, cols), dtype=np.float64)
                alphabet = ''
                row = 0
                while row < rows:
                    line = fh.readline()
                    field = line.strip().split()
                    alphabet += field[0]
                    col = 0
                    for value in field[1:]:
                        if value == '|':
                            self.offset = col
                        else:
                            matrix[row, col] = float(value)
                            col += 1
                    row += 1
                self.matrix = matrix

            elif line.find(':') >= 0:
                # line contains : read it as a value
                tag, value = line.split(':')
                tag = tag.strip()
                value = value.strip()
                setattr(self, tag, value)
            else:
                # something goes wrong
                sys.stderr.write(f'PSSM::read - unparseable line: {line}')

        return self.uid

    def set_n(self):
        """-------------------------------------------------------------------------------------------------------------
        sets self.n to be the largest sum across the columns

        :return: self.n
        -------------------------------------------------------------------------------------------------------------"""
        self.n = max(np.sum(self.matrix, axis=0))
        return self.n

    def copy(self):
        """-------------------------------------------------------------------------------------------------------------

        :return:
        -------------------------------------------------------------------------------------------------------------"""
        frequency = deepcopy(self)
        frequency.uid = PSSM.uid()
        frequency.comment += f'({PSSM.timestamp()}) copied from {self.uid} ({self.title})\n'
        return frequency

    def add_counts_sequence(self, sequence):
        """-------------------------------------------------------------------------------------------------------------
        Add sequence counts to PSSM. be careful that your pssm has not been converted to frequency. increments pssm.n,
        the count of bases at each position

        :param sequence: str        DNA or protein sequence
        :return: int                count of bases at each position
        -------------------------------------------------------------------------------------------------------------"""
        sarr = np.array(list(sequence))
        count = np.zeros((4, len(sequence)), dtype=int)
        a2i = self.a2i
        for base in a2i:
            mask = sarr == base
            count[a2i[base], mask] += 1

        self.matrix += count
        self.n += 1
        return self.n

    def frequency(self):
        """------------------------------------------------------------------------------------------------------------
        Return a new PSSM object with counts converted to frequencies.

        :return: PSSM     positional probability of base | splice_site
        -------------------------------------------------------------------------------------------------------------"""
        frequency = self.copy()
        frequency.comment += f'({PSSM.timestamp()}) converted to frequency\n'

        total = np.sum(self.matrix, axis=0)
        for base in range(len(self.rows)):
            frequency.matrix[base, :] /= total

        frequency.set_n()

        return frequency

    def information(self):
        """-------------------------------------------------------------------------------------------------------------
        Calculate positional Shannon information for self.matrix

        :return: list     positional information for pssm
        -------------------------------------------------------------------------------------------------------------"""
        # h = np.zeros((len(self.rows),self.cols), dtype=float)
        base = np.log2(len(self.rows))
        with np.errstate(divide='ignore'):
            # ignore 'divide by zero encountered in log2' errors
            h = self.matrix * np.nan_to_num(np.log2(self.matrix), nan=0.0)
        # h = np.nan_to_num(h, nan=0.0)
        return base + np.sum(h, axis=0)

    def sharpen(self, exponent, renormalize=True):
        """-------------------------------------------------------------------------------------------------------------
        Sharpen/flatten values by taking values to an exponent. Exponent > 1 sharpens, exponent < 1 flattens
        Sharpening updates the uid

        :param exponent: float      each value is replaced by value ** exponent
        :param renormalize: bool    whether to normalize the values by dividing by sum
        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        uid = self.uid
        self.uid = PSSM.uid()
        self.comment += f'({PSSM.timestamp()}) sharpened with exponent={exponent}, source= {uid} ({self.title})\n'
        matrix = self.matrix
        matrix **= exponent

        if renormalize:
            matrix /= np.sum(self.matrix, axis=0)

        return True


# ======================================================================================================================
# Testing
# ======================================================================================================================
if __name__ == '__main__':
    # create a small pssm
    donor = PSSM(title='donor', rows='ACGT', columns=5)
    donor.offset = 2
    donor.comment += 'comment1\n\ncomment2\n'
    print(donor)

    # fill with sequential values and print to file
    pssm = np.array([i for i in range(20)])
    pssm = np.resize(pssm, (4, 5))
    donor.matrix = pssm
    out = open('test.pssm', 'w')
    print(donor, file=out)
    out.close()

    # reread file
    pssm_from_file = PSSM()
    pssm_from_file.read('test.pssm')
    print(pssm_from_file)

    # frequency matrix
    frequency = pssm_from_file.frequency()
    frequency.fieldwidth = 7
    frequency.precision = 3
    print(frequency)

    # positional information content
    I = frequency.information()
    print(I)
    fmt = f'{frequency.fieldwidth}.{frequency.precision}f'
    fI = '\t '
    for i in I:
        fI += f'{i:{fmt}}'
    print(fI)

    # frequency sharpening
    frequency.sharpen(1.3, renormalize=True)
    print(frequency)
    I = frequency.information()
    print(I)
    fmt = f'{frequency.fieldwidth}.{frequency.precision}f'
    fI = '\t '
    for i in I:
        fI += f'{i:{fmt}}'
    print(fI)

    exit(0)
