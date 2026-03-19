"""=====================================================================================================================
pssm.py
class for storing and manipulating pssm (position specific scoring matrix)

2026-03-19 gribskov
====================================================================================================================="""
import datetime
import sys

import numpy as np


# import yaml

class PSSM:
    """=================================================================================================================

    ================================================================================================================="""

    def __init__(self, uid='pssm', rows='ACGT', columns=15, offset=0, fwidth=4, fprecision=0):
        """-------------------------------------------------------------------------------------------------------------
        uid: str             uid for this pssm; default='pssm'
        comment: str        any descriptive information
        creation: str       creation date: YYYY-MM-DD HH:MM:SS
        rows: str           alphabet for pssm, usually 'ACGT' or 'ACDEFGHIKLMNPQRSTVWY'; default='ACGT'
        columns: int        number of sequence positions in pssm; default = 15
        offset: int         offset for position to mark for occurrence
        n: int              number of observations; expect 1.0 for a frequency matrix
        matrix: np array    the actual data, rows are the letters of the alphabet, columns are sequence positions
        -------------------------------------------------------------------------------------------------------------"""
        self.uid = uid
        self.comment = ''
        self.creation = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
        self.rows = rows
        self.cols = columns
        self.offset = offset
        self.n = 0
        self.matrix = np.zeros((len(rows), columns), dtype=np.float64)

        # __str__ does not take arguments so define fieldwidth and precision within the object
        self.fieldwidth = fwidth
        self.precision = fprecision

    def __str__(self):
        """-------------------------------------------------------------------------------------------------------------
        formatted version of site. since __str__() does not accept parameters, the formatting parameters are  included
        in the object itself (self.fieldwidth, self.precision)

        :return: str    formatted string pssm metadata and values
        -------------------------------------------------------------------------------------------------------------"""
        divider = '|'
        fmt = f'{self.fieldwidth}.{self.precision}f'
        divider = f'{'|':>{self.fieldwidth}s}'
        dividerpos = self.offset

        self.set_n()
        outstr = f'uid: {self.uid}\ncreation: {self.creation}\nn: {self.n}\n'
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
                    outstr += f'{divider:{self.fieldwidth}}'

                outstr += f'{value:{fmt}}'
            outstr += '\n'

        return outstr

    def read(self, filename):
        """-------------------------------------------------------------------------------------------------------------
        read a pssm written to a file by __str__ as in print(pssm, file=fh). Format:
        donor:
            A   0.279 0.294 0.374 0.385 0.149     | 0.000 0.000 0.668 0.550 0.053 0.107 0.302 0.240 0.275 0.210
            C   0.294 0.275 0.252 0.172 0.195     | 0.000 0.019 0.042 0.183 0.027 0.156 0.275 0.248 0.275 0.302
            G   0.229 0.244 0.221 0.183 0.477     | 1.000 0.000 0.248 0.057 0.893 0.069 0.198 0.179 0.137 0.229
            T   0.198 0.187 0.153 0.260 0.179     | 0.000 0.981 0.042 0.210 0.027 0.668 0.225 0.332 0.313 0.260

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

# ======================================================================================================================
# Testing
# ======================================================================================================================
if __name__ == '__main__':
    # create a small pssm
    donor = PSSM(uid='donor', rows='ACGT', columns=5)
    donor.offset = 2
    donor.comment += 'comment1\n\ncomment2\n'
    print(donor)

    pssm = np.array([i for i in range(20)])
    pssm = np.resize(pssm, (4, 5))
    donor.matrix = pssm
    # donor.set_n()
    out = open('test.pssm', 'w')
    print(donor, file=out)
    out.close()

    pssm_from_file = PSSM()
    pssm_from_file.read('test.pssm')
    print(pssm_from_file)

    exit(0)
