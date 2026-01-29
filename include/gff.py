import sys
import re


class Gff:
    """#############################################################################################
    Manipulate GFF3 and GTF annotation files.
    For GFF3 use
        anno = Gff(file=gtf)
        anno.attr_sep = '='
    For GTF (and probably GFF2)
        anno = Gff(file=gtf)
        anno.attr_sep = ' ' ; <space>, space is the default so normally you won't need to change

    attr_sep    separator between tag and value in attributes column; gtf:<space>, gff3:'='
    del_attr    remove 'attribute' key after parsing tag/value pairs (feature_parse())

    9 October 2019    Michael Gribskov
    #############################################################################################"""

    # predefined columns 0-7 in the file
    column = ['sequence', 'method', 'feature', 'begin', 'end', 'score', 'strand', 'frame',
              'attribute']

    def __init__(self, file=""):
        self.data = []
        self.gff_in = None
        self.attr_sep = ' '
        self.del_attr = None

        if file:
            fh = self.open(file)
            if fh:
                self.gff_in = fh

    def open(self, file):
        """-----------------------------------------------------------------------------------------
        open a file for reading

        :param self:
        :param file: str, path to a GFF file
        :return: fh or False if unsuccessful
        -----------------------------------------------------------------------------------------"""
        try:
            fh = open(file, 'r')
            return fh
        except IOError:
            sys.stderr.write("gff.open unable to open GFF file ({})".format(file))
            return False

    def read(self):
        """-----------------------------------------------------------------------------------------
        Read a line from gff_in and delegate to the proper parsing function

        :return:
        -----------------------------------------------------------------------------------------"""
        self.line = self.gff_in.readline().rstrip()

        if self.line:
            if self.line.startswith('#'):
                self.comment_parse()
            else:
                self.data.append(self.feature_parse())

            return True

        else:
            # EOF
            return False

    def read_all(self, limit):
        """-----------------------------------------------------------------------------------------
        read the entire file into memory

        :return: int, number of lines read
        -----------------------------------------------------------------------------------------"""
        line = 0
        while gff.read():
            line += 1

        return line

    def read_feature(self, feature_list):
        """-----------------------------------------------------------------------------------------
        Read the whole file and store only features in feature_list

        :param feature_list:
        :return: int, number of reatures read
        -----------------------------------------------------------------------------------------"""
        count = 0
        for self.line in self.gff_in:
            if self.line:
                if self.line.startswith('#'):
                    self.comment_parse()
                else:
                    parsed = self.feature_parse()
                    if parsed['feature'] in feature_list:
                        self.data.append(parsed)
                        count += 1
        return count

    def feature_parse(self):
        """-----------------------------------------------------------------------------------------
        parse a feature line. This works for gtf  if self.attr_sep is ';' should work for GFF3 if
        self.attr_sep is '='. the latter is not well tested

        :return:
        -----------------------------------------------------------------------------------------"""
        field = self.line.split(maxsplit=8)
        parsed = {}
        for i in range(len(field)):
            # extract the 9 defined columns
            if i in [3, 4]:
                # change begin, end, and score to int
                # should score be float?
                parsed[Gff.column[i]] = int(field[i])
            else:
                parsed[Gff.column[i]] = field[i]

        # split the attributs on ; and restore as a hash

        field = parsed['attribute'].rstrip().split(';')
        # attribute may end in; so last field may be blank
        if not field[-1]:
            field.pop()

        for f in field:
            (key, value) = f.strip().split(self.attr_sep, maxsplit=1)
            parsed[key] = value.strip('"')

        if self.del_attr:
            del parsed['attribute']

        return parsed

    def comment_parse(self):
        """-----------------------------------------------------------------------------------------
        intended for parsing comments; not implemented
        :return:
        -----------------------------------------------------------------------------------------"""
        pass
        return True

    def get_by_feature(self, key, start=0):
        """-----------------------------------------------------------------------------------------
        Generator for lines that match a feature tag

        :param key: str, string matching a feature in column 3
        :param start: int, line on which to start
        :return:
        -----------------------------------------------------------------------------------------"""
        data = self.data
        for n in range(start, len(data)):
            if data[n]['feature'] == key:
                yield n, data[n]

        raise StopIteration

    def get_by_value(self, column, key, start=0, stop=0):
        """-----------------------------------------------------------------------------------------
        A generator that returns rows where the specified column matches the specified value.
        Rows missing columns, e.g., those generated from attributes, are skipped

        :param column: str, predefined or attribute column
        :param key: str, column value to match
        :param start: int, beginning row
        :param stop: int, ending row + 1
        :return: row, content_hash
        -----------------------------------------------------------------------------------------"""
        data = self.data
        if stop == 0:
            stop = len(data)

        for n in range(start, stop):
            if column not in data[n]:
                continue

            if data[n][column] == key:
                yield n, data[n]

        return

    def get_by_sequence_range(self, targetfeature, targetseq, begin, end):
        """-----------------------------------------------------------------------------------------

        :param targetseq: list      list of valid sequences
        :param begin: int           begin of feature must be >= begin
        :param end: int             end of feature must be <= end
        :return:
        -----------------------------------------------------------------------------------------"""
        count = 0
        for self.line in self.gff_in:
            if self.line:
                if self.line.startswith('#'):
                    self.comment_parse()
                else:
                    parsed = self.feature_parse()
                    if parsed['feature'] not in targetfeature:
                        continue
                    if parsed['sequence'] not in targetseq:
                        continue
                    if parsed['begin'] < begin or parsed['begin'] > end:
                        continue
                    # if parsed['end'] > end:
                    #     continue

                    self.data.append(parsed)
                    count += 1
        return count

    def get_gene_id_groups_by_generic_selector(self, selector, groupby):
        """-----------------------------------------------------------------------------------------
        selector is a dict that specifies the columns of self.data that must match for a selection
        selectable fields:
            sequence    list of one or more sequences (field 0)
            feature     list of one of more features (field 2)
            begin       int, begin of feature must be >= begin
            end         int, begin of feature must be <= end

        :param selector: dict   see above
        :param groupby: str     indicated attribute must exist
        :return: int            number of regions selected
        -----------------------------------------------------------------------------------------"""
        count = 0
        n_group = 0
        self.data = {}
        for self.line in self.gff_in:
            if self.line:
                if self.line.startswith('#'):
                    continue
                else:
                    parsed = self.feature_parse()
                    selected = False
                    for s in selector:
                        if s == 'feature':
                            if not parsed['feature'] in selector[s]:
                                break
                        if s == 'sequence':
                            # print(f'sequence:{parsed['sequence']}\tselector:{selector[s]}')
                            if not parsed['sequence'] in selector[s]:
                                break
                        if s == 'begin':
                            if parsed['begin'] < selector[s]:
                                break
                        if s == 'end':
                            if parsed['begin'] > selector[s]:
                                break

                        # if you reach here you are selected
                        count += 1
                        group = parsed[groupby]
                        if group in self.data:
                            self.data[group].append(parsed)
                        else:
                            self.data[group] = [parsed]
                            n_group += 1

        return count

    def replace_by_column(self, column, find, replace):
        """-----------------------------------------------------------------------------------------
        Replace all of find by replace in column

        :param column: str, name of a predifined or attribute column
        :param find: str, string to replace
        :param replace: string to substitute for find
        :return: int, rows examined
        -----------------------------------------------------------------------------------------"""
        data = self.data
        n = 0
        for d in data:
            if column in d:
                d[column] = d[column].replace(find, replace)
                n += 1

        return n

    def replace_columns_re(self, column_list, search, replace=''):
        """-----------------------------------------------------------------------------------------
        replace strings in a list of columns using a regular expression
        -----------------------------------------------------------------------------------------"""
        data = self.data
        n = 0

        query = re.compile(search)
        for d in data:
            for column in column_list:
                if column in d:
                    d[column] = query.sub(replace, d[column])
                    n += 1

        return n

    def rename_key(self, old_key, new_key):
        """-----------------------------------------------------------------------------------------
        changes old_key to new_key in data hash.  Does nothing if old_key does not exist.  Useful
        for changing names of attribute columns when comparing GFF files from different sources.

        :param old_key: str, current key
        :param new_key: str, new _key
        :return: True if old_key exists
        -----------------------------------------------------------------------------------------"""
        data = self.data
        if old_key in data[0]:
            for d in data:
                d[new_key] = d.pop(old_key)

            return True

        return False

    def position_to_int(self):
        """-----------------------------------------------------------------------------------------
        convert begin and end positions from str to int
        :return: 
        -----------------------------------------------------------------------------------------"""
        for d in self.data:
            d['begin'] = int(d['begin'])
            d['end'] = int(d['end'])

        return True


# ==================================================================================================
# test
# ==================================================================================================
if __name__ == '__main__':

    # gff = Gff(file='stringtie.gff')
    gff = Gff(file='genome.gff')
    gff.attr_sep = '='

    # line = 0
    # while gff.read():
    #     line += 1

    # line = gff.read_all()
    line = transcripts = gff.read_feature(['mRNA'])
    sys.stdout.write('{} lines read\n'.format(line))

    # remove the string 'lcl|' in the sequence names
    # gff.replace_by_column('sequence', 'lcl|', '')
    query = r'(maker|augustus|masked|processed|gene)(-|_)'
    gff.replace_columns_re(['Parent', 'ID'], query, r'')

    # transcripts is a generator function
    # transcripts = gff.get_by_feature('transcript')
    transcripts = gff.get_by_feature('mRNA')
    (begin, transcript) = next(transcripts)

    for (end, transcript) in transcripts:
        # sys.stdout.write('{}\t{}\t{}\n'.format(transcript['gene_id'],
        #                                        transcript['transcript_id'],
        #                                        transcript['sequence']
        #                                        ))
        sys.stdout.write('{}\t{}\t{}\n'.format(transcript['Parent'],
                                               transcript['ID'],
                                               transcript['sequence']
                                               ))

        for (exon_n, exon) in gff.get_by_value('feature', 'exon', begin, end):
            print('\t{}\t{}\t{}\t{}'.
                  # format(exon['exon_number'], exon['begin'], exon['end'], exon['strand']))
                  format(exon['ID'], exon['begin'], exon['end'], exon['strand']))

        begin = end

    # the final transcript
    sys.stdout.write('{}\t{}\t{}\n'.format(transcript['Parent'],
                                           transcript['ID'],
                                           transcript['sequence']
                                           ))
    for (exon_n, exon) in gff.get_by_value('feature', 'exon', begin, len(gff.data)):
        print('\t{}\t{}\t{}\t{}'.
              # format(exon['exon_number'], exon['begin'], exon['end'], exon['strand']))
              format(exon['ID'], exon['begin'], exon['end'], exon['strand']))

    # print(flist)

    exit(0)
