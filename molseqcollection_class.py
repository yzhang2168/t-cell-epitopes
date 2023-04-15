import re
from molseq_class import MolSeq
from fasta_parser import parse_fasta_file
from substring_search_rabin_karp import rabin_karp_substring_search


class MolSeqCollection(object):
    ID_RE = re.compile(">\s*(.+)")
    GAP_RE = re.compile("[-\s]+")
    WS_RE = re.compile("[\s]+")

    """
    A collection of MolSeq objects, either aligned or unaligned, stored in a list
    """
    def __init__(self, is_aligned=False):
        self.__molseq_list = []
        self.__label_list = set()
        self.__is_aligned = is_aligned

    def is_aligned(self):
        return self.__is_aligned

    def __str__(self):
        return self.to_fasta()

    def __len__(self):
        return self.get_length()

    @staticmethod
    def parse(obj, format="fasta", is_aligned=False, remove_gaps=False, allow_duplicate_labels=False):
        """
        static elements, bound to the class. defined outside of __init__. Use MSA.name_list or object.name_list.
        :param obj:
        :param format: only fasta is supported right now
        :param is_aligned:
        :param remove_gaps:
        :param allow_duplicate_labels:
        :return:
        """
        if format == "fasta":
            MolSeq_list, count_seqs_with_gap = parse_fasta_file(obj, remove_gaps)
            new = MolSeqCollection(is_aligned)
            for seq in MolSeq_list:
                new.add_molseq(seq, allow_duplicate_labels)
            if remove_gaps:
                return new, count_seqs_with_gap
            else:
                return new
        #elif format is "phylip":
        else:
            raise ValueError("unknown format: " + format)

    def add_molseq(self, molseq, allow_duplicate_labels=False):
        # NOTE: By maintaining a set (as instance variable) of all
        # labels encountered, you _could_ enforce unique labels.
        if self.__is_aligned:
            if self.__molseq_list != [] and len(molseq) != len(self.__molseq_list[0]):
                raise ValueError('The input sequence is not of the same length as the alignment')
        if ((not allow_duplicate_labels) and (molseq.get_label() in self.__label_list)):
            #print(self.__label_list)
            print('The label already exists:', molseq.get_label(), "skipping it")
            #raise ValueError('The label already exists:', molseq.get_label())
        else:
            self.__molseq_list.append(molseq)
            self.__label_list.add(molseq.get_label())

    def add_seq_from_strings(self, name, seq_str):
        self.add_molseq(MolSeq(name, seq_str))

    def add_molseq_list(self, molseq_list):
        for s in molseq_list:
            self.add_molseq(molseq=s, allow_duplicate_labels=False)

    def get_labels(self):
        return self.__label_list

    def get_molseq_list(self):
        return self.__molseq_list

    def get_seqs_by_labels(self, label_list):
        print("getting seqs by label list")
        new = MolSeqCollection()
        count = 0
        for label in label_list:
            found = False
            for seq in self.__molseq_list:
                if str(label) == seq.get_label():
                    new.add_molseq(seq)
                    count += 1
                    found = True
            if not found:
                    print(label + " was not found")
        print('# sequences wrote ', count)
        return new

    def get_seqs_by_label_substring(self, delimiter, substring_index, label_list, substring_delimiter=''):
        """
        select seqs based on label substring
        :param delimiter: # '|' in >gb:KY417142|Organism:Bat SARS-like coronavirus|Strain Name:As6526|Segment:null|Host:Bat
        :param substring_index: 0-indexed
        :param substring_delimiter: # ':' in ncbiId:AGA37377.1
        :param label_list:
        :return: a new MolSeqCollection object
        """
        new = MolSeqCollection()
        count = 0
        for molseq in self.__molseq_list:
            if substring_index in range(len(molseq.get_label().split(delimiter))):
                substring = molseq.get_label().split(delimiter)[substring_index]
                if substring_delimiter:
                    # ncbiId:AGA37377.1: get AGA37377.1
                    substring = substring.split(substring_delimiter)[1]
                if substring in label_list:
                    new.add_molseq(molseq)
                    count += 1
            else:
                raise ValueError('substring index is out of range')
        print('# sequences wrote ', count)
        return new

    def get_seqs_by_label_find_substring(self, substring):
        """
        select seqs based on a match of input string in any part of the label
        :param label_list:
        :return: a new MolSeqCollection object
        """
        new = MolSeqCollection(is_aligned=self.is_aligned())
        count = 0
        for molseq in self.__molseq_list:
            label = molseq.get_label()
            if label.lower().find(substring.lower()) != -1:
                new.add_molseq(molseq)
                count += 1
        print('# sequences wrote ', count)
        return new

    def get_seqs_by_label_substr_list(protein_fasta, label_list, outfile):
        seqs = MolSeqCollection.parse(protein_fasta)
        new = MolSeqCollection()
        count = 0
        for seq_label in label_list:
            found = False
            for molseq in seqs.get_molseq_list():
                label = molseq.get_label()
                if label.lower().find(seq_label.lower()) != -1:
                    new.add_molseq(molseq)
                    count += 1
                    found = True
            if not found:
                print(seq_label + " was not found")
        print('# sequences wrote ', count)
        new.write_to_file(outfile)
        return new

    def get_seqs_by_label_find_keyword_in_field(self, delimiter, field_index, keyword_list):
        """
        select seqs based on label substring
        :param delimiter: # '|' in >gb:KY417142|Organism:Bat SARS-like coronavirus|Strain Name:As6526|Segment:null|Host:Bat
        :param substring_index: 0-indexed
        :param substring_delimiter: # ':' in ncbiId:AGA37377.1
        :param label_list:
        :return: a new MolSeqCollection object
        """
        new = MolSeqCollection(is_aligned=self.is_aligned())
        count = 0
        for molseq in self.__molseq_list:
            if field_index in range(len(molseq.get_label().split(delimiter))):
                field = molseq.get_label().split(delimiter)[field_index]
                for keyword in keyword_list:
                    if field.lower().find(keyword.lower()) != -1:
                        new.add_molseq(molseq, allow_duplicate_labels=False)
                        count += 1
            else:
                print('substring index is out of range:', molseq.get_label())
                #raise ValueError('substring index is out of range')
        print('# sequences wrote ', count)
        return new

    def get_seqs_excluding_labels(self, label_list):
        new = MolSeqCollection(is_aligned=self.is_aligned())
        count = 0
        for seq in self.__molseq_list:
            # ncbiId:AGA37377.1: get AGA37377.1
            #if seq.get_label().split(':')[1] not in label_list:
            if seq.get_label() not in label_list:
                new.add_molseq(seq)
                count += 1
        print('# sequences wrote ', count)
        return new

    def get_seqs_excluding_label_substring(self, delimiter, substring_index, label_list, substring_delimiter=''):
        """
        select seqs based on label substring
        :param delimiter: # '|' in >gb:KY417142|Organism:Bat SARS-like coronavirus|Strain Name:As6526|Segment:null|Host:Bat
        :param substring_index: 0-indexed
        :param substring_delimiter: # ':' in ncbiId:AGA37377.1
        :param label_list:
        :return: a new MolSeqCollection object
        """
        new = MolSeqCollection(is_aligned=self.is_aligned())
        count = 0
        for molseq in self.__molseq_list:
            if substring_index in range(len(molseq.get_label().split(delimiter))):
                substring = molseq.get_label().split(delimiter)[substring_index]
                if substring_delimiter:
                    # ncbiId:AGA37377.1: get AGA37377.1
                    substring = substring.split(substring_delimiter)[1]
                if substring not in label_list:
                    new.add_molseq(molseq)
                    count += 1
            else:
                raise ValueError('substring index is out of range')
        print('# sequences wrote ', count)
        return new

    def set_label_w_map(self, map_dict):
        new = MolSeqCollection(is_aligned=self.is_aligned())
        for molseq in self.__molseq_list:
            if molseq.get_label() in map_dict:
                molseq.set_label(map_dict[molseq.get_label()])
                new.add_molseq(molseq)
            else:
                print(molseq.get_label(), 'not found in mapping dictionary')
                new.add_molseq(molseq)
                #raise ValueError(molseq.get_label(), 'not found in mapping dictionary')
        return new

    def set_label_w_substring(self, delimiter, substring_index):
        new = MolSeqCollection(is_aligned=self.is_aligned())
        for molseq in self.__molseq_list:
            if substring_index in range(len(molseq.get_label().split(delimiter))):
                substring = molseq.get_label().split(delimiter)[substring_index]
                molseq.set_label(substring)
                new.add_molseq(molseq)
            else:
                raise ValueError('substring index is out of range')
        return new

    def get_seqs_by_indices(self, index_list):
        new = MolSeqCollection(is_aligned=self.is_aligned())
        for index in index_list:
            index -= 1
            if index in range(self.get_number_sequences()):
                new.add_molseq(self.__molseq_list[index])
            else:
                raise ValueError('invalid index')
        return new

    def get_number_sequences(self):
        return len(self.__molseq_list)

    def get_length(self):
        # NOTE: If collection is empty, 0 should be returned as length, not an exception.
        # NOTE: AttributeError is not really a good choice if unaligned, maybe use GeneralException?
        if self.__molseq_list: # empty list evaluates to False
            if self.is_aligned():
                return len(self.__molseq_list[0])
            else:
                raise AttributeError('unaligned sequences')
        else: # empty list evaluates to False
            return 0

    def get_region(self, start, end):
        if self.is_aligned():
            if (start > 0) and (end <= self.get_length()) and (start <= end):
                new_alignment = MolSeqCollection(is_aligned=True)
                for seq in self.__molseq_list:
                    new_seq = MolSeq(seq.get_label(), seq.get_seq()[start-1:end])
                    new_alignment.add_molseq(new_seq)
                return new_alignment
            else:
                raise ValueError('invalid region')
        else:
            raise AttributeError('unaligned sequences')

    def get_seqs_by_length(self, min_length):
        count = 0
        new = MolSeqCollection(is_aligned=self.is_aligned())
        for seq in self.__molseq_list:
            if seq.get_length() >= min_length:
                new.add_molseq(seq)
                count += 1
        print("after filtering", count)
        return new

    def get_seqs_by_pattern(self, pattern):
        """

        :param pattern: substring in seq string
        :return: a new MolSeqCollection object
        """
        new = MolSeqCollection(is_aligned=self.is_aligned())
        ambig_list = []
        for seq in self.__molseq_list:
            substring_indices = rabin_karp_substring_search(pattern, seq.get_seq())
            if substring_indices:
                ambig_list.append(seq.get_label())
                new.add_molseq(seq)
        print(ambig_list)
        return new

    def get_seqs_excluding_pattern(self, pattern):
        """

        :param pattern: substring in seq string
        :return: a new MolSeqCollection object
        """
        new = MolSeqCollection(is_aligned=self.is_aligned())
        ambig_list = []
        for seq in self.__molseq_list:
            substring_indices = rabin_karp_substring_search(pattern, seq.get_seq())
            if substring_indices:
                ambig_list.append(seq.get_label())
            else:
                new.add_molseq(seq)
        print(ambig_list)
        return new

    def get_gap_only_columns(self):
        """
        :return: a list of coordinates for gap_only_columns, 1-indexed
        """
        if self.is_aligned():
            gap_only_cols = []
            length = self.get_length()
            for i in range(length):
                gap_only = True
                for seq in self.__molseq_list:
                    if seq.get_seq()[i] != '-':
                        gap_only = False
                        break
                if gap_only:
                    gap_only_cols.append(i + 1)
        else:
            raise AttributeError('unaligned sequences')
        return gap_only_cols

    def delete_gap_only_columns(self):
        """
        :return: a MolSeqCollection object with gap_only_columns deleted
        """
        if self.is_aligned():
            if self.get_length() > 0:
                # gap_only_cols: 1-based indexing
                gap_only_cols = self.get_gap_only_columns()
                #print(gap_only_cols)

                new_alignment = MolSeqCollection(is_aligned=True)
                for seq in self.__molseq_list:
                    seq_char_list = list(seq.get_seq())
                    for i in range(len(gap_only_cols), 0, -1):
                        # gap_only_cols: -1 to convert to 0-based indexing
                        seq_char_list.pop(gap_only_cols[i - 1] - 1)
                    new_seq = MolSeq(seq.get_label(), "".join(seq_char_list))
                    new_alignment.add_molseq(new_seq)
                return new_alignment
            else:
                raise ValueError('invalid region')
        else:
            raise AttributeError('unaligned sequences')

    def get_aln_coords_for_substring(self, seq_label, query):
        """

        :param seq_label:
        :param seq_substring:
        :return: 1-based index (start, end) in alignment with gaps if there is any,
        (-1, -1) if seq_label not found or seq_substring not found

        print(get_aln_coords_for_substring('---A-----', ' -     ')) # (-1, -1)
        print(get_aln_coords_for_substring('---A-A---', ' a -   ')) # (4, 4)
        print(get_aln_coords_for_substring('ABCDE----', ' ABCDE ')) # (1, 5)
        print(get_aln_coords_for_substring('--A--BCDE', ' ABCDE-')) # (3, 9)
        print(get_aln_coords_for_substring('--A-B-CDE', ' ABCDE-')) # (3, 9)
        print(get_aln_coords_for_substring('--ABCDE--', ' ABCDE-')) # (3, 7)
        print(get_aln_coords_for_substring('--A--BCDE', ' ABCDEE')) # (-1, -1)
        print(get_aln_coords_for_substring('--A--BCDE', ' ABXDE ')) # (-1, -1)
        print(get_aln_coords_for_substring('-ABAAABCD', ' AABCDE')) # (-1, -1)
        print(get_aln_coords_for_substring('--ABABCDE', ' ABCDE-')) # (5, 9)
        print(get_aln_coords_for_substring('--AAAACA-', ' AAAC -')) # (4, 7)
        """

        # case 0. query is None or empty
        if query is None:
            return -1, -1

        # case 1. query is None or empty
        query = re.sub(MolSeqCollection.GAP_RE, '', query)
        if len(query) == 0:
            return -1, -1

        # case 2. sequence is not in aln
        molseqs = self.get_seqs_by_label_find_substring(seq_label)
        if molseqs.get_number_sequences() == 0:
            print("target sequence not found")
            return -1, -1

        # case 3. aln_str is empty or query length > aln_str length
        aln_str = molseqs.get_molseq_list()[0].get_seq()
        if len(aln_str) == 0 or len(query) > len(aln_str):
            return -1, -1

        # case 4. passed input validation
        # to lowercase for char comparison
        query = query.lower()
        aln_str = aln_str.lower()

        i = 0 # query
        j = 0 # aln_str
        aln_start = -1
        aln_end = -1
        '''
        while j <= len(aln_str) - len(query) + i and i <= len(query) - 1:
            # find the 1st char
            if i == 0:
                while j <= len(aln_str) - len(query) + i and aln_str[j] != query[i]:
                    j += 1

                # j is out of search space or 1st char found
                if j > len(aln_str) - len(query) + i:
                    break
                else:
                    # aln_str[j] == query[i]
                    aln_start = j

                    # case 4.1 substring len = 1
                    if len(query) == 1:
                        aln_end = j
                        break
                    else:
                        i += 1
                        j += 1

            # case 4.2 substring len >= 2
            # i = 1, j = aln_start
            # skip '-' in aln_str
            while j <= len(aln_str) - len(query) + i and aln_str[j] == '-':
                j += 1

            # j is out of search space
            if j > len(aln_str) - len(query) + i:
                aln_start = -1
                break
            elif i <= len(query) - 1 and aln_str[j] == query[i]:
                # j is in search space and aln_str[j] != '-'
                if i == len(query) - 1:
                    aln_end = j
                    break
                else:
                    i += 1
                    j += 1
            else:
                i = 0
                # to handle overlapping substrings, restart at next position after previous start
                j = aln_start + 1
                aln_start = -1
        '''
        while i <= len(query) - 1 and j <= len(aln_str) - len(query) + i:
            # look for 1st char
            if i == 0:
                while j <= len(aln_str) - len(query) + i and query[i] != aln_str[j]:
                    j += 1
            else:
                # 2nd char onward
                while j <= len(aln_str) - len(query) + i and aln_str[j] == '-':
                    j += 1

            # j is out of search space
            if j > len(aln_str) - len(query) + i:
                break
            elif i == 0 or aln_str[j] == query[i]:
                if i == 0:
                    # 1st char found
                    aln_start = j
                if i == len(query) - 1:
                    # done with query
                    aln_end = j
                    break
                else:
                    i += 1
                    j += 1
            else:
                # reset
                i = 0
                j = aln_start + 1
                aln_start = -1

        # end of outer while
        if i < len(query) - 1:
            aln_start = -1
            return -1, -1
        else:
            return aln_start + 1, aln_end + 1

    def get_non_gap_coordinates_for_seq_substring(self, seq_label, query):
        # case 0. query is None or empty
        if query is None:
            return -1, -1

        # case 1. query is None or empty
        query = re.sub(MolSeqCollection.GAP_RE, '', query)
        if len(query) == 0:
            return -1, -1

        # case 2. sequence is not in aln
        molseqs = self.get_seqs_by_label_find_substring(seq_label)
        if molseqs.get_number_sequences() == 0:
            print("target sequence not found")
            return -1, -1

        # case 3. aln_str is empty or query length > aln_str length
        aln_str = molseqs.get_molseq_list()[0].get_seq()
        if len(aln_str) == 0 or len(query) > len(aln_str):
            return -1, -1

        # case 4. passed input validation
        # to lowercase for char comparison
        query = query.lower()
        aln_str = aln_str.lower()

        # get start, end in target seq itself without gaps
        seq_string_no_gaps = MolSeq(seq_label, aln_str.replace('-', ''))
        no_gap_start, no_gap_end = seq_string_no_gaps.get_coordinates(query)
        return no_gap_start, no_gap_end

    def to_fasta(self, chars_per_line=60):
        fasta = ''
        for molseq in self.__molseq_list:
            fasta += molseq.to_fasta(chars_per_line) + '\n'
        return fasta

    def write_to_file(self, outfile, chars_per_line=60):
        with open(outfile, 'w') as f:
            f.write(self.to_fasta(chars_per_line))

    def write_lowercase_md5_to_file(self, outfile):
        with open(outfile, 'w') as f, open(outfile + '.labels.tsv', 'w') as f2:
            f.write('Label\tLowercase MD5\n')
            seqs = 0
            digest_to_labels = {}
            max_count = 0
            max_count_digest = ''
            for molseq in self.__molseq_list:
                seqs += 1
                digest = molseq.get_lowercase_md5()
                f.write(molseq.get_label() + '\t' + digest + '\n')

                if digest in digest_to_labels:
                    digest_to_labels[digest].append(molseq.get_label())
                else:
                    digest_to_labels[digest] = [molseq.get_label()]

                if len(digest_to_labels[digest]) > max_count:
                    max_count = len(digest_to_labels[digest])
                    max_count_digest = digest

            f2.write('MD5\tCount\tLabels\n')
            for key, value in digest_to_labels.items():
                f2.write(key + '\t' + str(len(value)) + '\t' + ','.join(value) + '\n')

        print('# sequences:', seqs)
        print('# MD5 values:', len(digest_to_labels))
        print('MD5 max count:', max_count)
        print('MD5 with max count:', max_count_digest)

    def write_hash_to_file(self, hash_function, outfile):
        #if hash_function.lower() != 'md5' and hash_function.lower() != 'sha256' and hash_function.lower() != 'sha512':
        #    raise ValueError('invalid hash function name')

        with open(outfile, 'w') as f, open(outfile + '.labels.tsv', 'w') as f2:
            f.write('Label\t' + hash_function + '\n')
            seqs = 0
            digest_to_labels = {}
            max_count = 0
            max_count_digest = ''
            for molseq in self.__molseq_list:
                seqs += 1
                digest = ''
                if hash_function.lower() == 'md5':
                    digest = molseq.get_md5()
                elif hash_function.lower() == 'sha256':
                    digest = molseq.get_sha256()
                    #digest = molseq.get_lowercase_sha256()
                elif hash_function.lower() == 'sha512':
                    digest = molseq.get_sha512()

                f.write(molseq.get_label() + '\t' + digest + '\n')

                if digest in digest_to_labels:
                    digest_to_labels[digest].append(molseq.get_label())
                else:
                    digest_to_labels[digest] = [molseq.get_label()]

                if len(digest_to_labels[digest]) > max_count:
                    max_count = len(digest_to_labels[digest])
                    max_count_digest = digest

            f2.write(hash_function + '\tCount\tLabels\n')
            for key, value in digest_to_labels.items():
                f2.write(key + '\t' + str(len(value)) + '\t' + ','.join(value) + '\n')

        print('# sequences:', seqs)
        print('# hash values:', len(digest_to_labels))
        print('Hash with max count:', max_count_digest)
        print('Hash max count:', max_count)
        #print(self.test_collison(digest_to_labels[max_count_digest]))

    def test_collison(self, labels):
        molseqs = self.get_seqs_by_labels(labels)
        seq_set = set()
        for seq in molseqs.__molseq_list:
            seq_set.add(seq.get_seq())
        print('# unique seqs with the same hash value', len(seq_set))
        if len(seq_set) > 1:
            print(seq_set)
        return len(seq_set) > 1

    def count_unique_seqs(self):
        seq_to_labels = {}
        seqs = 0
        max_count = 0

        for molseq in self.__molseq_list:
            seqs += 1
            seq_str = molseq.get_seq()
            if seq_str in seq_to_labels:
                seq_to_labels[seq_str].append(molseq.get_label())
            else:
                seq_to_labels[seq_str] = [molseq.get_label()]

            if len(seq_to_labels[seq_str]) > max_count:
                max_count = len(seq_to_labels[seq_str])

        with open('test_count_unique_seqs.tsv', 'w') as f2:
            f2.write('Sequence\tCount\tLabels\n')
            for key, value in seq_to_labels.items():
                f2.write(key + '\t' + str(len(value)) + '\t' + ','.join(value) + '\n')

        print('\n# sequences:', seqs)
        print('# unique seq strings:', len(seq_to_labels))
        print('seq string with max count:', max_count)
