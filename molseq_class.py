"""
Author: Yun Zhang
Update: 11/30/2020
"""

import hashlib


class MolSeq(object):
    """
    for molecular sequences
        seq_label: a string representing the sequence.
                    Input seq_label could be an integer, which needs to be converted to a str.
        seq: a string representing the molecular sequence
    """

    """
    CONSTANTS
    """
    GAP              = '-'
    TERMINATE        = '*'

    AA_UNSPECIFIED   = 'X'
    NT_UNSPECIFIED   = 'N'

    # https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
    AA_AMBIGUOUS     = "XBZ"
    NT_AMBIGUOUS     = "NRYMKWSBDHV"

    # 22 standard amino acids, https://proteopedia.org/wiki/index.php/Amino_Acids
    AA_NONAMBIGUOUS  = "ARNDCQEGHILKMFPOUSTWYV-arndcqeghilkmfpoustwyv"
    DNA_NONAMBIGUOUS = "ACGT-acgt"
    RNA_NONAMBIGUOUS = "ACGU-acgu"

    AA_REGEXP        = "[^ARNDBCQEZGHILKMFPSTWYVUOXBZ\\-\\*]"
    DNA_REGEXP       = "[^ACGTRYMKWSNBDHV\\-\\*]"
    RNA_REGEXP       = "[^ACGURYMKWSNBDHV\\-\\*]"

    AA_SNP_STATES    = "ARNDCQEGHILKMFPOUSTWYVBZ-arndcqeghilkmfpoustwyvbz"
    DNA_SNP_STATES   = "ACGTRYMKWSBDHV-acgtrymkwsbdhv"
    RNA_SNP_STATES   = "ACGURYMKWSBDHV-acgurymkwsbdhv"

    def __init__(self, label, seq):
        self.__label = str(label).strip()
        self.__seq = str(seq).strip()

    def __str__(self):
        return self.to_fasta(60)

    def __len__(self):
        return self.get_length()

    def to_fasta(self, chars_per_line=60):
        """
        :param chars_per_line: an optional argument to limit numbers of chars per line.
                                False: no limit on chars_per_line
        :return: a string of fasta formatted sequence
        """
        if chars_per_line == False:
            result = ">" + str(self.__label) + "\n" + str(self.__seq)
        else:
            self.__fasta_line = ''
            for i in range(self.get_length()):
                if (i + 1) % chars_per_line == 0:
                    self.__fasta_line += self.__seq[i] + '\n'
                else:
                    self.__fasta_line += self.__seq[i]
            result = ">" + str(self.__label) + "\n" + str(self.__fasta_line)
        return result

    def get_label(self):
        return self.__label

    def set_label(self, new_label):
        self.__label = new_label

    def get_seq(self):
        return self.__seq

    def set_seq(self, new_seq):
        self.__seq = new_seq

    def get_length(self):
        return len(self.__seq)

    def get_residue_at(self, col):
        """
        :param col: 1-indexed
        :return: 0-indexed
        """
        if (col < 0 or col > self.get_length()):
            raise IndexError("invalid index")
        return self.__seq[col - 1]

    def get_coordinates(self, substring):
        """
        :param substring: a substring in __seq
        :return: 1-indexed, 1st occurrence
        """
        start = self.get_seq().find(substring)
        if start != -1:
            end = start + len(substring)
            start += 1
        else:
            end = -1
        return start, end

    def get_GC_percentage(self):
        residue_count = {}
        seq_string = self.get_seq()
        for residue in seq_string:
            count = residue_count.get(residue, 0)
            residue_count[residue] = count + 1

        gc_count = 0
        total_count = 0
        for key, value in residue_count.items():
            if key in ['C', 'c', 'G', 'g']:
                gc_count += value
            total_count += value
        return gc_count / total_count

    def trim_terminal_Ns(self):
        found_terminal_Ns = False
        original_length = self.get_length()
        self.__seq = self.__seq.strip('nN')
        if self.get_length() < original_length:
            found_terminal_Ns = True
        return found_terminal_Ns

    def trim_terminal_Xs(self):
        found_terminal_Xs = False
        original_length = self.get_length()
        self.__seq = self.__seq.strip('xX')
        if self.get_length() < original_length:
            found_terminal_Xs = True
        return found_terminal_Xs

    def has_excessive_consecutive_Ns(self, n_count):
        """
        :param n_count: min number of consecutive Ns or n's
        :return: True for passes, False otherwise
        """
        if n_count < 0:
            return False
        elif n_count == 0:
            return self.has_excessive_Ns(n_count)
        else:
            pattern = 'n' * n_count
            return pattern in self.get_seq().lower()

    def has_excessive_Ns(self, n_count):
        """
        :param n_count: min number of Ns or n's in sequence
        :return: True for count(N/n) != 0 and count(N/n) >= n_count, False otherwise
        """
        if n_count < 0:
            return False
        else:
            count = 0
            for c in self.__seq:
                if c == "N" or c == "n":
                    count += 1
                    if count >= n_count: # count >= 1
                        return True
            # count == 0 or count <= n_count
            return False

    def has_excessive_consecutive_Xs(self, x_count):
        """
        :param x_count: min number of consecutive Xs or x's
        :return: True for passes, False otherwise
        """
        if x_count < 0:
            return False
        elif x_count == 0:
            return self.has_excessive_Xs(x_count)
        else:
            pattern = 'x' * x_count
            return pattern in self.get_seq().lower()

    def has_excessive_Xs(self, x_count):
        """
        :param x_count: min number of Xs or x's in sequence
        :return: True for count(X/x) != 0 and count(X/x) >= x_count, False otherwise
        """
        if x_count < 0:
            return False
        else:
            count = 0
            for c in self.__seq:
                if c == "X" or c == "x":
                    count += 1
                    if count >= x_count: # count >= 1
                        return True
            # count == 0 or count <= n_count
            return False

    def has_min_length(self, min_length):
        """
        :param min_length: length requirement
        :return: True for passes, False otherwise
        """
        return self.get_length() >= min_length

    def passes_QC(self, min_length, n_count):
        """
        :param min_length: length requirement
        :param n_count: min number of consecutive Ns or n's
        :return: True for passes, False otherwise
        """
        return self.has_min_length(min_length) and not self.has_excessive_Ns(n_count)

    def get_lowercase_md5(self):
        return hashlib.md5(self.get_seq().lower().encode()).hexdigest()

    def get_md5(self):
        return hashlib.md5(self.get_seq().encode()).hexdigest()

    def get_lowercase_sha256(self):
        return hashlib.sha256(self.get_seq().lower().encode()).hexdigest()

    def get_sha256(self):
        return hashlib.sha256(self.get_seq().encode()).hexdigest()

    def get_sha512(self):
        return hashlib.sha512(self.get_seq().encode()).hexdigest()

    def get_reverse_complement(self):
        complement = {'A': 'T', 'T': 'A', 'U': 'A',
                      'C': 'G', 'G': 'C',
                      'M': 'K', 'K': 'M',
                      'R': 'Y', 'Y': 'R',
                      'V': 'B', 'B': 'V',
                      'H': 'D', 'D': 'H',
                      'N': 'N', 'W': 'W', 'S': 'S'}

        reverse_complement = []
        i = self.get_length() - 1
        while i >= 0:
            reverse_complement.append(complement.get(self.get_seq()[i]))
            i -= 1
        return "".join(reverse_complement)

    def get_no_gap_seq(self):
        seq_string = self.get_seq()
        seq_string_no_gap = ''
        for char in seq_string:
            if char != '-' and char != ' ':
                seq_string_no_gap += char
        return seq_string_no_gap