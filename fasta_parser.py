import re
from molseq_class import MolSeq

# using uppercase to indicate constant
LABEL_RE = re.compile(">\s*(.+)")
GAP_RE = re.compile("[-\s]+")
#GAP_RE = re.compile("[-\h]+")
WS_RE = re.compile("[\s]+")


def parse_fasta_file(file, remove_gaps=False):
    with open(file, 'r') as infile:
        # memorize what has been read in current_label and current_seq
        current_label = ''
        current_seq = ''
        MolSeq_list = []
        count_seqs_with_gap = 0
        found_gap = False

        for line in infile:
            line = line.strip()

            if line:
                if line.startswith('>'):
                    # if this is a new seq but not the 1st seq, create a MolSeq object of the previous seq
                    if current_seq:
                        MolSeq_list.append(MolSeq(current_label, current_seq))
                        if found_gap:
                            count_seqs_with_gap += 1

                        # then reset current_label and current_seq
                        current_label = LABEL_RE.search(line).group(1)
                        current_seq = ''

                    # for the 1st seq in the file
                    else:
                        current_label = LABEL_RE.search(line).group(1)

                    # reset state of found_gap
                    found_gap = False

                else:
                    if remove_gaps:
                        original_line = line
                        original_length = len(line)
                        line = re.sub(GAP_RE, "", line)
                        if len(line) < original_length:
                            found_gap = True
                            print("gap found in:", current_label)#, original_line)
                    else:
                        line = re.sub(WS_RE, "", line)
                    current_seq += line

        # for the last seq in the file
        MolSeq_list.append(MolSeq(current_label, current_seq))
        if found_gap:
            count_seqs_with_gap += 1

        if count_seqs_with_gap > 0:
            print("\nseqs with gap found", count_seqs_with_gap)
        return MolSeq_list, count_seqs_with_gap


if __name__ == '__main__':
    MolSeq_list1 = parse_fasta_file('test.fasta')
    for seq in MolSeq_list1:
        print(seq.to_fasta())

    MolSeq_list2 = parse_fasta_file('test.fasta', remove_gaps=True)
    for seq in MolSeq_list2:
        print(seq.to_fasta())
