import os
import subprocess

from molseqcollection_class import MolSeqCollection


def peptides_to_mapped_aln_n_15mer(peptide_file, path, aln_files, all_accessions):
    # read in peptide file
    peptide_file_prefix = peptide_file.split('tsv')[0]
    with open(peptide_file, 'r') as infile, \
        open(peptide_file_prefix + 'mapped_aln_n_15mers.tsv', 'w') as outfile1, \
        open(peptide_file_prefix + 'mapped_aln_matrix.tsv', 'w') as outfile2, \
        open(peptide_file_prefix + 'mapped_15mer_matrix.tsv', 'w') as outfile3, \
        open(peptide_file_prefix + 'mapped_aln_heatmap.tsv', 'w') as outfile4, \
        open(peptide_file_prefix + 'mapped_15mer_heatmap.tsv', 'w') as outfile5:

        count = 0
        count_kmer = 0
        count_kmer_len_eq_15 = 0
        count_kmer_len_lt_15 = 0
        count_kmer_len_gt_15 = 0

        for line in infile:
            line = line.strip()

            if line:
                if line.startswith('#'):
                    # flat file
                    outfile1.write(line + '\tMapped_taxon\tAln_coords\tAligned_epitope\tMapped_genome_accession\tAligned_target\tAln_score\tBest_15mer\tBest_15mer_score\tAln_file\n')

                    # matrix files
                    outfile2.write(line + '\t' + 'alpha\t' * 16 + 'beta\t' * 15 + 'sarbeco\t' * 10 + '\n')
                    outfile2.write('\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\n')
                    outfile3.write(line + '\t' + 'alpha\t' * 16 + 'beta\t' * 15 + 'sarbeco\t' * 10 + '\n')
                    outfile3.write('\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\n')

                    # heatmap files
                    outfile4.write(line + '\talpha' * 16 + '\tbeta' * 15 + '\tsarbeco' * 10 + '\n')
                    outfile4.write('\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\t\n')
                    outfile5.write(line + '\talpha' * 16 + '\tbeta' * 15 + '\tsarbeco' * 10 + '\n')
                    outfile5.write('\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\t\n')
                else:
                    # non-header lines
                    CD4, Megapool, Start, Epitope_sequence, Peptide_ID, Num_times_tested, Num_times_positive, Total_magnitude, Genome_accession = line.split('\t')

                    # clean sequence str
                    Epitope_sequence.replace('\s', '')

                    # get mapped aln files
                    if Megapool.startswith('nsp'):
                        mapped_files = aln_files[CD4 + '_orf1ab']
                    else:
                        mapped_files = aln_files[CD4 + '_' + Megapool]

                    # each row/peptide_ID has 3 mapped aln files
                    peptide_to_seqs_identities = {}
                    for file in mapped_files:
                        mapped_taxon = file.split('reps')[1].split('.')[1]

                        # read in aln file
                        mapped_aln_file = os.path.join(path, file)
                        aln = MolSeqCollection.parse(mapped_aln_file, is_aligned=True, allow_duplicate_labels=False)
                        mapped_aln_coords = aln.get_aln_coords_for_substring(seq_label=Genome_accession, query=Epitope_sequence)

                        # get matched aln region
                        if mapped_aln_coords != (-1, -1):
                            print(mapped_aln_file)
                            print(Genome_accession)
                            mapped_aln_region = aln.get_region(mapped_aln_coords[0], mapped_aln_coords[1])
                            aligned_epitope = mapped_aln_region.get_seqs_by_label_substring(delimiter='|', substring_index=2, label_list=[Genome_accession]).get_molseq_list()[0].get_seq()

                            # get mapped sequence in each sequence in aln file
                            for seq in mapped_aln_region.get_molseq_list():
                                matched_label = seq.get_label()
                                matched_accession = matched_label.split('|')[2]

                                # exclude KX512809 which is missing orf1ab
                                if 'KX512809' not in matched_label:
                                    # get target seq
                                    aligned_target = seq.get_seq()
                                    aligned_target_full = aln.get_seqs_by_labels([matched_label]).get_molseq_list()[0].get_seq()
                                    kmer_seed = seq.get_no_gap_seq()

                                    # calculate alignment score
                                    aln_score = calculate_aln_score(aligned_epitope, aligned_target)

                                    # find best_15mer with max 15mer_score

                                    if len(kmer_seed) == 15:
                                        # case 1
                                        # search space: kmer_seed
                                        # best_15mer_left = 1
                                        best_15mer_score = calculate_kmer_score(Epitope_sequence, kmer_seed)
                                        best_15mer = kmer_seed
                                        count_kmer_len_eq_15 += 1
                                    elif len(kmer_seed) > 15:
                                        # case 2
                                        # search space: kmer_seed
                                        best_15mer, best_15mer_score = find_best_15mer(Epitope_sequence, kmer_search_space=kmer_seed)
                                        count_kmer_len_lt_15 += 1
                                    else:
                                        # case 3
                                        # kmer_seed length < 15:
                                        # search space: (15-k) padding + kmer_seed + (15-k) padding
                                        # 0-based in str
                                        kmer_search_space = find_kmer_search_space(aligned_target_full, mapped_aln_coords, kmer_seed, k=15)
                                        best_15mer, best_15mer_score = find_best_15mer(Epitope_sequence, kmer_search_space)
                                        count_kmer_len_gt_15 += 1

                                    count_kmer += 1
                                    # save seqs and identities to map
                                    peptide_to_seqs_identities[matched_accession] = [aligned_target,
                                                                                     best_15mer,
                                                                                     aln_score,
                                                                                     best_15mer_score]

                                    outfile1.write(line + '\t' + mapped_taxon + '\t' + str(
                                        mapped_aln_coords) + '\t' + aligned_epitope + '\t' + matched_accession + '\t' + aligned_target + '\t' + str(aln_score) + '\t' + best_15mer + '\t' + str(best_15mer_score) + '\t' + file + '\n')

                        else:
                            outfile1.write(line + '\t\t\t\t\t\t\t\t\n')

                    # done with all files for this line
                    values2 = ''
                    values3 = ''
                    values4 = ''
                    values5 = ''
                    for accession in all_accessions:
                        if accession in peptide_to_seqs_identities.keys():
                            # mapped aln seq
                            values2 += '\t' + peptide_to_seqs_identities[accession][0]
                            # 15-mer
                            values3 += '\t' + peptide_to_seqs_identities[accession][1]

                            if len(peptide_to_seqs_identities[accession][1]) < 15:
                                print('< 15-mer', CD4, Megapool, Start, Epitope_sequence, Peptide_ID)

                            # mapped aln seq identity
                            values4 += '\t' + str(peptide_to_seqs_identities[accession][2])
                            # 15-mer identity
                            values5 += '\t' + str(peptide_to_seqs_identities[accession][3])
                        else:
                            values2 += '\t'
                            values3 += '\t'
                            values4 += '\t'
                            values5 += '\t'

                    outfile2.write(line + values2 + '\n')
                    outfile3.write(line + values3 + '\n')
                    outfile4.write(line + values4 + '\n')
                    outfile5.write(line + values5 + '\n')

                    # done with this line
                    count += 1
    print('# epitopes read    :', count)
    print("# kmer_len_eq_15   :", count_kmer_len_eq_15)
    print("# kmer_len_lt_15   :", count_kmer_len_lt_15)
    print("# kmer_len_gt_15   :", count_kmer_len_gt_15)
    print("# kmer             :", count_kmer)
    return


def peptides_to_mapped_aln_n_kmer(peptide_file, aln_files, all_accessions, mat_peptide=True):
    """

    :param peptide_file: sent by John Sidney
    #Sequence_in_source	Sequence_in_NC_045512	Len	Antigen	Start	End	Mapped_start	Mapped_end	Identity_w_Wuhan_ref	CD4_cluster_ID	Restriction	Genome_accession
    :param aln_files_path:
    :param aln_files:
    :param all_accessions:
    :return:
    """
    # read in peptide file
    peptide_file_prefix = peptide_file.split('tsv')[0]
    with open(peptide_file, 'r') as infile, \
        open(peptide_file_prefix + 'mapped_aln_n_kmers.tsv', 'w') as outfile1, \
        open(peptide_file_prefix + 'mapped_aln_matrix.tsv', 'w') as outfile2, \
        open(peptide_file_prefix + 'mapped_kmer_matrix.tsv', 'w') as outfile3, \
        open(peptide_file_prefix + 'mapped_aln_heatmap.tsv', 'w') as outfile4, \
        open(peptide_file_prefix + 'mapped_kmer_heatmap.tsv', 'w') as outfile5:

        count = 0
        count_kmer = 0
        count_seed_len_eq_epitope_len = 0
        count_seed_len_lt_epitope_len = 0
        count_seed_len_gt_epitope_len = 0

        for line in infile:
            line = line.strip()

            if line:
                if line.startswith('#'):
                    # flat file
                    outfile1.write(line + '\tMapped_taxon\tAln_coords\tAligned_epitope\tMapped_genome_accession\tAligned_target\tAln_score\tBest_kmer\tBest_kmer_score\tAln_file\n')

                    # matrix files
                    outfile2.write(line + '\t' + 'alpha\t' * 16 + 'beta\t' * 15 + 'SARS\t' * 9 + 'SARS\n')
                    outfile2.write('\t\t\t\t\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\n')
                    outfile3.write(line + '\t' + 'alpha\t' * 16 + 'beta\t' * 15 + 'SARS\t' * 9 + 'SARS\n')
                    outfile3.write('\t\t\t\t\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\n')

                    # heatmap files
                    outfile4.write(line + '\talpha' * 16 + '\tbeta' * 15 + '\tSARS' * 10 + '\n')
                    outfile4.write('\t\t\t\t\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\n')
                    outfile5.write(line + '\talpha' * 16 + '\tbeta' * 15 + '\tSARS' * 10 + '\n')
                    outfile5.write('\t\t\t\t\t\t\t\t\t\t\t\t\t' + '\t'.join(all_accessions) + '\n')
                else:
                    # non-header lines
                    Source_virus, Sequence_in_source, Sequence_in_NC_045512, Len, Antigen, Start, End, Mapped_start, Mapped_end, Identity_w_Wuhan_ref, CD4_cluster_ID, Restriction, Genome_accession = line.split('\t')

                    # clean sequence str
                    Sequence_in_NC_045512.replace('\s', '')

                    # convert len
                    epitope_length = int(Len)

                    # get mapped aln files
                    if not mat_peptide:
                        if Antigen.lower().startswith('nsp'):
                            mapped_files = aln_files[Source_virus + '_orf1ab']
                        elif Antigen.lower() in ['e', 'm', 'n', 's']:
                            mapped_files = aln_files[Source_virus + '_' + Antigen]
                        else:
                            continue
                    else:
                        if Antigen.lower().startswith('nsp') or Antigen.lower() in ['e', 'm', 'n', 's']:
                            mapped_files = aln_files[Source_virus + '_' + Antigen]
                        else:
                            continue # ignore this line

                    # each row/peptide_ID has 3 mapped aln files
                    peptide_to_seqs_identities = {}
                    for file in mapped_files:
                        mapped_taxon = file.split('reps')[1].split('.')[1]

                        # read in aln file
                        mapped_aln_file = os.path.join(path, file)
                        aln = MolSeqCollection.parse(mapped_aln_file, is_aligned=True, allow_duplicate_labels=False)
                        mapped_aln_coords = aln.get_aln_coords_for_substring(seq_label=Genome_accession, query=Sequence_in_NC_045512)

                        # get matched aln region
                        if mapped_aln_coords != (-1, -1):
                            mapped_aln_region = aln.get_region(mapped_aln_coords[0], mapped_aln_coords[1])
                            aligned_epitope = mapped_aln_region.get_seqs_by_label_substring(delimiter='|', substring_index=2, label_list=[Genome_accession]).get_molseq_list()[0].get_seq()

                            # get mapped sequence in each sequence in aln file
                            for seq in mapped_aln_region.get_molseq_list():
                                matched_label = seq.get_label()
                                matched_accession = matched_label.split('|')[2]

                                # exclude KX512809 which is missing orf1ab
                                if 'KX512809' not in matched_label:
                                    # get target seq
                                    aligned_target = seq.get_seq()
                                    aligned_target_full = aln.get_seqs_by_labels([matched_label]).get_molseq_list()[0].get_seq()
                                    seed = seq.get_no_gap_seq()

                                    # calculate alignment score
                                    aln_score = calculate_aln_score(aligned_epitope, aligned_target)

                                    # find best_15mer with max 15mer_score

                                    if len(seed) == epitope_length:
                                        # case 1
                                        # search space: seed
                                        # best_kmer_left = 1
                                        best_kmer_score = calculate_kmer_score(Sequence_in_source, seed)
                                        best_kmer = seed
                                        count_seed_len_eq_epitope_len += 1
                                    elif len(seed) > epitope_length:
                                        # case 2
                                        # search space: kmer_seed
                                        best_kmer, best_kmer_score = find_best_kmer(Sequence_in_source, kmer_search_space=seed, k=epitope_length)
                                        count_seed_len_lt_epitope_len += 1
                                    else:
                                        # case 3
                                        # seed length < epitope_length:
                                        # search space: (epitope_length - seed_length) padding + seed + (epitope_length - seed_length) padding
                                        # 0-based in str
                                        kmer_search_space = find_kmer_search_space(aligned_target_full, mapped_aln_coords, seed, k=epitope_length)
                                        best_kmer, best_kmer_score = find_best_kmer(Sequence_in_source, kmer_search_space, k=epitope_length)
                                        count_seed_len_gt_epitope_len += 1

                                    count_kmer += 1
                                    # save seqs and identities to map
                                    peptide_to_seqs_identities[matched_accession] = [aligned_target,
                                                                                     best_kmer,
                                                                                     aln_score,
                                                                                     best_kmer_score]

                                    outfile1.write(line + '\t' + mapped_taxon + '\t' + str(
                                        mapped_aln_coords) + '\t' + aligned_epitope + '\t' + matched_accession + '\t' + aligned_target + '\t' + format(aln_score, '.3f') + '\t' + best_kmer + '\t' + format(best_kmer_score, '.3f') + '\t' + file + '\n')

                        else:
                            outfile1.write(line + '\t\t\t\t\t\t\t\t\n')

                    # done with all files for this line
                    values2 = ''
                    values3 = ''
                    values4 = ''
                    values5 = ''
                    for accession in all_accessions:
                        if accession in peptide_to_seqs_identities.keys():
                            # mapped aln seq
                            values2 += '\t' + peptide_to_seqs_identities[accession][0]
                            # k-mer
                            values3 += '\t' + peptide_to_seqs_identities[accession][1]

                            # mapped aln seq identity
                            values4 += '\t' + format(peptide_to_seqs_identities[accession][2], '.3f')
                            # k-mer identity
                            values5 += '\t' + format(peptide_to_seqs_identities[accession][3], '.3f')
                        else:
                            values2 += '\t'
                            values3 += '\t'
                            values4 += '\t'
                            values5 += '\t'

                    outfile2.write(line + values2 + '\n')
                    outfile3.write(line + values3 + '\n')
                    outfile4.write(line + values4 + '\n')
                    outfile5.write(line + values5 + '\n')

                    # done with this line
                    count += 1
    print('# epitopes read                :', count)
    print("# count_seed_len_eq_epitope_len:", count_seed_len_eq_epitope_len)
    print("# count_seed_len_lt_epitope_len:", count_seed_len_lt_epitope_len)
    print("# count_seed_len_gt_epitope_len:", count_seed_len_gt_epitope_len)
    print("# kmer                         :", count_kmer)
    return


def get_no_gap_length(str):
    count = 0
    for c in str:
        if c not in ['-', ' ']:
            count += 1
    return count


def get_no_gap_seq(str):
    result = ''
    for c in str:
        if c not in ['-', ' ']:
            result += c
    return result


def find_kmer_search_space(aligned_target_full, mapped_aln_coords, seed, k):
    """
    search space: [left, right] inclusive
    left:  can be real residue or '-' (string index 0)
    right: can be real residue or '-' (string index n-1)
	to_pad = k - seed_length, seed excludes '-'
	left = mapped_aln_coords[0] - 1 (0-based), move left (k - seed_length) positions skipping ‘-’, or left = 0
	right= mapped_aln_coords[1] - 1 (0-based), move right (k - seed_length) positions skipping ‘-’, or right = n-1
	stop condition: to_pad = 0 or left = 0 / right = n - 1 (whether to_pad is 0 or not does not matter)

    :param aligned_target_full:
    :param mapped_aln_coords:
    :param seed:
    :param k: resulting k-mer size
    :return: search space of sequence string, with '-' removed
    """
    # search space: (k - seed_length) padding + seed + (k - seed_length) padding
    # 0-based in str
    left = mapped_aln_coords[0] - 1
    right = mapped_aln_coords[1] - 1

    # find left boundary
    to_pad = k - len(seed)
    while to_pad > 0 and left > 0:
        if aligned_target_full[left - 1] != '-':
            to_pad -= 1
        left -= 1
    # to_pad = 0 or left = 0 (to_pad is 0 or < 0 does not matter, because left bound is reached)

    to_pad = k - len(seed)
    while to_pad > 0 and right < len(aligned_target_full) - 1:
        if aligned_target_full[right + 1] != '-':
            to_pad -= 1
        right += 1
    # to_pad = 0 or right = n - 1

    return get_no_gap_seq(aligned_target_full[left: right + 1])


def calculate_aln_score(aligned_epitope, aligned_target):
    # count matched residues in mapped_aln_region
    matched_cols = 0
    for i in range(len(aligned_epitope)):
        if aligned_epitope[i] == aligned_target[i] and aligned_epitope[i] != '-':
            matched_cols += 1
    aln_score = matched_cols / get_no_gap_length(aligned_epitope)
    return aln_score


def calculate_kmer_score(epitope_sequence, target_sequence):
    """
    assumption: 2 sequences are of equal length

    **************012345678901234**************
    abcdefghijklmno
    l             r
                                abcdefghijklmno
                                l             r

    # deprecated
    left boundary
    epitope          0123456789
    target1 abcdefghij
    columns          + 0:9

    offset in target
    epitope  0123456789
    target1 abcdefghij
    columns  +++++++++ 0-8:1-9

    no offset
    epitope 0123456789
    target1 abcdefghij
    columns ++++++++++ 0-9:0-9

    offset in epitope
    epitope 0123456789
    target1  abcdefghij
    columns  +++++++++ 1-9:0:8

    right boundary
    epitope 0123456789
    target1          abcdefghij
    columns          + 9:0

    :return:
    """
    if len(epitope_sequence) != len(target_sequence):
        print("expected two input sequences to be of equal length")
        return -1
    else:
        padding = '*' * (len(epitope_sequence) - 1)
        epitope_sequence_padded = padding + epitope_sequence + padding

        max_matched_cols = 0
        max_matched_offset = 0 # offset in epitope_sequence_padded
        offset = 0

        while offset < len(epitope_sequence_padded) - len(target_sequence) + 1:
            matched_cols = 0
            # a window of target sequence size
            for i in range(len(target_sequence)):
                if target_sequence[i] == epitope_sequence_padded[offset + i]:
                    matched_cols += 1
            # done with this window
            if matched_cols > max_matched_cols:
                max_matched_cols = matched_cols
                max_matched_offset = offset
            offset += 1

        '''
        # index in string
        max_matched_epitope_left = 0
        for offset in range(0, len(epitope_sequence)): # [0, len(epitope_sequence) - 1]
            matched_cols = 0
            # window range
            for i in range(offset, len(epitope_sequence)): # [offset, len(epitope_sequence) - 1]
                if epitope_sequence[i] == target_sequence[i - offset]: #
                    matched_cols += 1
            # done with this window
            if matched_cols > max_matched_cols:
                max_matched_cols = matched_cols
                max_matched_offset = offset

        # offset in target
        for offset in range(len(target_sequence) - 1, 0, -1): # [len(target_sequence) - 1, 1]
            matched_cols = 0
            # window range
            for i in range(offset, len(target_sequence)): # [offset, len(target_sequence) - 1]
                if target_sequence[i] == epitope_sequence[i - offset]: #
                    matched_cols += 1
            # done with this window
            if matched_cols > max_matched_cols:
                max_matched_cols = matched_cols
                max_matched_offset = -offset
        '''
    #print("max_matched_offset:", max_matched_offset)
    return max_matched_cols / len(epitope_sequence)


def find_best_kmer(epitope_sequence, kmer_search_space, k):
    # 0-based in str
    left = 0
    best_kmer_left = 0
    best_kmer_score = 0
    best_kmer = ""

    while left + k - 1 < len(kmer_search_space):
        curr_kmer_score = calculate_kmer_score(epitope_sequence, kmer_search_space[left: left + k]) # str[left to right]
        if curr_kmer_score > best_kmer_score:
            best_kmer_score = curr_kmer_score
            best_kmer_left = left
        left += 1

    best_kmer = kmer_search_space[best_kmer_left: best_kmer_left + k]
    return best_kmer, best_kmer_score


def find_best_15mer(epitope_sequence, kmer_search_space):
    # 0-based in str
    left = 0
    best_15mer_left = 0
    best_15mer_score = 0
    best_15mer = ""

    while left + 14 < len(kmer_search_space):
        curr_15mer_score = calculate_kmer_score(epitope_sequence, kmer_search_space[left: left + 15]) # str[left to right]
        if curr_15mer_score > best_15mer_score:
            best_15mer_score = curr_15mer_score
            best_15mer_left = left
        left += 1

    best_15mer = kmer_search_space[best_15mer_left: best_15mer_left + 15]
    return best_15mer, best_15mer_score


def get_nsp_seq_strings(wuhan_mat_peptides_to_nsps, protein_file='/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_rep_viruses/reps_selected/Alphabetacoronavirus_Protein20210617.LJI_reps20210910.fasta',
                   accession='NC_045512'):
    wuhan_nsp_to_seq_str = {}
    seqs = MolSeqCollection.parse(protein_file, is_aligned=False)
    wuhan_proteins = seqs.get_seqs_by_label_find_substring(accession)
    for key, value in wuhan_mat_peptides_to_nsps.items():
        wuhan_nsp_to_seq_str[value] = wuhan_proteins.get_seqs_by_label_find_substring(key).get_molseq_list()[
            0].get_seq()
    return wuhan_nsp_to_seq_str


def orf1ab_aln_to_nsp_alns(orf1ab_einsi_files, wuhan_nsp_to_seq_str):
    # orf1ab aln to nsp alns
    path = '/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_rep_viruses/SARS2_Sidney_aligned/'

    #orf1ab_einsi_files = [path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.fasta',
    #                path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.fasta',
    #                path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.fasta']

    for file in orf1ab_einsi_files:
        aln = MolSeqCollection.parse(file, is_aligned=True)
        for label, seq_str in wuhan_nsp_to_seq_str.items():
            aln_start, aln_end = aln.get_aln_coords_for_substring(seq_label='NC_045512', query=seq_str)
            # end of nsp16, non-wuhan may have insertions
            if label == 'nsp16':
                aln_end = aln.get_length()
            print(file, label, aln_start, aln_end)

            if aln_start != -1 and aln_end != -1:
                aln.get_region(aln_start, aln_end).write_to_file(file.split('fasta')[0] + label + '.fasta')
    return


if __name__ == '__main__':
    peptide_file = '/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_peptides/LJI_CCC_Epitopes20210913.tsv'
    path = '/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_reps/aligned/'
    aln_files_einsi = {'NL63_orf1ab': [
        'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.einsi.fasta',
        'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_NL63_NC_005831.einsi.fasta',
        'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.plus_NL63_NC_005831.einsi.fasta'],
                       'NL63_M': [
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.alpha.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.beta.plus_NL63_NC_005831.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.sarbeco.plus_NL63_NC_005831.einsi.fasta'],
                       'NL63_N': [
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.alpha.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.beta.plus_NL63_NC_005831.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.sarbeco.plus_NL63_NC_005831.einsi.fasta'],
                       'NL63_S': [
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.alpha.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.beta.plus_NL63_NC_005831.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.sarbeco.plus_NL63_NC_005831.einsi.fasta'],
                       'OC43_orf1ab': [
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_OC43_NC_006213.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.plus_OC43_NC_006213.einsi.fasta'],
                       'OC43_M': [
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.alpha.plus_OC43_NC_006213.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.beta.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.sarbeco.plus_OC43_NC_006213.einsi.fasta'],
                       'OC43_N': [
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.alpha.plus_OC43_NC_006213.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.beta.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.sarbeco.plus_OC43_NC_006213.einsi.fasta'],
                       'OC43_S': [
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.alpha.plus_OC43_NC_006213.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.beta.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.sarbeco.plus_OC43_NC_006213.einsi.fasta'],
                       'OC43_HE': ['LJI_reps_proteins_HE.reps.beta.einsi.fasta'],
                       'SARS2_orf1ab': ['Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.fasta'],
                       'SARS2_E': ['Alphabetacoronavirus_Protein20210617.metadata_filtered.E_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.E_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.E_passed_QC.plus_must_include.reps.sarbeco.einsi.fasta'],
                       'SARS2_M': ['Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.sarbeco.einsi.fasta'],
                       'SARS2_N': ['Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.sarbeco.einsi.fasta'],
                       'SARS2_S': ['Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.fasta',
                           'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.sarbeco.einsi.fasta'],
                       'SARS2_NSP1': ['Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp1.fasta',
                                      'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp1.fasta',
                                      'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp1.fasta'],
                        'SARS2_NSP2': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp2.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp2.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp2.fasta'],
                        'SARS2_NSP3': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp3.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp3.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp3.fasta'],
                        'SARS2_NSP4': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp4.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp4.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp4.fasta'],
                        'SARS2_NSP5': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp5.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp5.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp5.fasta'],
                        'SARS2_NSP6': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp6.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp6.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp6.fasta'],
                        'SARS2_NSP7': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp7.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp7.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp7.fasta'],
                        'SARS2_NSP8': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp8.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp8.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp8.fasta'],
                        'SARS2_NSP9': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp9.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp9.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp9.fasta'],
                        'SARS2_NSP10': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp10.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp10.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp10.fasta'],
                        'SARS2_NSP12': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp12.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp12.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp12.fasta'],
                        'SARS2_NSP13': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp13.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp13.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp13.fasta'],
                        'SARS2_NSP14': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp14.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp14.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp14.fasta'],
                        'SARS2_NSP15': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp15.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp15.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp15.fasta'],
                        'SARS2_NSP16': [
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.nsp16.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.nsp16.fasta',
                            'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.nsp16.fasta']
    }

    alpha20210910 = ['NC_005831', 'NC_002645', 'NC_028752', 'NC_009657', 'NC_018871', 'MH687935', 'NC_009988', 'NC_028824', 'NC_010437', 'NC_010438', 'KJ473798', 'NC_048216', 'NC_022103', 'NC_046964', 'NC_028814', 'MK720945']
    beta = ['NC_019843', 'NC_038294', 'NC_003045', 'NC_006213', 'NC_006577', 'NC_039207', 'MK211374', 'MT121216', 'NC_009019', 'NC_009020', 'MT337386', 'KC869678', 'MG596802', 'HM211100', 'HM211101']
    sarbeco = ['NC_045512', 'MW636737', 'NC_004718', 'MT706389', 'MW868471', 'MT745698', 'MW731141', 'MW848086', 'MW725757', 'MT952134']
    all_accessions20210910 = alpha20210910 + beta + sarbeco
    #peptides_to_mapped_aln_n_15mer(peptide_file, path, aln_files_einsi, all_accessions20210910)


    ########
    # append rep seq to other groups
    path = '/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_rep_viruses/not_aligned/'

    protein_files = [
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.fasta'
    ]
    
    protein_files = [
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.E_passed_QC.plus_must_include.reps.alpha.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.E_passed_QC.plus_must_include.reps.beta.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.E_passed_QC.plus_must_include.reps.sarbeco.fasta'
    ]
    
    protein_files = [
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.alpha.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.beta.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.M_passed_QC.plus_must_include.reps.sarbeco.fasta'
    ]

    protein_files = [
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.alpha.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.beta.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.N_passed_QC.plus_must_include.reps.sarbeco.fasta'
    ]
    
    protein_files = [
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.alpha.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.beta.fasta',
        path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC.plus_must_include.reps.sarbeco.fasta'
    ]
    
    # list of molseqcollection objects [alpha, beta, sarbeco]
    seqs_list = []

    # parse seq files
    for file in protein_files:
        seqs = MolSeqCollection.parse(file, is_aligned=False)
        seqs_list.append(seqs)

    # extract refs
    # alpha
    NL63_NC_005831 = seqs_list[0].get_seqs_by_label_substring(delimiter='|', substring_index=2,
                                                              label_list=['NC_005831']).get_molseq_list()[0]
    # beta
    OC43_NC_006213 = seqs_list[1].get_seqs_by_label_substring(delimiter='|', substring_index=2,
                                                              label_list=['NC_006213']).get_molseq_list()[0]
    # sars2
    SARS2_NC_045512 = seqs_list[2].get_seqs_by_label_substring(delimiter='|', substring_index=2,
                                                              label_list=['NC_045512']).get_molseq_list()[0]
    # append sars2 ref to alpha and beta groups
    seqs_list[0].add_molseq(SARS2_NC_045512)
    seqs_list[0].write_to_file(protein_files[0].split('.fasta')[0] + ".plus_SARS2_NC_045512.fasta")
    seqs_list[1].add_molseq(SARS2_NC_045512)
    seqs_list[1].write_to_file(protein_files[1].split('.fasta')[0] + ".plus_SARS2_NC_045512.fasta")
    

    ########
    # align reps + SARS2_NC_045512 refs
    directory = os.fsencode(path)

    aln_mode = 'einsi'
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('SARS2_NC_045512.fasta') or filename.endswith('sarbeco.fasta'):
            print(file) # b'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.sarbeco.plus_must_include.fasta'
            infile = os.path.join(path, filename)
            outfile = open(infile.split('fasta')[0] + aln_mode + '.fasta', 'w')
            subprocess.call([aln_mode, '--thread', '-1', '--reorder', '--quiet', infile], stdout=outfile)

            if filename.startswith('Alphabetacoronavirus_Protein20210617.metadata_filtered.S_passed_QC'):
                outfile = open(infile.split('fasta')[0] + 'dash.fasta', 'w')
                subprocess.call(['mafft', '--dash', '--localpair', '--thread', '-1', '--reorder', '--quiet', infile], stdout=outfile)


    ########
    # protein names in fasta to nsps in epitopes.tsv
    wuhan_mat_peptides_to_nsps = {
        "2'_O_ribose_methyltransferase": 'nsp16',
        "nsp11": 'nsp11',
        "nsp4": 'nsp4',
        "3'_to_5'_exonuclease": 'nsp14',
        "nsp10": 'nsp10',
        "endoRNAse": 'nsp15',
        "3C_like_proteinase": 'nsp5',
        "nsp9": 'nsp9',
        "nsp2": 'nsp2',
        "nsp3": 'nsp3',
        "leader_protein": 'nsp1',
        "nsp8": 'nsp8',
        "nsp7": 'nsp7',
        "helicase": 'nsp13',
        "nsp6": 'nsp6',
        "RNA_dependent_RNA_polymerase": 'nsp12'
    }

    ########
    # orf1ab to nsp files

    # get nsp seq strs
    wuhan_nsp_to_seq_str = get_nsp_seq_strings(wuhan_mat_peptides_to_nsps,
                                              protein_file='/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_rep_viruses/reps_selected/Alphabetacoronavirus_Protein20210617.LJI_reps20210910.fasta',
                                              accession='NC_045512')
    print(wuhan_nsp_to_seq_str)

    # orf1ab aln to nsp alns
    path = '/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_rep_viruses/SARS2_Sidney_aligned/'
    orf1ab_einsi = [path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.alpha.plus_SARS2_NC_045512.einsi.fasta',
                    path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.beta.plus_SARS2_NC_045512.einsi.fasta',
                    path + 'Alphabetacoronavirus_Protein20210617.metadata_filtered.orf1ab_passed_QC.plus_must_include.reps.sarbeco.einsi.fasta']

    orf1ab_aln_to_nsp_alns(orf1ab_einsi, wuhan_nsp_to_seq_str)

    ########
    peptide_file = '/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_epitopes/SARS-CoV-2_epitopes.tsv'
    path = '/Users/yuzhang/OneDrive - J. Craig Venter Institute/2019nCov/LJI_Pan-Coronavirus_Vaccines/LJI_rep_viruses/SARS2_Sidney_aligned/'
    peptides_to_mapped_aln_n_kmer(peptide_file, aln_files_einsi, all_accessions20210910, mat_peptide=True)