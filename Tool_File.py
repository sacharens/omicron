import pandas as pd
import os
import re
from Bio import SeqIO
import concurrent.futures
import numpy as np
from scipy.stats import entropy
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import subprocess
import os
import re
import time
import random
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# ----------------------------------------------------------------------------------------------------------------------
# clean up fasta file from unknowns and devied into months
# ----------------------------------------------------------------------------------------------------------------------

def clean_x_from_fasta_nad_month_seperation(file_and_path, save_path, file_name_for_save, ref_seq_f):
    """
    this function will remove all sequences with amino acid x (unknown)
    :param file_and_path: full path with folder and it must contain the file extension
    :return: new file to the path sent in the second arg
    """
    seq_mot_in_range = 0
    all_seq_count = 0
    more_than_one_x = 0
    more_than_five = 0
    list_for_save_befor_ago = []
    list_for_save_ago = []
    list_for_save_sep = []
    list_for_save_oct = []
    list_for_save_nov = []
    count_X_seq = 0
    seq_with_Illegal_character = 0
    range_val = len(ref_seq_f) * 0.025
    low_limit = len(ref_seq_f) - range_val
    hight_limit = len(ref_seq_f) + range_val

    input_handle = open(file_and_path, "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        all_seq_count += 1
        date = find_date(record.description)
        NoneType = type(None)

        if record.seq.find('X') == -1 and low_limit < len(record.seq) < hight_limit:
            if type(date) != NoneType:

                if date < datetime.datetime.strptime("2021-08-01".strip(), "%Y-%m-%d"):
                    list_for_save_befor_ago.append(record)
                if datetime.datetime.strptime("2021-08-01".strip(), "%Y-%m-%d") < date < datetime.datetime.strptime(
                        "2021-09-01".strip(), "%Y-%m-%d"):
                    list_for_save_ago.append(record)
                if datetime.datetime.strptime("2021-09-01".strip(), "%Y-%m-%d") < date < datetime.datetime.strptime(
                        "2021-10-01".strip(), "%Y-%m-%d"):
                    list_for_save_sep.append(record)
                if datetime.datetime.strptime("2021-10-01".strip(), "%Y-%m-%d") < date < datetime.datetime.strptime(
                        "2021-11-01".strip(), "%Y-%m-%d"):
                    list_for_save_oct.append(record)
                if datetime.datetime.strptime("2021-11-01".strip(), "%Y-%m-%d") < date < datetime.datetime.strptime(
                        "2021-12-01".strip(), "%Y-%m-%d"):
                    list_for_save_nov.append(record)
                else:
                    continue
        else:
            if len(record.seq) < low_limit or len(record.seq) > hight_limit:
                seq_mot_in_range += 1
        if record.seq.find('X') != -1:
            count_X_seq += 1
            if record.seq.count('X') < 2:
                more_than_one_x += 1
            if 1 < record.seq.count('X') < 6:
                more_than_five += 1

    input_handle.close()

    # discarding x
    SeqIO.write(list_for_save_befor_ago, save_path + file_name_for_save + '_befor_ago.fa', 'fasta')
    SeqIO.write(list_for_save_ago, save_path + file_name_for_save + 'ago.fa', 'fasta')
    SeqIO.write(list_for_save_sep, save_path + file_name_for_save + 'sep.fa', 'fasta')
    SeqIO.write(list_for_save_oct, save_path + file_name_for_save + 'oct.fa', 'fasta')
    SeqIO.write(list_for_save_nov, save_path + file_name_for_save + 'nov.fa', 'fasta')

    print('Number of sequences in fasta file :', all_seq_count)
    print('Number of sequences that include X :', count_X_seq)
    print('Number of sequences that are not in range:', seq_mot_in_range)
    print('Number of sequences that include only X one :', more_than_one_x)
    print('Number of sequences that includemore than one x and less than five :', more_than_five)


#                           !!!!!!!!!!!!!!!!!!!!!!!! for states now !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# ----------------------------------------------------------------------------------------------------------------------
# deal with every country's mutations by month ###### for states now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ----------------------------------------------------------------------------------------------------------------------


def thread_for_state(date_list, path_and_file, country, save_path_perant, ref_sequence):
    # for country in list_of_countries:
    print('starting threading  {}....'.format(country))  # --------------------------------------------------------start
    try:
        os.mkdir(save_path_perant + country)
    except OSError as error:
        print(error)

    save_path = save_path_perant + country + '/'
    print(save_path)
    ##########################  creating fasta file with seq from country  ################################
    list_for_save = []
    input_handle = open(path_and_file, "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        # find the date of the sequence
        location = find_seq_states(record.description)
        if location == country:
            list_for_save.append(record)

    input_handle.close()
    SeqIO.write(list_for_save, save_path + country + '_sequences.fa', 'fasta')

    column_list = ['position', 'mutated aa']
    # column_list.extend(date_list)
    the_mane_df = pd.DataFrame(columns=column_list)

    for month in date_list:
        print(month)
        list_for_save = []
        input_handle = open(save_path + country + '_sequences.fa', "r")
        for record in SeqIO.parse(input_handle, "fasta"):
            # find the date of the sequence
            dt = find_date(record.description)
            # print(type(dt))
            NoneType = type(None)
            if type(dt) != NoneType:
                dt = dt.strftime("%Y-%b")
                # print(month, dt)
                if dt == month:
                    list_for_save.append(record)
            else:
                continue

        input_handle.close()
        if len(list_for_save) == 0:
            continue

        try:
            os.mkdir(save_path + month)
        except OSError as error:
            print(error)
        num_of_seq_in_month = len(list_for_save)

        ####################################################################
        ####### if the month has less than 10 sequences we dont care ######

        if num_of_seq_in_month < 10:
            continue
        ####################################################################

        SeqIO.write(list_for_save, save_path + month + '/' + month + '_sequences.fa', 'fasta')

        pfm, count_matrix = simple_PFM(file_and_path=save_path + month + '/' + month + '_sequences.fa',
                                       ref_sequence=ref_sequence)
        mutation_df = simple_mutation_finder(count_df=count_matrix, protein='spike', ref_sequence=ref_sequence)

        pfm.to_excel(save_path + month + '/' + month + '_PFM.xlsx')
        mutation_df.to_excel(save_path + month + '/' + month + '_mutation.xlsx')

        mutation_df.reset_index(inplace=True)
        mutation_df.drop(mutation_df[mutation_df['number of occurrences'] == 0].index, inplace=True)
        mutation_df['number of occurrences'] = (mutation_df['number of occurrences'] / num_of_seq_in_month) * 100
        mutation_df['mutation aa'] = mutation_df['new_index_no'].apply(lambda x: x[-1:])
        filtered_df = mutation_df[mutation_df['number of occurrences'] > 1].copy(
            deep=True)  # droping all sequense that are less the 1 percent
        num_of_seq_in_month = str(num_of_seq_in_month)
        filtered_df['amino acid position'] = filtered_df['amino acid position'].apply(lambda x: int(x))
        filtered_df.reset_index(inplace=True, drop=True)
        the_mane_df[month + ' ' + num_of_seq_in_month] = 0
        filtered_df.to_excel(save_path + month + '/' + month + '_filtered.xlsx')
        for idx in filtered_df.index:
            a = int(filtered_df.loc[idx, 'amino acid position'])
            temp_list = list(the_mane_df['position'])
            if a in temp_list:

                # if filtered_df.at[idx,'amino acid position'] in the_mane_df['position']:
                the_mane_df.reset_index(inplace=True, drop=True)
                for mane_idx in the_mane_df.index:

                    b = int(the_mane_df.loc[mane_idx, 'position'])
                    x = the_mane_df.loc[mane_idx, 'mutated aa']
                    y = filtered_df.loc[idx, 'mutation aa']

                    if x != y and b == a:
                        if temp_list.count(a) > 1:
                            continue
                        else:
                            new_mane_idx = len(the_mane_df.index) + 1
                            the_mane_df.loc[new_mane_idx, 'mutated aa'] = filtered_df.loc[idx, 'mutation aa']
                            the_mane_df.loc[new_mane_idx, 'position'] = filtered_df.loc[idx, 'amino acid position']
                            the_mane_df.loc[new_mane_idx, month + ' ' + num_of_seq_in_month] = filtered_df.loc[
                                idx, 'number of occurrences']
                            break

                    if x == y and b == a:
                        the_mane_df.loc[mane_idx, month + ' ' + num_of_seq_in_month] = filtered_df.loc[
                            idx, 'number of occurrences']
                        break

                    # if the_mane_df.at[mane_idx, 'mutated aa'] == filtered_df.at[idx, 'mutation aa'] and the_mane_df.at[mane_idx, 'position'] == filtered_df.at[idx, 'amino acid position']:
                    #     the_mane_df.at[mane_idx, month] = filtered_df.at[idx, 'number of occurrences']
                    #     break
                    # else:
                    #     if the_mane_df.at[mane_idx, 'position'] == filtered_df.at[idx, 'amino acid position']and the_mane_df.at[mane_idx, 'mutated aa'] != filtered_df.at[idx, 'mutation aa']:
                    #         new_mane_idx = len(the_mane_df.index) + 1
                    #         the_mane_df.at[new_mane_idx, 'mutated aa'] = filtered_df.at[idx, 'mutation aa']
                    #         the_mane_df.at[new_mane_idx, 'position'] = filtered_df.at[idx, 'amino acid position']
                    #         the_mane_df.at[new_mane_idx, month] = filtered_df.at[idx, 'number of occurrences']
                    #         break


            else:
                the_mane_df.reset_index(inplace=True, drop=True)
                new_mane_idx = len(the_mane_df.index) + 1
                the_mane_df.loc[new_mane_idx, 'mutated aa'] = filtered_df.loc[idx, 'mutation aa']
                the_mane_df.loc[new_mane_idx, 'position'] = filtered_df.loc[idx, 'amino acid position']
                the_mane_df.loc[new_mane_idx, month + ' ' + num_of_seq_in_month] = filtered_df.loc[
                    idx, 'number of occurrences']

            the_mane_df.fillna(0, inplace=True)
        the_mane_df.sort_values(by=['position'], inplace=True)
    the_mane_df.to_excel(save_path + 'the_df.xlsx')

    print(
        'threading {} ended....'.format(country))  # ----------------------------------------------------------------end


# ----------------------------------------------------------------------------------------------------------------------
# new location finder 25.1.21
# ----------------------------------------------------------------------------------------------------------------------

def find_seq_location(record_discription):
    country_list = ['China', 'USA', 'Australia', 'Singapore', 'Vietnam', 'Hong Kong', 'Japan', 'Germany', 'Thailand',
                    'Malaysia', 'Taiwan', 'South Korea', 'Cambodia', 'France', 'Nepal', 'Canada', 'Italy',
                    'United Arab Emirates', 'Czech Republic', 'United Kingdom', 'Philippines', 'India', 'Viet Nam',
                    'Sweden', 'Switzerland', 'Spain', 'Saudi Arabia', 'Oman', 'Mexico', 'Netherlands', 'Georgia',
                    'Iran', 'Iceland', 'Senegal', 'Belgium', 'Ireland', 'Panama', 'Denmark', 'Austria', 'Brazil',
                    'Lithuania', 'Israel', 'Nigeria', 'Luxembourg', 'Lebanon', 'New Zealand', 'Morocco', 'Norway',
                    'Saint Barth\u00e9lemy', 'Greece', 'Slovenia', 'Jordan', 'Chile', 'Jamaica', 'Russia', 'Mongolia',
                    'Turkey', 'Costa Rica', 'Palestine', 'Tunisia', 'Mali', 'Portugal', 'Bosnia and Herzegovina',
                    'Colombia', 'Bahrain', 'Latvia', 'Peru', 'South Africa', 'Guadeloupe', 'Indonesia', 'Finland',
                    'Poland', 'Pakistan', 'Uganda', 'Hungary', 'Ecuador', 'Faroe Islands', 'Reunion', 'Kazakhstan',
                    'Uruguay', 'Slovakia', 'Zambia', 'Benin', 'Guam', 'Aruba', 'Democratic Republic of the Congo',
                    'Croatia', 'Estonia', 'Cyprus', 'Kenya', 'Qatar', 'Bulgaria', 'Moldova', 'Egypt', 'Cuba',
                    'Madagascar', 'Curacao', 'Gabon', 'Belize', 'Gambia', 'Romania', 'Argentina', 'Botswana',
                    'Guatemala', 'Bolivia', 'Crimea', 'Kuwait', 'Andorra', 'Belarus', 'Serbia', 'Sri Lanka',
                    'Saint Martin', 'Montenegro', 'Sierra Leone', 'Brunei', 'Ghana', 'Zimbabwe', 'Dominican Republic',
                    'Algeria', 'Bangladesh', 'Ukraine', 'Equatorial Guinea', 'Myanmar', 'Venezuela', 'Suriname',
                    'North Macedonia', 'Iraq', 'Republic of Congo', 'Rwanda', 'Guinea', 'Cameroon', 'Gibraltar',
                    'Papua New Guinea', 'Mozambique', 'Malta', 'Burkina Faso', 'Albania', 'El Salvador', 'St Eustatius',
                    'Bermuda', 'Mayotte']

    record_str = str(record_discription)
    record_list = record_str.split('|')
    presumed_location = record_list[-1]
    if presumed_location not in country_list:
        # print(record_discription)
        try:
            match = re.search(r'hCoV-19/\w+/', record_discription)
            location = match.group()
            location = location[8:]
            location = location[:-1]
            if location in ['England', 'Scotland', 'Wales', 'Northern Ireland']:
                return 'United Kingdom'
            else:
                return location
        except:
            try:
                match = re.search(r'hCoV-19/\w+\s\w+/', record_discription)
                location = match.group()
                location = location[8:]
                location = location[:-1]
                if location in ['England', 'Scotland', 'Wales', 'Northern Ireland']:
                    return 'United Kingdom'
                else:
                    return location
            except:
                try:
                    match = re.search(r'hCoV-19/\w+\s\w+\s\w+/', record_discription)
                    location = match.group()
                    location = location[8:]
                    location = location[:-1]
                    if location in ['England', 'Scotland', 'Wales', 'Northern Ireland']:
                        return 'United Kingdom'
                    else:
                        return location
                except:
                    try:
                        match = re.search(r'hCoV-19/\w+\s\w+\s\w+\s\w+/', record_discription)
                        location = match.group()
                        location = location[8:]
                        location = location[:-1]
                        if location in ['England', 'Scotland', 'Wales', 'Northern Ireland']:
                            return 'United Kingdom'
                        else:
                            return location
                    except:
                        return None
    else:
        return presumed_location


# ----------------------------------------------------------------------------------------------------------------------
# find the states
# ----------------------------------------------------------------------------------------------------------------------


def find_seq_states(record_f):
    try:
        match = re.search(r'\|hCoV-19\^\^\w+\||\|hCoV-19\^\^\w+\s\w+\||\|hCoV-19\^\^\w+\s\w+\s\w+\|', record_f)
        location = match.group()
        return location[10:-1]
    except:
        print(record_f)


# ----------------------------------------------------------------------------------------------------------------------
# reading in the reference  sequences and aligning them

# Do a global alignment. Identical characters are given 2 points,
# 1 point is deducted for each non-identical character. >>>
# for a in pairwise2.align.globalmx("ACCGT", "ACG", 2, -1):

# Same as above, except now 0.5 points are deducted when opening a
# gap, and 0.1 points are deducted when extending it. >>>
# for a in pairwise2.align.globalms("ACCGT", "ACG", 2, -1, -.5, -.1):
# ----------------------------------------------------------------------------------------------------------------------
def deletions(folder_path, file, ref_sequence):
    ref_sequence = ref_sequence + '*'
    print(len(ref_sequence))
    list_to_save = []
    count = 0
    start = time.time()
    seq_discarded_count = 0
    ref_cahnge = 0
    input_handle = open(folder_path + file + '.fa', "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        if len(record.seq) < (len(ref_sequence) - 15):
            seq_discarded_count += 1
            continue
        if len(record.seq) == (len(ref_sequence)):
            list_to_save.append(record)
        else:
            alignments = pairwise2.align.globalms(ref_sequence, record.seq, 2, -0.5, -3, -2)
            if alignments[0][1][-1] == '-':
                continue
            # if len(alignments[0][0]) > 1273:
            #     ref_cahnge += 1
            #     continue
            if len(alignments[0][0]) != 1274:
                seq_discarded_count += 1
                print(len(record.seq))
                print('len of ref', len(alignments[0][0]))
                print('len of new align', len(alignments[0][1]))
                continue
            n_record = SeqRecord(
                Seq(alignments[0][1]),
                id=record.id,
                name=record.name,
                description=record.description,
            )
            list_to_save.append(n_record)

        count += 1
        end = time.time()
        delta = end - start
        # print(int(delta/5))
        if int(delta % 10) == 0:
            print(" %.2f sequences read in" % count)
    input_handle.close()

    SeqIO.write(list_to_save, folder_path + file + '_with_alignment' + '.fa', 'fasta')
    print('number of seq discarded', seq_discarded_count)


# ----------------------------------------------------------------------------------------------------------------------
# clean up fasta file from unknowns
# ----------------------------------------------------------------------------------------------------------------------

def clean_x_from_fasta(file_and_path, save_path, file_name_for_save, ref_seq_f):
    """
    this function will remove all sequences with amino acid x (unknown)
    :param file_and_path: full path with folder and it must contain the file extension
    :return: new file to the path sent in the second arg
    """
    seq_mot_in_range = 0
    all_seq_count = 0
    more_than_one_x = 0
    more_than_five = 0
    list_for_save = []
    count_X_seq = 0
    seq_with_Illegal_character = 0
    range_val = len(ref_seq_f) * 0.025
    low_limit = len(ref_seq_f) - range_val
    hight_limit = len(ref_seq_f) + range_val

    input_handle = open(file_and_path, "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        all_seq_count += 1
        for aa in record.seq:
            if aa.islower():
                seq_with_Illegal_character += 1
                break
        # print(record.seq)
        if record.seq.find('X') == -1 and low_limit < len(record.seq) < hight_limit:
            list_for_save.append(record)
        else:
            if len(record.seq) < low_limit or len(record.seq) > hight_limit:
                seq_mot_in_range += 1
        if record.seq.find('X') != -1:
            count_X_seq += 1
            if record.seq.count('X') < 2:
                more_than_one_x += 1
            if 1 < record.seq.count('X') < 6:
                more_than_five += 1

    input_handle.close()

    # discarding x
    SeqIO.write(list_for_save, save_path + file_name_for_save + '.fa', 'fasta')
    print('Number of sequences in fasta file :', all_seq_count)
    print('Number of sequences that include X :', count_X_seq)
    print('Number of sequences that are not in range:', seq_mot_in_range)
    print('Number of sequences that include only X one :', more_than_one_x)
    print('Number of sequences that includemore than one x and less than five :', more_than_five)
    print('Number of sequences after discarding X', len(list_for_save))
    print('Number of sequnces with Illegal_character', seq_with_Illegal_character)


# ----------------------------------------------------------------------------------------------------------------------
# finding if a mutation is inside a peptide
# ----------------------------------------------------------------------------------------------------------------------


def is_mutation_in_peptide(mut_pos_f, pep_f, ref_sequence_f):
    """
    this function will cheak if the position is in the peptide or not
    :param mut_pos_f: position number (int)
    :param pep_f: peptide sequence (aminoaacid)
    :return: boll
    """
    mut_pos_f_int = int(mut_pos_f)
    start_pos_f = find_peptide_position(pep_f, ref_seq=ref_sequence_f)
    end_pos_f = start_pos_f + len(pep_f) - 1
    if start_pos_f <= mut_pos_f_int <= end_pos_f:
        return True
    else:
        return False


# ----------------------------------------------------------------------------------------------------------------------
# creating a mutated peptide
# ----------------------------------------------------------------------------------------------------------------------


def creat_the_mutated_peptide(peptide_f, mutation_pos_f, mutation_aa_f, ref_sequence_f):
    """
    this function will get a peptide and a mutation
    and the rturn the mutaded peptide
    :return: the peptide with the mutation
    """
    peptide_start = find_peptide_position(peptide_f, ref_seq=ref_sequence_f)
    peptide_end = peptide_start + len(peptide_f)
    # print(peptide,peptide_start)

    if mutation_pos_f == peptide_end:
        mut_idx = len(peptide_f) - 1
    else:
        mut_idx = mutation_pos_f - peptide_start
    pep_list = list(peptide_f)
    pep_list[mut_idx] = mutation_aa_f
    return ''.join(pep_list)


# ----------------------------------------------------------------------------------------------------------------------
# ccreat a file with a list of peptides
# ----------------------------------------------------------------------------------------------------------------------


def creat_pep_file_from_list(pep_list_f, file_path_and_name_f):
    """
    the net MHC predictor needs to get a peptide file as an input so this function will creat that file
    from a given list of peptides
    :param pep_list_f:
    :param file_path_and_name_f:
    """
    epfile = open(file_path_and_name_f, 'w')

    # epetope = 'QYIKWPWYI'
    for p_f in pep_list_f:
        s = p_f + "\n"

        # Writing a string to file
        epfile.write(s)

    # Closing file
    epfile.close()


# ----------------------------------------------------------------------------------------------------------------------
# creating the PFM with list of sequences
# ----------------------------------------------------------------------------------------------------------------------

def simple_list_PFM(list_of_recoreds, ref_sequence):
    ####################################################################################################################
    ###############################   need to take into  counsiduration - b z ##########################################

    mane_dict = {}
    aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
    for d in aa_list:
        mane_dict[d] = subdict = {}

    max = len(ref_sequence)

    # index_list = list(range(1, max + 1))
    # pfm_df = pd.DataFrame(0, index=index_list, columns=aa_list)

    count = 0
    start = time.time()
    for record in list_of_recoreds:

        ################## stupid asterisk at the end of the sequence ##############
        if record.seq[-1:] == '*':
            record_seq = record.seq[:-1]
        else:
            record_seq = record.seq

        if len(record_seq) != max:
            print('One of the sequences in the file is not the length of the reference length')
            break
        # print(record.seq)
        for pos in range(len(record_seq)):
            record_aa = record_seq[pos]
            index = pos + 1
            # print(mane_dict)
            if index in mane_dict[record_aa].keys():
                mane_dict[record_aa][index] = mane_dict[record_aa][index] + 1
            else:
                mane_dict[record_aa][index] = 1

        count += 1

    pfm_df = pd.DataFrame.from_dict(mane_dict)
    pfm_df.sort_index(inplace=True)
    pfm_df.fillna(0, inplace=True)
    count_matrix = pfm_df.copy(deep=True)

    pfm_df['sum'] = pfm_df.sum(axis=1)
    for col in pfm_df.columns:
        pfm_df[col] = pfm_df[col] / pfm_df['sum']
    pfm_df.drop('sum', axis=1, inplace=True)
    end = time.time()
    # delta = end - start
    # print(delta)
    # # print(int(delta/5))
    # if int(delta % 10) == 0:
    #     print(" %.2f sequences read in" % count)
    return pfm_df, count_matrix


# ----------------------------------------------------------------------------------------------------------------------
# creating the PFM
# ----------------------------------------------------------------------------------------------------------------------

def simple_PFM(file_and_path, ref_sequence):
    ####################################################################################################################
    ###############################   need to take into  counsiduration - b z ##########################################

    mane_dict = {}
    aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
    for d in aa_list:
        mane_dict[d] = subdict = {}

    max = len(ref_sequence)

    # index_list = list(range(1, max + 1))
    # pfm_df = pd.DataFrame(0, index=index_list, columns=aa_list)

    count = 0
    start = time.time()
    input_handle = open(file_and_path, "r")
    for record in SeqIO.parse(input_handle, "fasta"):

        ################## stupid asterisk at the end of the sequence ##############
        if record.seq[-1:] == '*':
            record_seq = record.seq[:-1]
        else:
            record_seq = record.seq

        if len(record_seq) != max:
            print('One of the sequences in the file is not the length of the reference length')
            break
        # print(record.seq)
        for pos in range(len(record_seq)):
            record_aa = record_seq[pos]
            index = pos + 1
            # print(mane_dict)
            if index in mane_dict[record_aa].keys():
                mane_dict[record_aa][index] = mane_dict[record_aa][index] + 1
            else:
                mane_dict[record_aa][index] = 1

        count += 1

    input_handle.close()
    pfm_df = pd.DataFrame.from_dict(mane_dict)
    pfm_df.sort_index(inplace=True)
    pfm_df.fillna(0, inplace=True)
    count_matrix = pfm_df.copy(deep=True)

    pfm_df['sum'] = pfm_df.sum(axis=1)
    for col in pfm_df.columns:
        pfm_df[col] = pfm_df[col] / pfm_df['sum']
    pfm_df.drop('sum', axis=1, inplace=True)
    end = time.time()
    # delta = end - start
    # print(delta)
    # # print(int(delta/5))
    # if int(delta % 10) == 0:
    #     print(" %.2f sequences read in" % count)
    return pfm_df, count_matrix


# ----------------------------------------------------------------------------------------------------------------------
# simple entropy finder
# ----------------------------------------------------------------------------------------------------------------------


def simple_entropy_function(path_to_save, save_file_name, pfm_df):
    """
    this function takes in a pfm data frame and the calculates for every position its entropy
    then it returns the result as a data frame that has one column for the position and a
    corresponding one for the entropy score

    """
    dict_entropy = {}
    df_transposed = pfm_df.transpose()
    for col in list(df_transposed):
        dict_entropy[col] = entropy(df_transposed[col], base=2)

    df = pd.DataFrame(list(dict_entropy.items()))
    df.rename({1: 'entropy'}, axis=1, inplace=True)

    df.set_index(0, inplace=True)
    df.to_excel(path_to_save + save_file_name + '.xlsx')
    return df


# ----------------------------------------------------------------------------------------------------------------------
# simple mutation finder
# ----------------------------------------------------------------------------------------------------------------------


def simple_mutation_finder(count_df, protein, ref_sequence):
    # df_with_f = 'C:/Users/User/Dropbox/HLA/only_spike/Spike/Spike_COUNT_matrix_with_repetitions.xlsx'
    p = protein
    df_with = count_df.copy(deep=True)
    # dropping wild type from df
    for position in range(len(ref_sequence)):
        df_with.loc[(position + 1), ref_sequence[position]] = 0
    # find wild type for df_with

    stacked = df_with.stack().to_frame()
    stacked.reset_index(inplace=True)
    stacked['level_0'] = stacked['level_0'].apply(lambda x: str(x))
    stacked['new_index_no'] = p + '_' + stacked['level_0'] + '_' + stacked['level_1']
    # stacked.drop(['level_0'], axis=1, inplace=True)
    stacked.drop(['level_1'], axis=1, inplace=True)
    stacked = stacked.set_index('new_index_no')
    # change col 0 to 1
    stacked.rename({0: 'number of occurrences'}, axis=1, inplace=True)
    stacked.rename({'level_0': 'amino acid position'}, axis=1, inplace=True)

    # add the wild type amino acid for every mutation
    stacked['wild_type_aa'] = 0
    for i in list(stacked.index):
        wild_type_position = int(stacked.loc[i, 'amino acid position'])
        stacked.loc[i, 'wild_type_aa'] = ref_sequence[wild_type_position - 1]
    return stacked


# ----------------------------------------------------------------------------------------------------------------------
# peptide position finder
# ----------------------------------------------------------------------------------------------------------------------


def find_peptide_position(peptide, ref_seq):
    """
    this function will find the place of the peptide in the protein(ref_seq)
    :param seq: peptide
    :ref_seq: protein sequence
    :return: start position
    """
    if ref_seq.find(peptide) != -1:
        return (ref_seq.find(peptide) + 1)
    else:
        if "X" in peptide:
            return -1
        else:
            print('epitope {} not found'.format(peptide))
            return -1


# ----------------------------------------------------------------------------------------------------------------------
# paired_mutations finder
# ----------------------------------------------------------------------------------------------------------------------


def paired_mutations(folder_path, file_name, list_of_mutation):
    """
    this function take a list of mutations and the fasta file that contains all the sequences of interest
    in the following format aa _ position _ mut aa

    :param folder_path: path of folder were file is at
    :param file_name: the name of fasta file with proper ending
    :param list_of_mutation: list of all mutations in the following format wt
    :return: a xlsx file with the id counts
    """
    df = pd.DataFrame(columns=['mutation', 'record ID'])

    idx = 0
    for mut in list_of_mutation:
        mut_list = mut.split('_')
        position = int(mut_list[1])
        wt_aa = mut_list[0]
        mut_aa = mut_list[2]
        input_handle = open(folder_path + file_name, "r")
        count = 0
        for record in SeqIO.parse(input_handle, "fasta"):
            count += 1
            aa = record.seq[position - 1]
            if aa == mut_aa:
                df.loc[idx, 'mutation'] = mut
                df.loc[idx, 'record ID'] = record.description
                idx += 1

    df = df['record ID'].value_counts()

    df.to_excel(folder_path + 'moutation_groupby.xlsx')


# ----------------------------------------------------------------------------------------------------------------------
# for every protein the most common length is calculated
# ----------------------------------------------------------------------------------------------------------------------

def count_df_for_protein(p_file, p_file_no_red):
    """
    this function will find:
    1. the len of every sequence
    2. find the number of repetition for every sequence
    then creat a df that contains all the data for protein p
    :param p_file: fasta file for protein p
    :param p_file_no_red: fasta file for protein p with no redundancies
    :return: df for further analysis
    """

    list_of_id = []
    list_of_seq = []
    input_handle = open(p_file_no_red, "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        list_of_id.append(record.id)
        list_of_seq.append(record.seq)
    input_handle.close()
    data = {'ID': list_of_id, 'seq': list_of_seq}
    df = pd.DataFrame(data)
    df['seq_len'] = df['seq'].apply(lambda x: len(x))
    df['repetitions'] = df['seq'].apply(lambda x: 0)

    # finding the number of repetition
    df.set_index('seq', inplace=True)
    input_handle = open(p_file, "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        df.loc[record.seq, 'repetitions'] += 1
    input_handle.close()
    df.reset_index(inplace=True)
    return df


# ----------------------------------------------------------------------------------------------------------------------
# for every protein a fasta file is made for furthered analysis
# ----------------------------------------------------------------------------------------------------------------------

def p_to_fasta(p, no_x_file_path, p_folder):
    """
    this function takes in a big fasta file that has different proteins and saves the selected ones in
    the same file
    :param p: the protein name
    :param no_x_file_path: the path and file name of the file that has no sequences with x in it
    :param p_folder: the path and the folder that was created for said protein
    :return: Number of sequences in file and it creates the fasta file with the selected protein
    """

    list_for_saved_protein = []
    input_handle = open(no_x_file_path, "r")
    count_file_seq = 0

    for record in SeqIO.parse(input_handle, "fasta"):
        if record.id.startswith(p + '|'):
            list_for_saved_protein.append(record)
            count_file_seq += 1

    input_handle.close()
    SeqIO.write(list_for_saved_protein, p_folder + p + '.fa', 'fasta')
    # print('Number of sequences in file - ' + p + ' is :', count_file_seq)
    return count_file_seq


# ----------------------------------------------------------------------------------------------------------------------
# for every protein fasta file a new fasta file is created without redundancies
# ----------------------------------------------------------------------------------------------------------------------

def fasta_no_redundanc(p, p_folder):
    """
    this function takes a fasta file and reads it in , finds the redundancies and then creats a new
    fasta file with no redundancies
    :param p: the protein name
    :param p_folder: the folder path and file name of said protein
    :return: saved file in p_folder directory that has no redundancies and there count
    """
    number_of_seq = 0
    no_redundancies_list = []
    no_redundancies_list_to_save = []
    count_file_seq = 0
    input_handle = open(p_folder + p + '.fa', "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        number_of_seq += 1
        if record.seq not in no_redundancies_list:
            no_redundancies_list.append(record.seq)
            no_redundancies_list_to_save.append(record)
            count_file_seq += 1
    input_handle.close()
    print('The number of sequences in protein is : ', number_of_seq)
    print('The number of sequences in protein ' + p + ' that are not redundant is :', count_file_seq)
    SeqIO.write(no_redundancies_list_to_save, p_folder + p + '_no_redundancies.fa', 'fasta')
    return count_file_seq


# ----------------------------------------------------------------------------------------------------------------------
# for every protein a PFM matrix is calculated with or without including repetitions
# ----------------------------------------------------------------------------------------------------------------------


def find_pfm(p, p_folder, p_df, repetitions=False, msa=False, BandZ=False):
    """
    finding the pfm is possible when doing the fallowing steps:
    1 - finding the max length of the protien sequence
    2 - creating a empty matrix with the length a the aa positions (roes)
        and all known amino acids as the columns
    3 - iterating over every sequence in the data frame and assigning its amino acid to the
        corresponding place in the matrix, the row is the position in the sequence.
    4 - if repetitions is true so the assigment will be multiplied by it.
    5 - the sum of evey row is calculated and then every row is divided by it.

    :param p: protein name
    :param p_folder: protein folder location
    :param p_df: protein df that i created in function "count_df_for_protein"
    :return: matrix df as a parameter and a saved file
        if count true then function returns only count matrix
        if count false , function returns frequency matrix
    """
    if msa:
        if BandZ:
            aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                       'V', '-', 'B', 'Z']
        else:
            aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                       'V', '-']
    else:
        if BandZ:
            aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                       'V', 'B', 'Z']
        else:
            aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                       'V']
    max = p_df['seq_len'].max()
    index_list = list(range(1, max + 1))
    pfm_df = pd.DataFrame(0, index=index_list, columns=aa_list)
    for i in p_df.index:
        current_seq = p_df.loc[i, 'seq']
        for j in range(len(current_seq)):
            aa = current_seq[j]
            if repetitions:
                pfm_df.loc[j + 1, aa] = pfm_df.loc[j + 1, aa] + (1 * p_df.loc[i, 'repetitions'])
            else:
                pfm_df.loc[j + 1, aa] = pfm_df.loc[j + 1, aa] + 1
    count = pfm_df.copy(deep=True)
    pfm_df['sum'] = pfm_df.sum(axis=1)
    for col in pfm_df.columns:
        pfm_df[col] = pfm_df[col] / pfm_df['sum']
    pfm_df.drop('sum', axis=1, inplace=True)
    if repetitions:
        pfm_df.to_excel(p_folder + p + '_PFM_with_repetitions.xlsx')
        count.to_excel(p_folder + p + '_COUNT_matrix_with_repetitions.xlsx')
        return pfm_df, count
    else:
        pfm_df.to_excel(p_folder + p + '_PFM.xlsx')
        count.to_excel(p_folder + p + '_COUNT_matrix.xlsx')
        return pfm_df, count


# ----------------------------------------------------------------------------------------------------------------------
# The function for finding the entropy for evey protein
# ----------------------------------------------------------------------------------------------------------------------


def entropy_function(p, p_folder, pfm_df_with, pfm_df):
    """
    this function takes in a pfm data frame and the calculates for every position its entropy
    then it returns the result as a data frame that has one column for the position and a
    corresponding one for the entropy score
    :param p: protein name
    :param p_folder: protein folder location
    :param p_df: protein df that i created in function "count_df_for_protein"
    # :return: graph
    """
    dict_entropy = {}
    df_transposed = pfm_df.transpose()
    for col in list(df_transposed):
        dict_entropy[col] = entropy(df_transposed[col], base=2)

    dict_entropy_with = {}
    df_transposed_with = pfm_df_with.transpose()
    for col in list(df_transposed_with):
        dict_entropy_with[col] = entropy(df_transposed_with[col], base=2)

    df = pd.DataFrame(list(dict_entropy.items()))
    df.rename({1: 'no'}, axis=1, inplace=True)
    df_with = pd.DataFrame(list(dict_entropy_with.items()))
    df_with.rename({1: 'with'}, axis=1, inplace=True)
    # result = pd.merge([df, df_with])
    entropy_df = df_with.merge(df, on=0, )
    entropy_df.set_index(0, inplace=True)
    entropy_df.to_excel(p_folder + p + '_entropy_df.xlsx')
    # xlim= len(entropy_df)
    # plt.subplots(figsize=(xlim/3,20))
    # plt.bar(entropy_df[0], entropy_df['no'], label='no repetitions ')
    # plt.bar(entropy_df[0]+ 0.25, entropy_df['with'], label='with repetitions', )
    # plt.xticks(entropy_df[0])
    # plt.xlabel('aa position', fontweight='bold')
    # plt.ylabel('score', fontweight='bold')
    # # plt.xlim(0,xlim)
    # plt.legend()
    # # plt.show()
    # plt.savefig(p_folder+'entropy.jpg')
    # plt.close()


# ----------------------------------------------------------------------------------------------------------------------
# mutation mutation_find
# ----------------------------------------------------------------------------------------------------------------------


def mutation_finder(df_f, df_with_f, p, p_folder, ref):
    """
This function will take the count matrix and turn it into a pandas series, the series index will contain the protein
name_position in sequence amino acid, while the value will be the number of occurrences.
This step is done for both the count matrix that includes the repetition and the one that does not.
For presentation purposes, the number of occurrences of the “wild type” amino acid is dropped according to the reference
protein.
    :param df: count df
    :param df_with: count df with repetitions
    :param p: the protein
    :param p_folder: folder file
    :return: result file, a data frame with all the mutations per position in the protein
    """

    df = df_f.copy(deep=True)
    df_with = df_with_f.copy(deep=True)
    sr_wild_type = pd.Series(ref)
    # dropping wild type from df
    for position in range(len(ref)):
        df.loc[(position + 1), ref[position]] = None
        df_with.loc[(position + 1), ref[position]] = None
    # find wild type for df_with

    stacked = df.stack().to_frame()
    stacked.reset_index(inplace=True)
    stacked['level_0'] = stacked['level_0'].apply(lambda x: str(x))
    stacked['new_index_no'] = p + '_' + stacked['level_0'] + '_' + stacked['level_1']
    # stacked.drop(['level_0'], axis=1, inplace=True)
    stacked.drop(['level_1'], axis=1, inplace=True)
    stacked = stacked.set_index('new_index_no')
    # change col 0 to 1
    stacked.rename({0: 'no'}, axis=1, inplace=True)
    stacked.rename({'level_0': 'amino acid position'}, axis=1, inplace=True)

    stacked2 = df_with.stack().to_frame()
    stacked2.reset_index(inplace=True)
    stacked2['level_0'] = stacked2['level_0'].apply(lambda x: str(x))
    stacked2['new_index_with'] = p + '_' + stacked2['level_0'] + '_' + stacked2['level_1']
    stacked2.drop(['level_0'], axis=1, inplace=True)
    stacked2.drop(['level_1'], axis=1, inplace=True)
    stacked2 = stacked2.set_index('new_index_with')
    # change col 0 to 2
    stacked2.rename({0: 'with'}, axis=1, inplace=True)
    result = pd.concat([stacked, stacked2], axis=1, sort=False)

    # add the wild type amino acid for every mutation
    result['wild_type_aa'] = 0
    for i in list(result.index):
        wild_type_position = int(result.loc[i, 'amino acid position'])
        result.loc[i, 'wild_type_aa'] = ref[wild_type_position - 1]

    result.to_csv(p_folder + p + '_mutations.csv', sep='\t')


# ----------------------------------------------------------------------------------------------------------------------
# find number of seq in fasta
# ----------------------------------------------------------------------------------------------------------------------

def find_nuber_of_seq_in_fasta(p, fasta_file):
    """
    this function is for a sanity check, it will find all the sequences for the said protien
    and then count them
    :param p: protein
    :param fasta_file: up to date fasta file from gisaid
    :return: number of sequences
    """
    count_p_seq = 0
    input_handle = open(fasta_file + '.fasta', "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        # print(record.seq)
        if record.id.startswith(p + '|'):
            count_p_seq += 1
    input_handle.close()
    return count_p_seq


# ----------------------------------------------------------------------------------------------------------------------
# dropping all seq that are not the common length
# ----------------------------------------------------------------------------------------------------------------------

def clean_the_non_common(file_and_path, common_length, save_file_name_and_path):
    """
    this function will take the fasta file of said protein including the repetitions, and then
    it will clear all sequences that are not the common length
    :param p: protein
    :param p_folder: protein folder
    :param common_length: protein most common seq length
    :return: saved fasta file with protiens that are only the common length
    """
    list_for_saved_protein = []
    input_handle = open(file_and_path, "r")
    count_file_seq = 0

    for record in SeqIO.parse(input_handle, "fasta"):
        if len(record.seq) == common_length:
            list_for_saved_protein.append(record)
            count_file_seq += 1

    input_handle.close()
    SeqIO.write(list_for_saved_protein, save_file_name_and_path, 'fasta')
    print('number of seq that were the common length: ', count_file_seq)


# ----------------------------------------------------------------------------------------------------------------------
# mutations over time
# ----------------------------------------------------------------------------------------------------------------------


def find_date(recordid):
    match = re.search(r'\d{4}-\d{2}-\d{2}|\d{4}-\d{2}-\d{1}|\d{4}-\d{1}-\d{2}|\d{4}-\d{1}-\d{1}', recordid)
    try:
        date_match = match.group()
    except:
        print('failed to retrieve date from:')
        print(recordid)
        exit()
    # checking if date is lagit
    date_list = date_match.split('-')
    if int(date_list[2]) == 0:
        date_match = date_match[:-2] + '01'
    if int(date_list[1]) == 0:
        return None
    if int(date_list[0]) != 2020 and int(date_list[0]) != 2019 and int(date_list[0]) != 2021 and int(
            date_list[0]) != 2022:
        return None
    try:
        dt = datetime.datetime.strptime(date_match.strip(), "%Y-%m-%d")
        return dt
    except:
        print('no datetime object was created')


# def find_location(recordid):
#     try:
#         match = re.search(r'hCoV-19/\w+/', recordid)
#         location = match.group()
#         location = location[8:]
#         location = location[:-1]
#         return location
#     except:
#         return None


def creat_df(mutation, folder_path, p):
    print('thread for ', mutation, 'has started')
    df = pd.DataFrame(columns=['position', 'amino acid', 'date', 'location', 'ID'])
    # print(df)
    idx = -1

    position = int(mutation.split('_')[1])
    # input_handle = open(folder_path + p + '/' + p +'_only_common_length.fa', "r")
    input_handle = open(folder_path + p + '/' + p + '_only_common_length.fa',
                        "r")  # ------------------------------reading the fasta file
    for record in SeqIO.parse(input_handle, "fasta"):
        idx += 1
        df.loc[idx, 'position'] = position
        df.loc[idx, 'amino acid'] = record.seq[position - 1]
        df.loc[idx, 'date'] = find_date(record.id)
        df.loc[idx, 'location'] = find_location(record.id)
        df.loc[idx, 'ID'] = record.id

    print('saving _mutations_location_and_date: ', mutation)
    df['date'] = pd.to_datetime(df['date'])
    df.to_csv(folder_path + p + '/' + mutation + '_location_and_date.csv', sep='\t')
    print('saved ', mutation + '_location_and_date.csv')


def monthly_data(mutation, folder_path, p):
    print('thread for ', mutation, 'monthly has started')

    monthly_resampled_data = pd.read_csv(folder_path + p + '/' + mutation + '_location_and_date.csv', sep='\t')
    monthly_resampled_data.set_index(['Unnamed: 0'], inplace=True)
    monthly_resampled_data = monthly_resampled_data.dropna(axis=0)
    monthly_resampled_data['date'] = pd.to_datetime(monthly_resampled_data['date'])
    # print(monthly_resampled_data.info())
    # monthly_resampled_data = df.set_index('date')
    # monthly_sum = pd.DatetimeIndex(monthly_resampled_data.date)
    # monthly_resampled_data = monthly_resampled_data.resample('M', on='date').count()
    monthly_sum = (
        monthly_resampled_data.groupby(['position', pd.Grouper(key='date', freq='M')]).agg({'position': 'count'}))
    monthly_sum.rename({'position': 'position count'}, axis=1, inplace=True)
    monthly_sum = monthly_sum.reset_index()

    # grouping the amino acids by month

    new_df = (monthly_resampled_data.groupby(['position', 'amino acid', pd.Grouper(key='date', freq='M')]).agg(
        {'amino acid': 'count'}))
    new_df.rename({'amino acid': 'amino acid count'}, axis=1, inplace=True)
    new_df = new_df.reset_index()
    # print(monthly_resampled_data.info())

    data = new_df.merge(monthly_sum, on=['position', 'date'])

    data['amino acid monthly frequency'] = data['amino acid count'] / data['position count']

    data.to_excel(folder_path + p + '/' + mutation + '_month_mutations.xlsx')
    print('saved ', mutation + '_month_mutations.xlsx')


def monthly_location(mutation, folder_path, p):
    # tou can run it once and it is enough
    monthly_resampled_data = pd.read_csv(folder_path + p + '/' + mutation + '_location_and_date.csv', sep='\t')
    monthly_resampled_data.set_index(['Unnamed: 0'], inplace=True)
    monthly_resampled_data = monthly_resampled_data.dropna(axis=0)
    monthly_resampled_data['date'] = pd.to_datetime(monthly_resampled_data['date'])
    monthly_sum = (monthly_resampled_data.groupby(['position', 'amino acid', pd.Grouper(key='date', freq='M')]).agg(
        {'amino acid': 'count'}))
    monthly_sum.rename({'amino acid': 'amino acid count'}, axis=1, inplace=True)
    monthly_sum = monthly_sum.reset_index()

    new_df = (monthly_resampled_data.groupby(['amino acid', 'location', pd.Grouper(key='date', freq='M')]).agg(
        {'location': 'count'}))
    new_df.rename({'location': 'location count'}, axis=1, inplace=True)
    new_df = new_df.reset_index()
    # data = new_df.merge(monthly_sum, on=['position','date'])

    new_df.to_excel(folder_path + p + '/' + mutation + '_month_location.xlsx')
    print('saved ', mutation + '_month_location.xlsx')


if __name__ == '__main__':
    print('Tool File used...')
