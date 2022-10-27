"""
##############################################################################################
###################################### THE GRAND FINALLY T CELL SCORE ########################
##############################################################################################

this scrip will give us the top peptides that are suspected to have mutation that escape hla binding,
or TCR recognition.

parameters fot the final ranking df:

1. IEDB suorce
2. seasonal corona similarity
3. similarity to self (blast)
4. is in ancher resedu
5. number of hlas that bind
6. has been seen for mor than 3 months
7. has been seen for mor than 11 months
8. mean entropy for the peptide

## filer all peptides that had no change

files needed:

1. sras IEDB data
2. HKU1 and OC43 full genome sequences
3. normolised by country df for all protiens
4. entropy for all the protiens
5. human proteome


by: Sinai Sacharen

"""

# ---------------------""" imports  & paths"""------------------------- <editor-fold>
from Tool_File import *
from utiles import *

from scipy.stats import fisher_exact
import statsmodels

starttime = time.time()

path_with_data = '/home/sacharen/Documents/omicron/'
path_to_save = '/home/sacharen/Documents/omicron/'

if os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()

# ------------------------------------------------------------------------ </editor-fold>

# --------------------------------"""files to import"""---------------- <editor-fold>
allele_data = pd.read_excel(path_with_data + "new_allele_data.xlsx", index_col=0).reset_index()
ref_df.set_index(['GISAID name'], inplace=True)
monthly_df = pd.read_csv(path_with_data + 'monthly_df.csv', sep=',')

sup_list = ['HLA-A01:01', 'HLA-A02:01', 'HLA-A03:01', 'HLA-A24:02', 'HLA-A26:01', 'HLA-B07:02', 'HLA-B08:01',
            'HLA-B27:05', 'HLA-B39:01', 'HLA-B40:01', 'HLA-B58:01', 'HLA-B15:01']

# ------------------------------------------------------------------------ </editor-fold>

# --------------------------------""" Functions """----------------------- <editor-fold>

def calc_pepe_mean_entropy(peptide_str, pf):
    """
    here we get a peptide and the spikes entropy per position as a global parameter
    so the mean entropy of a peptide can be calculated
    :param peptide_str: a peptide
    :return: the mean peptide Entropy score
    """
    global ref_df
    prot_seq = ref_df.loc[pf, 'sequence']
    entropy_df = pd.read_excel(path_with_data + 'entropy/' + pf + '_entropy_matrix.xlsx')
    peptide_start = find_peptide_position(peptide_str, ref_seq=prot_seq)
    peptide_end = peptide_start + len(peptide_str)
    # print(peptide,peptide_start)
    pep_entropy = entropy_df.loc[peptide_start:peptide_end, 'entropy'].mean()
    return pep_entropy


def find_number_of_monthly_occurrences(position, mut_aa, protein_f, freq=0.0, mean_months=False, max_month=False):
    """
    this function takes in a protein a position and an amino aced
    the function will return the number of months the query was grater than zero
    :param mean_months:
    :param freq:
    :param position:
    :param mut_aa:
    :param protein_f:
    :return: number of months
    """
    global monthly_df_melt

    # slice the df
    sliced_df = monthly_df_melt[
        (monthly_df_melt['Unnamed: 0'] == position) & (monthly_df_melt['variable'] == mut_aa) & (
                monthly_df_melt['protein'] == protein_f) & (monthly_df_melt['value'] > freq)]
    number_of_months = list(sliced_df['value'])
    if max_month == True:
        idx = sliced_df['value'].idxmax()
        max_mont = sliced_df.loc[idx, 'date']
        print(max_mont)
        return max_mont
    if mean_months == True:
        month_mean = np.mean(number_of_months)
        return month_mean
    else:
        return len(number_of_months)


# ------------------------------------------------------------------------ </editor-fold>

# ---------------------MAIN----------------------------------------------- <editor-fold>

# ------------concat all pritien dfs together to get a main df------------ <editor-fold>

text_files = [f for f in os.listdir(path_with_data+'/omicronHla/') if f.endswith('final_df.xlsx')]
all_prot_df = pd.DataFrame()
count = 0
for file in text_files:

    temp_df = pd.read_excel(path_with_data+'/omicronHla/' + file, index_col=0)
    temp_df['protein'] = temp_df['Peptide'].apply(lambda x: file[:-14])
    print(file[:-14])
    # shift column 'Name' to first position
    first_column = temp_df.pop('protein')
    temp_df.insert(0, 'protein', first_column)
    if count == 0:
        all_prot_df = temp_df.copy()
        count += 1
        continue
    all_prot_df = pd.concat([all_prot_df, temp_df])

all_prot_df['Delta'] = all_prot_df['MUT mean'] - all_prot_df['mean']
all_prot_df.to_csv(path_with_data + 'main_df.csv')
main_df = all_prot_df
# ------------------------------------------------------------------------ </editor-fold>

### add all the IEDB DATA

# sars t cell
iedb_t_cell = pd.read_csv(path_with_data + 'sars_t_cell.csv', skiprows=1)
t_cell_set = set(iedb_t_cell['Description'])
main_df['IEDB T cell'] = main_df['Peptide'].apply(lambda x: True if x in t_cell_set else False)
# sars hla
sars_hla = pd.read_csv(path_with_data + 'sars_hla.csv', skiprows=1)
sars_hla = set(sars_hla['Description'])
main_df['IEDB SARS HLA'] = main_df['Peptide'].apply(lambda x: True if x in sars_hla else False)

# allele dta to dict
allele_dict = {}
for idx in allele_data.index:
    key = allele_data.at[idx, 'Allele']
    allele_dict[key] = allele_data.at[idx, 'Frequency']

wt_hla_list = []
mut_hla_list = []
allele_freq_df = main_df.copy()
for col in main_df.columns:
    if col.startswith('HLA-'):
        allele_freq_df[col] = main_df[col] * allele_dict[col]
        wt_hla_list.append(col)
    if col.startswith('MUT HLA-'):
        allele_freq_df[col] = main_df[col] * allele_dict[col[4:]]
        mut_hla_list.append(col)

main_df['Allele freq mean'] = allele_freq_df[wt_hla_list].mean(axis=1)
main_df['Allele freq MUT mean'] = allele_freq_df[mut_hla_list].mean(axis=1)
main_df['Allele freq delta'] = main_df['Allele freq MUT mean'] - main_df['Allele freq mean']

### number of hlas that had a score grater than 2
main_df['Nuber of WT Binders'] = main_df[wt_hla_list].apply(count_number_of_binders, axis=1)
main_df['Nuber of MUT Binders'] = main_df[mut_hla_list].apply(count_number_of_binders, axis=1)

### adding peptide mean entropy no bias
main_df['Peptide mean entropy'] = main_df.apply(lambda x: calc_pepe_mean_entropy(x['Peptide'], x['protein']), axis=1)

### adding monthly data
monthly_df_melt = pd.melt(monthly_df, id_vars=['Unnamed: 0', 'date', 'protein'],
                          value_vars=['A', 'R', 'N', 'D', 'C', 'Q',
                                      'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-'])

main_df['Number of months not ZERO freq'] = main_df.apply(lambda x: find_number_of_monthly_occurrences(
    x['mutations position'], x['mutations aa'], x['protein']), axis=1)

main_df['Number of months gt 0.01'] = main_df.apply(lambda x: find_number_of_monthly_occurrences(
    x['mutations position'], x['mutations aa'], x['protein'], freq=0.01), axis=1)

main_df['Number of months gt 0.001'] = main_df.apply(lambda x: find_number_of_monthly_occurrences(
    x['mutations position'], x['mutations aa'], x['protein'], freq=0.001), axis=1)

main_df['Mutation monthly AVG'] = main_df.apply(lambda x: find_number_of_monthly_occurrences(
    x['mutations position'], x['mutations aa'], x['protein'], mean_months=True), axis=1)

main_df['Delta binders'] = main_df['Nuber of MUT Binders'] - main_df['Nuber of WT Binders']

main_df = main_df[main_df['start pep'] != -1]
main_df = main_df[main_df['mutations aa'] != '-']


# ------------------------------ anchor residues ------------------------- <editor-fold>


main_df['anchor state'] = main_df[['mutations position', 'Peptide', 'start pep', 'end pep']].apply(row_anchor_paprams,
                                                                                                   axis=1)

main_df['Hydropathicity Delta'] = main_df[['protein', 'mutations position', 'mutations aa']].apply(phisical_change,
                                                                                                   phis='H', axis=1)
main_df['Bulkiness Delta'] = main_df[['protein', 'mutations position', 'mutations aa']].apply(phisical_change, phis='B',
                                                                                              axis=1)

##### get max month #######
main_df['max month for mutation'] = main_df.apply(lambda x: find_number_of_monthly_occurrences(x['mutations position'], x['mutations aa'], x['protein'],max_month=True), axis=1)

# allele dta to dict
hla_rep_dict = {}
allele_dict = {}
sup_hla_re = {}
for idx in allele_data.index:
    key = allele_data.at[idx, 'Allele']
    allele_dict[key] = allele_data.at[idx, 'Frequency']
    hla_rep_dict[key] = [0, 0]
    sup_hla_re[key] = [0, 0]


def smart_delta(score_row):
    new_delta = 0
    new_exposure = 0
    hla_list = []
    sb_new_delta = 0
    sb_new_exposure = 0
    sb_hla_list = []
    hla_exposure_list = []
    sb_hla_exposure_list = []
    tapasin_delta = 0
    tapasin_hla_list = []
    my_new_exposure = 0
    my_new_delta = 0
    my_hla_list = []
    my_hla_exposure_list = []
    global allele_dict
    for allele in allele_dict.keys():

        if score_row[allele] > 2 and score_row['MUT ' + allele] < 2:  # find if an escape happend
            new_delta = new_delta + allele_dict[allele]
            hla_list.append(allele)
            hla_rep_dict[allele][0] += 1


        if score_row[allele] < 2 and score_row['MUT ' + allele] > 2:  # find if there was exposure
            new_exposure = new_exposure + allele_dict[allele]
            hla_exposure_list.append(allele)
            hla_rep_dict[allele][1] += 1

        if score_row[allele] > 9 and score_row['MUT ' + allele] < 9:  # find if an escape happend
            sb_new_delta = sb_new_delta + allele_dict[allele]
            sb_hla_list.append(allele)


        if score_row[allele] < 9 and score_row['MUT ' + allele] > 9:  # find if there was exposure
            sb_new_exposure = sb_new_exposure + allele_dict[allele]
            sb_hla_exposure_list.append(allele)


        if score_row[allele] > 2 and score_row['MUT ' + allele] < 2:  # find if an escape happend
            my_new_delta = my_new_delta + ((score_row[allele] - score_row['MUT ' + allele])) * allele_dict[
                allele]
            my_hla_list.append(allele)

        if score_row[allele] < 2 and score_row['MUT ' + allele] > 2:  # find if there was exposure
            my_new_exposure = my_new_exposure + (score_row['MUT ' + allele] - score_row[allele]) * \
                              allele_dict[allele]
            my_hla_exposure_list.append(allele)



    score_row['smart delta escape'] = new_delta if len(hla_list) != 0 else None
    score_row['Allele seen escape'] = hla_list
    score_row['smart delta exposure'] = new_exposure if len(hla_exposure_list) != 0 else None
    score_row['Alleles seen exposure'] = hla_exposure_list
    score_row['SB smart delta escape'] = sb_new_delta if len(sb_hla_list) != 0 else None
    score_row['SB Allele seen escape'] = sb_hla_list
    score_row['SB smart delta exposure'] = sb_new_exposure if len(sb_hla_exposure_list) != 0 else None
    score_row['SB Alleles seen exposure'] = sb_hla_exposure_list
    score_row['HLA with Tapasin score'] = tapasin_hla_list

    score_row['my smart delta escape'] = my_new_delta if len(my_hla_list) != 0 else None
    score_row['my Allele seen escape'] = my_hla_list
    score_row['my smart delta exposure'] = my_new_exposure if len(my_hla_exposure_list) != 0 else None
    score_row['my Alleles seen exposure'] = my_hla_exposure_list
    return score_row


# -----------------------------df per allele------------------------------ <editor-fold>


# ------------------------------------------------------------------------ </editor-fold>

score_df = main_df[main_df['mutations aa'] != 'C']



# score_df = score_df[score_df['Mutation monthly AVG']>0.01]
###### P& fix #######
score_df['end-mut'] = score_df['end pep'] - score_df['mutations position']
# score_df['anchor state'] = score_df['end-mut'].apply(lambda x: 'P7' if x==2)

score_df = score_df.apply(smart_delta, axis=1)
score_df.to_csv(path_to_save + 'END_GAME_omicron.csv')
