# ---------------------""" imports  & paths"""------------------------- <editor-fold>
from Tool_File import *
from utiles import *

from scipy.stats import fisher_exact
import statsmodels

path_with_data = 'C:/Users/User/Documents/COVID-19/Omicron/'
path_to_save = path_with_data

if os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()

# ------------------------------------------------------------------------ </editor-fold>

# --------------------------------"""files to import"""---------------- <editor-fold>
# score_df = pd.read_csv(path_with_data+'END_GAME_omicron.csv',sep=',',index_col=0)
score_df = pd.read_csv(path_with_data+'new final.csv',sep=',',index_col=0)

allele_data_omst = pd.read_csv(path_with_data + "Luo_SD_1_inferred_HLAallele_summary.txt", index_col=0, sep='\s+',
                               header=1)
allele_data_omst = allele_data_omst.filter(regex='^[A,B]/*', axis=0)
# allele_data_omst = allele_data_omst.filter(regex='G$',axis=0)
allele_data_omst = allele_data_omst[allele_data_omst['Freq'] > 0.01]
allele_data_omst.reset_index(inplace=True)


def get_allele_g_type(allel_str):
    # a = 'HLA-A02:02'
    allel_str = allel_str[:7]
    allel_str = allel_str.replace('*', '')
    return 'HLA-' + allel_str

allele_data_omst['Allele'] = allele_data_omst['Allele'].apply(get_allele_g_type)

# allele dta to dict
allele_high_dict = {}
for idx in allele_data_omst.index:
    key = allele_data_omst.at[idx, 'Allele']
    allele_high_dict[key] = allele_data_omst.at[idx, 'Freq']
# ------------------------------------------------------------------------ </editor-fold>

def smart_delta(score_row):

    my_new_exposure = 0
    my_new_delta = 0
    my_hla_list = []
    my_hla_exposure_list = []
    global allele_high_dict
    for allele in allele_high_dict.keys():

        if score_row[allele] > 2 and score_row['MUT ' + allele] < 2:  # find if an escape happend
            my_new_delta = my_new_delta + ((score_row[allele] - score_row['MUT ' + allele])) * allele_high_dict[
                allele]
            my_hla_list.append(allele)

        if score_row[allele] < 2 and score_row['MUT ' + allele] > 2:  # find if there was exposure
            my_new_exposure = my_new_exposure + (score_row['MUT ' + allele] - score_row[allele]) * \
                              allele_high_dict[allele]
            my_hla_exposure_list.append(allele)

    score_row['escape score prevalent alleles'] = my_new_delta if len(my_hla_list) != 0 else None
    score_row['prevalent alleles seen escape'] = my_hla_list
    score_row['exposure score prevalent alleles'] = my_new_exposure if len(my_hla_exposure_list) != 0 else None
    score_row['prevalent alleles seen exposure'] = my_hla_exposure_list
    return score_row


# -----------------------------df per allele------------------------------ <editor-fold>


# ------------------------------------------------------------------------ </editor-fold>


score_df = score_df.apply(smart_delta, axis=1)
score_df.to_csv(path_to_save + 'END_GAME_new_final_most_frequent.csv')
