"""
this script will take in the new_ref_df xlsx file and the all prot no X file
then it will creat a new all prot no X file with only BA2 related sequences
this will be sent to Mafft MSA and will then be used for the HLA landscape calculations
and then mutations peptide creation

"""

""" imports """
# ------------------------------------------------------------------------ <editor-fold>
from utiles import *
from Tool_File import *

import sys

import io

# ------------------------------------------------------------------------ </editor-fold>


""" paths """
# ------------------------------------------------------------------------ <editor-fold>

path_with_data = '/home/sacharen/Documents/' + '/only_wonted_prots/'
path_to_save = '/home/sacharen/Documents/omicron/'

if os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()

# ------------------------------------------------------------------------ </editor-fold>


# ----------------------------------IMPORTED FILES--------------------------- <editor-fold>

file_og_name = 'allprot0628'

ref_df = pd.read_excel(path_to_save + 'omicron_ref_df.xlsx', index_col=0)
mutation_file = pd.read_csv(path_to_save + 'mutations.txt', sep='\s+')
mutation_ref = pd.read_excel(path_to_save + 'mutation_df.xlsx', index_col=0)

ref_df.set_index(['GISAID name'], inplace=True)


# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------------FUNCTIONS--------------------------- <editor-fold>
def mutation_finder(spike_seq, mut_list):
    """
    this function will take in a og seq and a mutation list
    will returen True if all mutation are in the sequence False if not
    :param spike_seq: og spike sequence
    :param mut_list: a list of mutations
    :return: the mutated seq
    """
    flag = True
    spike_list = list(spike_seq)
    for m in mut_list:
        pos = int(m[1:-1])
        if m[0] != spike_list[pos - 1]:
            flag = False
            break
    return flag

# ------------------------------------------------------------------------ </editor-fold>





p_list = ['Spike', 'N']

for P in p_list:
    all_seq_list = []
    no_filter = False
    mut_list = list(mutation_file[mutation_file['nsp'] == P]['nsp mut'])
    if len(mut_list) == 0:
        no_filter = True
    input_handle = open(path_to_save + P + '_no_x' + '.fa', "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        if no_filter == True:
            all_seq_list.append(record)

        else:
            if mutation_finder(record.seq, mut_list):
                all_seq_list.append(record)

    input_handle.close()
    SeqIO.write(all_seq_list, path_to_save + P + '_no_x_omicron' + '.fa', 'fasta')

for P in p_list:

    input_file = P + '_no_x_omicron' + '.fa'
    output_file = path_to_save + P + "_maffet_alinged_omicron" + ".fa"
    cd_c = 'cd ' + path_to_save + '; '
    command = cd_c + "mafft --6merpair --keeplength --thread -1 --anysymbol --addfragments  " + input_file + " /home/sacharen/Documents/omicron/omicron_prot/" + P + "_fasta.fa > " + output_file
    sed_c = "'s/*//g'"
    comm = 'sed -i ' + sed_c + ' ' + path_to_save + input_file
    subprocess.check_output('%s' % 'cd', shell=True)

    subprocess.check_output('%s' % comm, shell=True)

    try:
        subprocess.check_output('%s' % command, shell=True)
    except subprocess.CalledProcessError:
        continue

    print(P,
          'done-----------------------------------------------------------------------------------------------------')

count = 0
for P in p_list:
    ref_sequence = ref_df.loc[P, 'sequence']

    list_of_recoreds = []

    input_handle = open(path_with_data + P + '_maffet_alinged_omicron'  '.fa', "r")
    for record in SeqIO.parse(input_handle, "fasta"):
        if '>' + P in record.description:
            continue
        list_of_recoreds.append(record)
    input_handle.close()

    ####################################################  PFM and entropy
    p, c = simple_list_PFM(list_of_recoreds, ref_sequence=ref_sequence)
    p.to_excel(path_to_save+'pfm/' + P + 'PFM.xlsx')
    c.to_excel(path_to_save+'pfm/' + P + 'COUNT.xlsx')
    simple_entropy_function(path_to_save=path_to_save+'pfm/', save_file_name=P + '_entropy_matrix', pfm_df=p)

    print('PFM done!')

    # -----------------------------------entropy no bias---------------------- <editor-fold>

    # when we will wont it to be up to date

    date_list = pd.date_range(start="3/3/2020", end=datetime.datetime.today()).tolist()

    # date_list = pd.date_range(start="12/12/2019",end="10/12/2020")
    new_date_list = []
    for i in date_list:
        month_year = str(i.month) + '-' + str(i.year)
        if month_year not in new_date_list:
            new_date_list.append(month_year)

    date_list = pd.to_datetime(new_date_list).strftime("%Y-%b").tolist()
    date_list.pop(-1)

    list_for_save = []
    protien_dict = {}
    for month in date_list:
        protien_dict[month] = []

    for rec in list_of_recoreds:
        try:
            dt = find_date(rec.description)
        except NameError:
            continue
        # print(type(dt))
        NoneType = type(None)
        if type(dt) != NoneType:
            dt = dt.strftime("%Y-%b")
            if dt in date_list:
                protien_dict[dt].append(rec)

    for key in protien_dict.keys():
        recored_temp_list = protien_dict[key]

        random.shuffle(recored_temp_list)
        list_for_save = list_for_save + recored_temp_list[:3901]

        ##### this part is for monthly entropy ####
        p, c = simple_list_PFM(recored_temp_list, ref_sequence=ref_sequence)
        p['date'] = p['A'].apply(lambda x: key)
        p['protein'] = p['A'].apply(lambda x: P)
        if count == 0:
            count += 1
            monthly_df = p.copy()
        monthly_df = pd.concat([monthly_df, p])

    print(P, 'monthly PFM done!!!!!!!!!!!!!')

    ####################################################  PFM and entropy
    p, c = simple_list_PFM(list_for_save, ref_sequence=ref_sequence)
    p.to_excel(path_to_save+'pfm/' + P + '_entropy no bias PFM.xlsx')
    c.to_excel(path_to_save+'pfm/' + P + '_entropy no bias COUNT.xlsx')
    simple_entropy_function(path_to_save=path_to_save, save_file_name=P + '_entropy no bias entropy_matrix',
                            pfm_df=p)

    print(P, 'entropy no bias PFM done!')

monthly_df.to_csv(path_to_save + 'monthly_df_all_prot.csv')

