"""

this script contains the main pipeline used for HLA escape and exposure analysis.


"""


""" imports """
# ------------------------------------------------------------------------ <editor-fold>

from Tool_File import *
import math

# ------------------------------------------------------------------------ </editor-fold>

""" paths """
# ------------------------------------------------------------------------ <editor-fold>

path_to_tool = '/home/sacharen/netMHCpan-4.1/'
path_with_data = '/home/sacharen/Documents/omicron/'
path_to_save = '/home/sacharen/Documents/omicron/omicronHla/'
if os.path.exists(path_to_tool) and os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()


# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------FUNctIONs--------------------------------- <editor-fold>

def connect_patient_to_pep(peptide_f):
    """
    this function takes in a peptide and finds out were it is lokated in the protien
    and then it can find the ralevant mutation and petient that are connected to it
    :param peptide_f: from the df_f
    global: orf1ab
    :return: the patient and the mutation
    """
    global orf1ab_seq
    global patiant_dict
    return_list = []

    pep_start = find_peptide_position(peptide=peptide_f, ref_seq=orf1ab_seq)
    pep_end = pep_start + len(peptide_f)
    for key, val in patiant_dict.items():
        for mut in val:
            mut_pos = int(mut[1:-1])
            if mut_pos in range(pep_start, pep_end + 1):
                return_list.append({key: mut})

    return return_list


def scale_function(number):
    e = math.e
    return 1 / (e ** (number - e))


def proces_out_MHC_xls(pro_df):
    """
    this function takes in the output xlsx file from net MHCpan and
    cleans it up so ut can be merged into the main df
    :param pro_df:
    :return: clean df
    """
    new_col_dict = {}

    current_hla = ''
    for col in pro_df.columns:
        print(col)
        if col[0].startswith('HLA'):
            current_hla = col[0]
        if current_hla == '':
            new_col_dict[col] = col[1]
        if col[1] == 'EL_Rank':
            new_col_dict[col] = current_hla

    pro_df.columns = pro_df.columns.to_flat_index()
    pro_df.rename(columns=new_col_dict, inplace=True)
    for col in pro_df.columns:
        if type(col) is tuple:
            pro_df.drop(col, axis=1, inplace=True)

    allele_cols_scale_part = []
    for dcs in pro_df.columns:
        if 'HLA-' in dcs:
            allele_cols_scale_part.append(dcs)
    print(pro_df)
    pro_df[allele_cols_scale_part] = pro_df[allele_cols_scale_part].apply(scale_function)
    pro_df['mean'] = pro_df[allele_cols_scale_part].mean(axis=1)

    return pro_df


def read_in_allels(prealent_allel_file_and_path):
    """
    tis function will read in the file and return all the allelels in it in the following format:
    'HLA-B40:01,HLA-B58:01,HLA-B15:01'
    so the netHmcpan predicter can read it in
    :param prealent_allel_file_and_path:
    :return: string
    """
    with open(prealent_allel_file_and_path) as file_in:
        lines = ''
        for line in file_in:
            temp_line = line.rstrip() + ','
            if len(temp_line) < 2:
                continue
            else:
                lines += temp_line
    return lines[:-1]


def is_not_og_aa(aa_id_f, ref_s):
    pos = int(aa_id_f[:-1])
    aa = aa_id_f[-1:]
    if ref_s[pos - 1] == aa:
        return False
    else:
        return True

# ------------------------------------------------------------------------ </editor-fold>


# ------------------------------------------------------------------------ </editor-fold>

""" global variables """
# ------------------------------------------------------------------------ <editor-fold>
prevalent_allel_file = path_with_data + 'prevalent_alleles'
# full_alele_file = pd.read_csv(path_with_data + 'all_netMHCpan_alleles.txt', sep='\s+', header=None)
path_with_all_protien_fa = '/home/sacharen/Documents/omicron/omicron_prot/'
ref_df = pd.read_excel('/home/sacharen/Documents/omicron/omicron_df.xlsx')
ref_df.set_index(['GISAID name'], inplace=True)
file_with_mut_peptides = path_with_data+'temp_mut_prp_file.txt'



# ------------------------------MAIN for all prevalent alleles------------- <editor-fold>
for p in ['NSP10', 'NSP9', 'NSP8', 'NSP7', 'NSP6', 'NSP5', 'NSP4', 'NSP3', 'NSP2', 'NSP1', 'NS7b', 'NSP11', 'NSP16',
          'NSP15', 'NSP14'
    , 'NSP13', 'NSP12', 'N', 'NS8', 'NS7a', 'NS6', 'M', 'E', 'NS3', 'Spike']:

    print('processing protein '+p+'...........................')
    file_and_path_with_fasta = path_with_all_protien_fa + p + "_fasta.fa"

    # creat temp xlsx file
    # with open(path_with_data + 'temp.xlsx', 'a') as fp:
    #     pass


    open(path_with_data + 'temp.xlsx', 'a').close()

    temp_xlsx_file_and_path_f = path_with_data + 'temp.xlsx'

    HLA_str = read_in_allels(prevalent_allel_file)

    # peptide_length = "8,9,10,11"  # -l
    peptide_length = "9"  # -l

    command = path_to_tool + "netMHCpan " + " -a " + HLA_str + " -f " + file_and_path_with_fasta + " -l " + "9" + " -xls " + " -xlsfile " + temp_xlsx_file_and_path_f

    subprocess.run(command, shell=True)

    df_f = pd.read_csv(temp_xlsx_file_and_path_f, sep='\t', header=[0, 1])

    df_f = proces_out_MHC_xls(pro_df=df_f)
    df_f.to_csv(path_to_save + p + 'hla_coverage.csv')
    # delete temp xlsx file
    os.remove(temp_xlsx_file_and_path_f)

    print('prtotien: ' + p + '   .................................')
    # P = 'M'
    prot_seq = ref_df.loc[p, 'sequence']

    # ------------------------------------------------------------------------ </editor-fold>

    # -------------------------------take all 9 mers from hla coverage file --- <editor-fold>

    df_f['len pep'] = df_f['Peptide'].apply(lambda x: len(x))
    df_f = df_f[df_f['len pep'] == 9]
    df_f.drop('len pep', axis=1, inplace=True)
    unique_pep_list = list(df_f['Peptide'])

    # ------------------------------------------------------------------------ </editor-fold>

    # -----------------------here a make a mut list from the pfm no bias file------ <editor-fold>
    mut_pdf_df = pd.read_excel(path_with_data + p + '_entropy no bias PFM.xlsx')
    mut_pdf_df = pd.melt(mut_pdf_df, id_vars=['Unnamed: 0'],
                         value_vars=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                                     'T',
                                     'W', 'Y', 'V', '-'])

    mut_pdf_df = mut_pdf_df[mut_pdf_df['value'] > 0.0001]  ###### the mut cut off######
    mut_pdf_df['mus_id'] = mut_pdf_df['Unnamed: 0'].astype(str) + mut_pdf_df['variable']

    mut_pdf_df['for list'] = mut_pdf_df['mus_id'].apply(is_not_og_aa, ref_s=prot_seq)
    mut_pdf_df = mut_pdf_df[mut_pdf_df['for list']]

    mut_list = list(mut_pdf_df['mus_id'])

    # ----------------------------------------------------------------------------- </editor-fold>

    # creat the main df with all data  ----------------------------------------------------------------------

    main_df = pd.DataFrame(columns=['mutations position', 'mutations aa', 'Peptide', 'start pep', 'end pep'])

    main_idx = 0
    for mut in mut_list:
        mut_position = int(mut[:-1])
        mutation_aa = mut[-1:]
        for pep in unique_pep_list:
            if is_mutation_in_peptide(mut_pos_f=mut_position, pep_f=pep, ref_sequence_f=prot_seq):
                start_pos = find_peptide_position(pep, ref_seq=prot_seq)
                end_pos = start_pos + len(pep) - 1
                main_df.at[main_idx, 'mutations position'] = mut_position
                main_df.at[main_idx, 'mutations aa'] = mutation_aa
                main_df.at[main_idx, 'Peptide'] = pep
                main_df.at[main_idx, 'start pep'] = start_pos
                main_df.at[main_idx, 'end pep'] = end_pos
                main_idx += 1
    main_df = pd.merge(main_df, df_f, how="left", on=["Peptide"])

    # creat all the mutated peptides in the merged df
    # ----------------------------------------------------------------------
    main_df['MUT peptides'] = main_df.apply(lambda x: creat_the_mutated_peptide(x['Peptide'], x['mutations position'], x['mutations aa'],ref_sequence_f=prot_seq), axis=1)

    up_dated_mut_list = list(main_df['MUT peptides'])
    creat_pep_file_from_list(pep_list_f=up_dated_mut_list, file_path_and_name_f=file_with_mut_peptides)

    HLA_str = read_in_allels(prevalent_allel_file)
    command = path_to_tool + "netMHCpan " + " -a " + HLA_str + " -p " + file_with_mut_peptides + " -xls " + " -xlsfile " + temp_xlsx_file_and_path_f

    subprocess.check_output('%s' % command, shell=True)

    mut_predict_df = pd.read_csv(temp_xlsx_file_and_path_f, sep='\t', header=[0, 1])

    mut_predict_df = proces_out_MHC_xls(pro_df=mut_predict_df)
    dict_mut_name = {}
    for col in mut_predict_df.columns:
        if 'HLA' in col:
            dict_mut_name[col] = 'MUT ' + col

    dict_mut_name['Peptide'] = "MUT peptides"
    dict_mut_name['mean'] = "MUT mean"

    mut_predict_df.rename(columns=dict_mut_name, inplace=True)

    mut_predict_df.to_excel(path_to_save + p + '_mut_df.xlsx')


    final_merged_df = pd.merge(main_df, mut_predict_df, how="left", on=["MUT peptides"])

    final_merged_df.to_excel(path_to_save + p + '_final_df.xlsx')

# ------------------------------------------------------------------------ </editor-fold>



