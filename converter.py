"""
this script will take in a list of all the BA2 mutations and find the coresponding mutations in other protien that are
in the ORF1a and ORF1b.

this will be don in a few steps:

1. make a dictinary that will the protiens in the orf as a key and the value is a tuppel with the stat and end indexies
as ites values
2. a df will be constructesd contaning all the infi for convertinf the mutations in the realavent protien
3. all the new proteins will be save in the refrece file format.
    * xlsx file named new_ref_df
    * a directory contaning all the fasta files for all the protiens



**********  NSP12 is a problematic one so it was delt with ocording to the overlping gap of 9 aa between orf1a and orf1b

by: sinai sacharen

"""


""" imports """
# ------------------------------------------------------------------------ <editor-fold>
from utiles import *

# ------------------------------------------------------------------------ </editor-fold>


""" paths """
# ------------------------------------------------------------------------ <editor-fold>

path_with_data = 'C:/Users/User/Documents/COVID-19/find_ref/'
path_to_save = 'C:/Users/User/Documents/COVID-19/omicron/'

if os.path.exists(path_to_save) and os.path.exists(path_with_data):
    print('found all paths  :) ')
else:
    print('cant find a path  :( ')
    exit()

# ------------------------------------------------------------------------ </editor-fold>


# ----------------------------------IMPORTED FILES--------------------------- <editor-fold>
ref_df = pd.read_excel(path_with_data+'new_ref_df.xlsx',index_col=0)
mutation_file = pd.read_csv(path_to_save+'mutations.txt',sep='\s+')
# ------------------------------------------------------------------------ </editor-fold>

# ------------------------------------FUNCTIONS--------------------------- <editor-fold>

# orf1a

dict_index = {}

input_handle = open(path_with_data +'ORF1b.fasta', "r")
for record in SeqIO.parse(input_handle, "fasta"):
    if 'ORF1a' in record.description:
        for idx in ref_df.index:
            if 'NS' in ref_df.at[idx,'GISAID name']:
                seq = ref_df.at[idx,'sequence']
                seq_index = record.seq.find(seq)
                if seq_index != -1:
                    dict_index[ref_df.at[idx,'GISAID name']] = [seq_index,seq_index+len(seq),'ORF1a']
    if 'ORF1b' in record.description:
        for idx in ref_df.index:
            if 'NS' in ref_df.at[idx,'GISAID name']:
                seq = ref_df.at[idx,'sequence']
                if 'NSP12' == ref_df.at[idx,'GISAID name']:
                    seq = seq[9:]
                seq_index = record.seq.find(seq)
                if seq_index != -1:
                    dict_index[ref_df.at[idx,'GISAID name']] = [seq_index,seq_index+len(seq),'ORF1b']
                    print(ref_df.at[idx,'GISAID name'])

input_handle.close()

conveter_df = pd.DataFrame.from_dict(dict_index).T


def transform_index(orfMutation_row,orf_index_dict):
    """
    this function will get the mutation from outbrakeinfo and will
    return the mutation in the NSP format
    :param orfStart:
    :return: mutations
    """
    orf = orfMutation_row[0]
    orfMutation = orfMutation_row[1]
    pos = int(orfMutation[1:-1])
    for key,p in orf_index_dict.items():
        if orf == 'ORF1a' and p[2] == 'ORF1a':
            if pos in range(p[0],p[1]+1):
                nsp_pos = pos-p[0]
                mut = orfMutation[-1]
                return key+'_'+orfMutation[0]+str(nsp_pos)+mut
        if orf == 'ORF1b' and p[2] == 'ORF1b':
            if pos in range(p[0]+9,p[1]+10):
                if key == 'NSP12':
                    nsp_pos = pos-p[0] + 9
                else:
                    nsp_pos = pos - p[0]
                mut = orfMutation[-1]
                return key+'_'+orfMutation[0]+str(nsp_pos)+mut

    return orf+'_'+orfMutation

mutation_file['nsp mutation'] = mutation_file.apply(transform_index,orf_index_dict=dict_index,axis=1)
mutation_file[['nsp', 'nsp mut']] = mutation_file['nsp mutation'].str.split('_', expand=True)

# ------------------------------------------------------------------------ </editor-fold>


## now we creat a new omicron based refrece directory
# find_ref/
# ├── new_ref_df.xlsx
# ├── sars_protien_fasta
    #       ├── E_fasta.fa
    #       ├── key doesn't exist_fasta.fa
    #       ├── M_fasta.fa
    #       ├── N_fasta.fa etc......


# ------------------------------------------------------------------------ <editor-fold>
omicron_ref_df = ref_df.copy()
for idx in ref_df.index:# for mutation
    protin = ref_df.at[idx,'alternate name']
    sequnce = ref_df.at[idx,'sequence']

    # mutation inducer
    mut_list = list(mutation_file[mutation_file['nsp']==protin]['nsp mut'])
    if len(mut_list) == 0:
        continue
    # update omicron df
    omicron_seq = mutation_iducer(sequnce,mut_list)
    omicron_seq = omicron_seq.replace('-','')
    # creat new fasta file on omicron directory
    omicron_ref_df.at[idx,'sequence'] = omicron_seq
    gisaid_name = ref_df.at[idx,'GISAID name']
    with open(path_to_save+'omicron_prot/'+gisaid_name+'_fasta.fa', "w") as file1:
        # Writing data to a file
        file1.write(">Omicron "+ gisaid_name+" \n")
        file1.write(omicron_seq)



omicron_ref_df.to_excel(path_to_save+'omicron_ref_df.xlsx')
mutation_file.to_excel(path_to_save+'mutation_df.xlsx')

# ------------------------------------------------------------------------ </editor-fold>


