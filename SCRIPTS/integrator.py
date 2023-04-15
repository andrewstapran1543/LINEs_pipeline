import pandas as pd
import os
import sys

list_args = list(sys.argv)
input_len_file = str(list_args[1])
experiment_id = str(list_args[2])
analysis_type = str(list_args[3])
mode = str(list_args[4])

input_len_df = pd.read_csv(input_len_file, sep=';')
subset1 = input_len_df[input_len_df['EXPERIMENT'] == experiment_id]
read_length = dict(zip(subset1.SAMPLE,subset1.READNUM))


home_dir = f'/home/andrewstap/LINEs/Data_Analysis/{experiment_id}/'
list_dirs = [f'{home_dir}FINAL_TABLES',f'{home_dir}FINAL_TABLES/FLANKS',f'{home_dir}FINAL_TABLES/FULL']
for analysis_each in ['FULL','FLANKS']:
    for mode_each in ['REG_PEP','FOOTPR']:
        list_dirs.append(f'{home_dir}FINAL_TABLES/{analysis_each}/{mode_each}')
for each in list_dirs:
    if os.path.isdir(each):
        pass
    else:
        os.mkdir(each)

def reader(SRA, len_dict):
    if mode == 'REG_PEP':
        path_to_read = f'{home_dir}{SRA}/{analysis_type}/regions_to_peptides/{SRA}_reg_to_pep_coverage.txt'
    elif mode == 'FOOTPR':
        path_to_read = f'{home_dir}{SRA}/{analysis_type}/peptides_to_regions/{SRA}_pep_to_reg_coverage.txt'
    print('Reading DF')
    first_df2 = pd.read_csv(path_to_read,sep='\t',header=None)
    first_df = first_df2[first_df2[2] != 0]
    print('Calculating coverage')
    first_df[3] = first_df[2] / read_length[SRA]
    print('Setting position names')

    first_df[4] = first_df[0].astype(str) + '_' + first_df[1].astype(str)
    print('Subsetting dataframe')
    input_df = first_df[[3,4]].copy()
    input_df.set_index(4,inplace=True)
    print(input_df)
    return input_df


for i in list(read_length.keys()):
    print(i)
    if i == list(read_length.keys())[0]:
        input_df = reader(i,read_length)
    else:
        new_df = reader(i,read_length)
        input_df = pd.concat([input_df,new_df],axis=1)
    print('--------------')

final_dataframe = input_df.fillna(0)
final_dataframe.columns = list(read_length.keys())
final_dataframe.to_csv(f'{home_dir}FINAL_TABLES/{analysis_type}/{mode}/final_table.csv',sep='\t')

