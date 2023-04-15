import pandas as pd
import numpy as np
import sys

list_args = list(sys.argv)
experiment_id = str(list_args[1])
analysis_type = str(list_args[2])
input_pairs_file = str(list_args[3])
mode = str(list_args[4])
pairs_threshold = float(list_args[5])
region_length_threshold = int(list_args[6])

home_dir = f'/home/andrewstap/LINEs/Data_Analysis/{experiment_id}/'
input_pairs_df = pd.read_csv(input_pairs_file,sep=';')
input_pairs_df_subset1 = input_pairs_df[input_pairs_df['EXPERIMENT'] == experiment_id]
normal_samples = list(input_pairs_df_subset1['SAMPLE_NORM'])
tumor_samples = list(input_pairs_df_subset1['SAMPLE_TUMOR'])
pair_names = list(input_pairs_df_subset1['PAIR_NAME'])

pairs_norm_tumor = {}
for index,element in enumerate(pair_names):
    pairs_norm_tumor[element] = [normal_samples[index],tumor_samples[index]]

print('Setting Variables')
path_to_files = f'{home_dir}FINAL_TABLES/{analysis_type}/{mode}'
file_initial_csv = f'{path_to_files}/final_table.csv'
#file_MEANs = f'{path_to_files}/MEAN.csv'
file_MEDIANs = f'{path_to_files}/MEDIAN.csv'
file_PAIRS = f'{path_to_files}/PAIRS_VALUES.csv'
#file_MEANs_towrite = f'{path_to_files}/MEAN_cont.txt'
file_MEDIANs_towrite = f'{path_to_files}/MEDIAN_cont.txt'
file_PAIRS_towrite = f'{path_to_files}/PAIRS_VALUES_cont.txt'
print('------------')

print('reading df')
final_table = pd.read_csv(file_initial_csv,sep='\t',index_col='4')
print('------------')


# print('calculating mean')
# final_table['NORM_mean'] = final_table[normal_samples].mean(axis=1)
# final_table['TUMOR_mean'] = final_table[tumor_samples].mean(axis=1)
# final_table[['NORM_mean','TUMOR_mean']].copy().to_csv(file_MEANs,sep='\t')
# print('------------')



print('calculating median')
final_table['NORM_median'] = final_table[normal_samples].median(axis=1)
final_table['TUMOR_median'] = final_table[tumor_samples].median(axis=1)
final_table[['NORM_median','TUMOR_median']].copy().to_csv(file_MEDIANs,sep='\t')
print('------------')



print('conducting pairwise comparison')
print('caculating difference at each position for each pair')
for each in list(pairs_norm_tumor.keys()):
    print(each)
    final_table[each] = 0
    norm = pairs_norm_tumor[each][0]
    tumor = pairs_norm_tumor[each][1]
    final_table.loc[final_table[norm] < final_table[tumor],each] = 1

print('subsetting the status columns')
paired_values = final_table[list(pairs_norm_tumor.keys())].copy()
print('calculating mean value')
paired_values['pairs_average'] = paired_values[list(pairs_norm_tumor.keys())].mean(axis=1)
print('subsetting based on average threshold')
subset = paired_values[paired_values['pairs_average'] > pairs_threshold]
subset[['pairs_average']].copy().to_csv(file_PAIRS,sep='\t')
print('------------')





print('setting function for identifying continuous regions')
def continuous_regions(list1,list2,threshold,file_to_write):
    file_writing = open(file_to_write,'w')
    regions = []
    copy_name = list2.iloc[0]
    start_pos = list1.iloc[0]
    for index,value in enumerate(list1):
        if int(index + 1) < len(list2):
            if int(int(list1.iloc[index+1]) - int(value)) == 1:
                continue 
            elif int(int(list1.iloc[index+1]) - int(value)) != 1:
                end_pos = value
                if int(end_pos) - int(start_pos) > threshold:
                    print(f'{copy_name}:{start_pos}-{end_pos}')
                    file_writing.write(f'{copy_name}:{start_pos}-{end_pos}\n')
                    regions.append(f'{copy_name}:{start_pos}-{end_pos}')
                copy_name = list2.iloc[index+1]
                start_pos = list1.iloc[index+1]
        elif int(index + 1) >= len(list2)-1:
            print('Complete')
    file_writing.close()
    return regions


print('setting function for returning the file with regions for mean/median')
def cont_regions_interpreter_MEAN_MEDIAN(file,threshold,file_to_write):
    print('reading table')
    DF = pd.read_csv(file,sep='\t')
    print('-------------')
    print('renaming columns')
    DF.columns = ['CopyNamePos','NORM','TUMOR']
    print('-------------')
    print('setting status on overexpression')
    DF['status'] = np.where(DF['TUMOR'] > DF['NORM'],1,0)
    print('-------------')
    print('splitting the name column')
# #for REGIONS TO PEPTIDES approach
#     if mode == 'REG_PEP':
#         DF[['copy_name','position']] = DF.CopyNamePos.str.split("\)_",expand=True)
#         DF['copy_name'] = DF['copy_name'] + ')'
# #for PEPTIDES REGIONS approach
#     elif mode == 'FOOTPR':
#         DF[['copy_name''position']] = DF.CopyNamePos.str.split("\)_",expand=True)
#         DF['copy_name'] = DF['copy_name'] + ')'
    DF[['copy_name','position']] = DF.CopyNamePos.str.split("\)_",expand=True)
    DF['copy_name'] = DF['copy_name'] + ')'
    DF = DF.dropna(axis=0)
    DF = DF.astype({'position':'int'})

    print('-------------')
    print('sorting the table by copy name and position')
    DF2 = DF.sort_values(['copy_name', 'position'],ascending = [True, True])
    print('-------------')
    print('creating lists for identifying continuous regions')
    DF_positions = DF2[DF2['status'] == 1]['position']
    DF_copies = DF2[DF2['status'] == 1]['copy_name']
    print('-------------')
    print('CONTINUOUS REGIONS ARE THE FOLLOWING:')
    return continuous_regions(DF_positions,DF_copies,threshold,file_to_write)


print('setting function for returning the file with regions for intervals')
def cont_regions_interpreter_PAIRS(file,threshold,file_to_write):
    print('reading table')
    DF = pd.read_csv(file,sep='\t')
    print('-------------')
    print('renaming columns')
    DF.columns = ['CopyNamePos','STATUS']
    print('-------------')
    print('splitting the name column')
# #for REGIONS TO PEPTIDES approach    
#     if mode == 'REG_PEP':
#         DF[['copy_name','position']] = DF.CopyNamePos.str.split("\)_",expand=True)
#         DF['copy_name'] = DF['copy_name'] + ')'
# #for PEPTIDES REGIONS approach
#     elif mode == 'FOOTPR':
#         DF[['copy_name','peptide','position']] = DF.CopyNamePos.str.split("\)_",expand=True)
#         DF['copy_name'] = DF['copy_name'] + ')_' + DF['peptide'] + ')'
    DF[['copy_name','position']] = DF.CopyNamePos.str.split("\)_",expand=True)
    DF['copy_name'] = DF['copy_name'] + ')'
    DF = DF.dropna(axis=0)
    DF = DF.astype({'position':'int'})

    print('-------------')
    print('sorting the table by copy name and position')
    DF2 = DF.sort_values(['copy_name', 'position'],ascending = [True, True])
    print('-------------')
    print('creating lists for identifying continuous regions')
    DF_positions = DF2['position']
    DF_copies = DF2['copy_name']
    print('-------------')
    print('CONTINUOUS REGIONS ARE THE FOLLOWING:')
    return continuous_regions(DF_positions,DF_copies,threshold,file_to_write)
print('------------')



print('running the functions, obtaining continuous regions')
# print('MEAN')
# cont_regions_interpreter_MEAN_MEDIAN(file_MEANs,region_length_threshold,file_MEANs_towrite)
print('MEDIAN')
cont_regions_interpreter_MEAN_MEDIAN(file_MEDIANs,region_length_threshold,file_MEDIANs_towrite)
print('PAIRS')
cont_regions_interpreter_PAIRS(file_PAIRS,region_length_threshold,file_PAIRS_towrite)
print('------------')
print('done')
