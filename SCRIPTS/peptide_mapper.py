import pandas as pd
import subprocess
import sys
import os

list_args = list(sys.argv)
experiment_id = list_args[1]
analysis_type = list_args[2]
mode = list_args[3]
blastp_search_mode = list_args[4]

print('setting variables')
home_dir = f'/home/andrewstap/LINEs/Data_Analysis/{experiment_id}/FINAL_TABLES/'
script_dir = '/home/andrewstap/LINEs/SCRIPTS/LINE_COPIES_RMSK/AUTOMATIZE/POST_ALIGN_ANALYSIS/EXTRA/'
reference_home_dir = '/home/andrewstap/LINEs/REFERENCES/GENOME/GCHR38/RMSK/'
input_peptides = f'{reference_home_dir}MAPPINGS/BLASTP/united_peptidome_fasta.txt'
if analysis_type == 'FLANKS':
    if mode == 'REG_PEP':
        reference_fasta = f'{reference_home_dir}FLANKS/UNMASKED/LINES_L1_5500_ELSE_1500_FLANKS.fasta'
    elif mode == 'FOOTPR':
        reference_fasta = f'{reference_home_dir}MAPPINGS/FLANKS_200/FLANKS_mapping_200_footprints.fasta'
elif analysis_type == 'FULL':
    if mode == 'REG_PEP':
        reference_fasta = f'{reference_home_dir}LINES_L1_5500_ELSE_1500.fasta'
    elif mode == 'FOOTPR':
        reference_fasta = f'{reference_home_dir}MAPPINGS/LINES_200/LINES_mapping_200_footprints.fasta'
print('----------------')

print('making MAPPING directory')
analysis_folder = f'{home_dir}{analysis_type}/{mode}/MAPPING/'
if os.path.isdir(analysis_folder):
	pass
else:
	os.mkdir(analysis_folder)
print('----------------')


print('setting the blastp interpreter function')
def blastp_interpreter(file,align_type,mode_analysis,mapping_coords):
    blastp = pd.read_csv(file,sep='\t', header=None)
    blastp.columns = ['Query_Name','Query_Len','Reference_Name','Reference_Len',
                     'alignment_len','%_ident_pos','mismatch#','gapopen#','evalue',
                     'bitscore','query_start_align','qend_start_align','ref_start_align',
                     'ref_end_align','aligned_part_query','aligned_part_ref']
    approach2 = blastp[blastp['Query_Name'] == blastp['aligned_part_ref']]
    if align_type == 'STRINGENT': 
        a = approach2
    elif align_type == 'RELAXED':
        approach3 = blastp[((blastp['mismatch#'] == 0) & (blastp['gapopen#'] != 0) & (blastp['Query_Len'] == blastp['alignment_len']))]
        approach4 = blastp[((blastp['mismatch#'].isin([1])) & (blastp['Query_Len'] == blastp['alignment_len']))]
        a = pd.concat([approach2,approach3,approach4])
    a[['FAMILY','CHR','GENOMIC_START','GENOMIC_END_C']] = a['Reference_Name'].str.split(':',expand=True)
    if mode_analysis == 'FOOTPR':
        a[['GENOMIC_END','PEPTIDE_PLUS_STRAND','ORIENT','ORF']] = a['GENOMIC_END_C'].str.split('_',expand=True)
        a[['FOOTR_PEPTIDE_MAPPED','STRAND']] = a['PEPTIDE_PLUS_STRAND'].str.split('\(',expand=True)
        a[['STRAND','EMPTY']] = a['STRAND'].str.split('\)',expand=True)
        a['REFERENCE'] = a['FAMILY'].astype(str) + ':' + a['CHR'].astype(str) + ':' + a['GENOMIC_START'].astype(str) + ':' + a['GENOMIC_END'].astype(str)+ '_' + a['FOOTR_PEPTIDE_MAPPED'].astype(str) + '(' + a['STRAND'].astype(str) + ')'
        a.drop(['GENOMIC_END_C','EMPTY','PEPTIDE_PLUS_STRAND'],axis = 1,inplace=True)

    elif mode_analysis == 'REG_PEP':
        a[['GENOMIC_END','ELSE']] = a['GENOMIC_END_C'].str.split('\(',expand=True)
        a[['STRAND','ELSE2']] = a['ELSE'].str.split('\)_',expand=True)
        a[['ORIENT','ORF']] = a['ELSE2'].str.split('_',expand=True)
        a['REFERENCE'] = a['FAMILY'].astype(str) + ':' + a['CHR'].astype(str) + ':' + a['GENOMIC_START'].astype(str) + ':' + a['GENOMIC_END'].astype(str) + '(' + a['STRAND'].astype(str) + ')'
        a.drop(['GENOMIC_END_C','ELSE','ELSE2'],axis = 1,inplace=True)
    print(a)
    a['OVEREXP_START'] = ''
    a['OVEREXP_END'] = ''
    for each in list(a.index):
        a.loc[each,'OVEREXP_START'] = mapping_coords[a.loc[each,'REFERENCE']][0]
        a.loc[each,'OVEREXP_END'] = mapping_coords[a.loc[each,'REFERENCE']][1]
        a.loc[each,'GENOMIC_START'] = mapping_coords[a.loc[each,'REFERENCE']][2]
        a.loc[each,'GENOMIC_END'] = mapping_coords[a.loc[each,'REFERENCE']][3]

    a['PEPTIDE_STRAND'] = ''
    a.loc[((a['STRAND'] == '+') & (a['ORIENT'] == 'RV')),'PEPTIDE_STRAND'] = '-'
    a.loc[((a['STRAND'] == '+') & (a['ORIENT'] == 'FW')),'PEPTIDE_STRAND'] = '+'
    a.loc[((a['STRAND'] == '-') & (a['ORIENT'] == 'FW')),'PEPTIDE_STRAND'] = '-'
    a.loc[((a['STRAND'] == '-') & (a['ORIENT'] == 'RV')),'PEPTIDE_STRAND'] = '+'

    a['GENOME_LOC_START_PEPTIDE'] = ''
    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['ORF'] == 'ORF1')),'GENOME_LOC_START_PEPTIDE'] = (a['OVEREXP_START'].astype('int64') + 0) + (a['ref_start_align'].astype('int64') - 1)*3
    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['ORF'] == 'ORF2')),'GENOME_LOC_START_PEPTIDE'] = (a['OVEREXP_START'].astype('int64') + 1) + (a['ref_start_align'].astype('int64') - 1)*3
    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['ORF'] == 'ORF3')),'GENOME_LOC_START_PEPTIDE'] = (a['OVEREXP_START'].astype('int64') + 2) + (a['ref_start_align'].astype('int64') - 1)*3
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['ORF'] == 'ORF1')),'GENOME_LOC_START_PEPTIDE'] = (a['OVEREXP_END'].astype('int64') - 0) - (a['ref_end_align'].astype('int64') - 1)*3
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['ORF'] == 'ORF2')),'GENOME_LOC_START_PEPTIDE'] = (a['OVEREXP_END'].astype('int64') - 1) - (a['ref_end_align'].astype('int64') - 1)*3
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['ORF'] == 'ORF3')),'GENOME_LOC_START_PEPTIDE'] = (a['OVEREXP_END'].astype('int64') - 2) - (a['ref_end_align'].astype('int64') - 1)*3
    a['GENOME_LOC_END_PEPTIDE'] = ''
    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['ORF'] == 'ORF1')),'GENOME_LOC_END_PEPTIDE'] = (a['OVEREXP_START'].astype('int64') + 0) + (a['alignment_len'].astype('int64'))*3
    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['ORF'] == 'ORF2')),'GENOME_LOC_END_PEPTIDE'] = (a['OVEREXP_START'].astype('int64') + 1) + (a['alignment_len'].astype('int64'))*3
    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['ORF'] == 'ORF3')),'GENOME_LOC_END_PEPTIDE'] = (a['OVEREXP_START'].astype('int64') + 2) + (a['alignment_len'].astype('int64'))*3
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['ORF'] == 'ORF1')),'GENOME_LOC_END_PEPTIDE'] = (a['OVEREXP_END'].astype('int64') - 0) - (a['alignment_len'].astype('int64'))*3
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['ORF'] == 'ORF2')),'GENOME_LOC_END_PEPTIDE'] = (a['OVEREXP_END'].astype('int64') - 1) - (a['alignment_len'].astype('int64'))*3
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['ORF'] == 'ORF3')),'GENOME_LOC_END_PEPTIDE'] = (a['OVEREXP_END'].astype('int64') - 2) - (a['alignment_len'].astype('int64'))*3    

    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['GENOME_LOC_START_PEPTIDE'].astype('int64') < a['GENOMIC_START'].astype('int64'))),'GENOME_LOC_START_PEPTIDE'] = a['GENOMIC_START']
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['GENOME_LOC_START_PEPTIDE'].astype('int64') > a['GENOMIC_END'].astype('int64'))),'GENOME_LOC_START_PEPTIDE'] = a['GENOMIC_END']
    a.loc[((a['PEPTIDE_STRAND'] == '+') & (a['GENOME_LOC_END_PEPTIDE'].astype('int64') > a['GENOMIC_END'].astype('int64'))),'GENOME_LOC_END_PEPTIDE'] = a['GENOMIC_END']
    a.loc[((a['PEPTIDE_STRAND'] == '-') & (a['GENOME_LOC_END_PEPTIDE'].astype('int64') < a['GENOMIC_START'].astype('int64'))),'GENOME_LOC_END_PEPTIDE'] = a['GENOMIC_START']

    final_table = a[['Query_Name','FAMILY','CHR','GENOMIC_START','GENOMIC_END','PEPTIDE_STRAND','GENOME_LOC_START_PEPTIDE','GENOME_LOC_END_PEPTIDE']]
    final_table.sort_values(by=['Query_Name'],ascending=True,inplace=True)
    peptide_list = list(set(list(final_table['Query_Name'])))
    return [final_table,peptide_list]
print('----------------')


print('reading the continuous regions')
df = pd.read_csv(f'{home_dir}{analysis_type}/{mode}/PAIRS_VALUES_cont.txt', header=None)	
df[['copy_name','startstop']] = df[0].str.split("\):",expand=True)
df['copy_name'] = df['copy_name'] + ')'
df[['start','stop']] = df.startstop.str.split("-",expand=True)

df[['trash1','strand']] = df['copy_name'].str.split('\(',expand=True)
df[['strand','empty']] = df['strand'].str.split('\)',expand=True)
if mode == 'REG_PEP':
    df[['FAMILY','CHR','GENOME_START','GENOME_END']] = df['trash1'].str.split(':',expand=True)
elif mode == 'FOOTPR':
    df[['FAMILY','CHR','GENOME_START','GENOME_END_C']] = df['trash1'].str.split(':',expand=True)
    df[['GENOME_END','PEPTIDE']] = df['GENOME_END_C'].str.split('_',expand=True)    
df['OVEREXP_START'] = 0
df['OVEREXP_END'] = 0
df.loc[(df['strand'] == '+'),'OVEREXP_START'] = df['GENOME_START'].astype(int) + df['start'].astype(int)
df.loc[(df['strand'] == '+'),'OVEREXP_END'] = df['GENOME_START'].astype(int) + df['stop'].astype(int)
df.loc[(df['strand'] == '-'),'OVEREXP_END'] = df['GENOME_END'].astype(int) - df['start'].astype(int)
df.loc[(df['strand'] == '-'),'OVEREXP_START'] = df['GENOME_END'].astype(int) - df['stop'].astype(int)
if mode == 'REG_PEP':
    df['copy_name_overexp'] = df['FAMILY'].astype(str) + ':' + df['CHR'].astype(str) + ':' + df['OVEREXP_START'].astype(str) + ':' + df['OVEREXP_END'].astype(str) + '(' + df['strand'].astype(str) + ')'
elif mode == 'FOOTPR':
    df['copy_name_overexp'] = df['FAMILY'].astype(str) + ':' + df['CHR'].astype(str) + ':' + df['OVEREXP_START'].astype(str) + ':' + df['OVEREXP_END'].astype(str) + '_' + df['PEPTIDE'].astype(str) + '(' + df['strand'].astype(str) + ')'
overexp_startend_dict = dict(zip(df['copy_name_overexp'], list(zip(df['OVEREXP_START'],df['OVEREXP_END'],df['GENOME_START'],df['GENOME_END']))))

df_final = df[['copy_name','start','stop','copy_name_overexp']]
df_final['quality'] = 0
df_final['strand'] = '+'
df_final.to_csv(f'{analysis_folder}PAIRS_VALUES_overexpressed.bed',header=False,index=False,sep='\t')

bed_file = f'{analysis_folder}PAIRS_VALUES_overexpressed.bed'
fasta_file = f'{analysis_folder}PAIRS_VALUES_overexpressed.fasta'
amino_acid_file = f'{analysis_folder}PAIRS_VALUES_peptide_seq.fasta'
print('----------------')

print('converting continuous bed to fasta')
bedtools_cmd = f'bedtools getfasta -fi {reference_fasta} -bed {bed_file} -nameOnly -fo {fasta_file}'
subprocess.run(bedtools_cmd.split(' '))
# with open(fasta_file, 'r') as file_fasta_wrong:
#     lines_wrong = file_fasta_wrong.readlines()
# file_fasta_wrong.close()
# file_fasta_correct = open(fasta_file,'w')
# for each in lines_wrong:
#     file_fasta_correct.write(each.replace('()',''))
# file_fasta_correct.close()
print('----------------')

print('translating conitnuous fasta into peptides')
translator_cmd = f'python3 {script_dir}translator.py {fasta_file} {amino_acid_file}'
subprocess.run(translator_cmd.split(' '))
print('----------------')

print('creating the database for peptide sequences')
DB_NAME = f'{experiment_id}_{analysis_type}_{mode}'
database_cmd = f'makeblastdb -in {amino_acid_file} -dbtype prot -title {DB_NAME}'
subprocess.run(database_cmd.split(' '))
print('----------------')

print('running the blastp against our database')
blastp_cmd = f'blastp -query {input_peptides} -db {amino_acid_file} -out {analysis_folder}{analysis_type}_mappings.txt -outfmt \"6 qseqid qlen sseqid slen length pident mismatch gapopen evalue bitscore qstart qend sstart send qseq sseq\"'
os.system(blastp_cmd)
#subprocess.run(blastp_cmd.split(' '))
print('----------------')

print('running the blastp_interpreter function')
output = blastp_interpreter(f'{analysis_folder}{analysis_type}_mappings.txt',blastp_search_mode,mode,overexp_startend_dict)
output[0].to_csv(f'{analysis_folder}{experiment_id}_{analysis_type}_{mode}.txt',sep='\t',index=False)
peptide_list_file = open(f'{analysis_folder}{experiment_id}_{analysis_type}_{mode}_peptides.txt','w')
for each_peptide in output[1]:
	if each_peptide != output[1][-1]:
		peptide_list_file.write(f'{each_peptide}\n')
	else:
		peptide_list_file.write(f'{each_peptide}')
peptide_list_file.close()

