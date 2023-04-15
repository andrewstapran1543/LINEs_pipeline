import os
import sys
import pandas as pd
import subprocess

list_args = list(sys.argv)
norm_tumor_file = list_args[1]
experiment_id = list_args[2]

norm_tumor_df = pd.read_csv(norm_tumor_file,sep=';')
norm_tumor_df_subset1 = norm_tumor_df[norm_tumor_df['EXPERIMENT'] == experiment_id]
normal_samples = list(norm_tumor_df_subset1['SAMPLE_NORM'])
tumor_samples = list(norm_tumor_df_subset1['SAMPLE_TUMOR'])
bam_file_loc_path = f'/home/andrewstap/LINEs/Data_Analysis/{experiment_id}/'
TE_tool = '/home/andrewstap/LINEs/TOOLS/TEtranscripts-2.2.3/bin/TEtranscripts'
TE_file_path = "/home/andrewstap/LINEs/REFERENCES/GENOME/TE_tool/GRCh38_GENCODE_rmsk_TE.gtf" 
GTF_file_path = "/home/andrewstap/LINEs/REFERENCES/GENOME/TE_tool/gencode.v38.annotation.gtf"
output_directory = f'{bam_file_loc_path}TETRANSCRIPT/'
if os.path.isdir(output_directory):
	print('Output folder exists already')
else:
	os.mkdir(output_directory)



TETranscript_cmd = f'{TE_tool} -t'
for each in tumor_samples:
	bam_file = f'{bam_file_loc_path}{each}/STAR_GENOME/{each}_Aligned.out.bam'
	TETranscript_cmd += ' '
	TETranscript_cmd += f'{bam_file}'
TETranscript_cmd += f' -c'
for each in normal_samples:
	bam_file = f'{bam_file_loc_path}{each}/STAR_GENOME/{each}_Aligned.out.bam'
	TETranscript_cmd += ' '
	TETranscript_cmd += f'{bam_file}'
TETranscript_cmd += f' --GTF {GTF_file_path} --TE {TE_file_path} --format BAM --mode multi --outdir {output_directory} --verbose'
subprocess.run(TETranscript_cmd.split(' '))
print('complete')