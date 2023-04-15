import subprocess
import sys

list_args = list(sys.argv)

input_len_file = str(list_args[1])
input_pairs_file = str(list_args[2])
experiment_id = str(list_args[3])
pairs_threshold = float(list_args[4])
region_length_threshold = int(list_args[5])
blastp_mode = str(list_args[6])


script_dir = '/home/andrewstap/LINEs/SCRIPTS/LINE_COPIES_RMSK/AUTOMATIZE/POST_ALIGN_ANALYSIS/EXTRA/'

def analysis_runner(analysis,mode):
	print('Integrating the tables into FINAL_CSV')
	integrator_cmd = f'python3 {script_dir}integrator.py {input_len_file} {experiment_id} {analysis} {mode}'
	subprocess.run(integrator_cmd.split(' '))
	print('-------------------')

	print('Calulating the continuous regions from FINAL_CSV')
	calculator_cmd = f'python3 {script_dir}calculator.py {experiment_id} {analysis} {input_pairs_file} {mode} {pairs_threshold} {region_length_threshold}'
	subprocess.run(calculator_cmd.split(' '))
	print('-------------------')

	print('Running the PEPTIDE MAPPING - getting the final peptide list')
	blastp_runner_cmd = f'python3 {script_dir}peptide_mapper.py {experiment_id} {analysis} {mode} {blastp_mode}'
	subprocess.run(blastp_runner_cmd.split(' '))
	print('-------------------')

for each in ['FULL','FLANKS']:
	for type_exp in ['REG_PEP','FOOTPR']:
		print(f'POST_ALIGNMENT analysis: {experiment_id}\nReference files: {each} LINE copies\nAnalysis mode: {type_exp}\n')
		analysis_runner(each, type_exp)
print("COMPLETE")


