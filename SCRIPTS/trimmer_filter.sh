#!/bin/bash

list_samples=$1
seq_type=$2
experiment_id=$3
length_filter=$4

samples_total=$(cat $list_samples)
samples_number=$(cat $list_samples | wc -l)

mkdir -p "/ngs/2023/02_rnaseq_datasets_Stapran/${experiment_id}/ANALYSIS"
analysis_path="/ngs/2023/02_rnaseq_datasets_Stapran/${experiment_id}/ANALYSIS/"
storage_path="/ngs/2023/02_rnaseq_datasets_Stapran/${experiment_id}/"


for i in $(seq 1 $samples_number); do
	SRR=$(echo $samples_total | cut -d' ' -f$i)
	echo "analyzing sample ${SRR} in experiment ${experiment_id}"
	mkdir ${analysis_path}${SRR}_trimmed_dedup
	folder_path="${analysis_path}${SRR}_trimmed_dedup/"

	if [ "$seq_type" == "PE" ]; then
		echo "running fastp trimming and filtering for sample ${SRR}"
		fastp --trim_poly_g --trim_poly_x --length_required $length_filter --average_qual 30 \
			-i ${storage_path}${SRR}_1.fastq.gz \
			-I ${storage_path}${SRR}_2.fastq.gz \
			-o ${folder_path}${SRR}_trimmed_1.fastq.gz -O ${folder_path}${SRR}_trimmed_2.fastq.gz -V \
			-h ${folder_path}fastp.html -j ${folder_path}fastp.json

		echo "unzipping the trimmed and filtered files for sample ${SRR} for prinseq"
		gunzip ${folder_path}${SRR}_trimmed_1.fastq.gz
		gunzip ${folder_path}${SRR}_trimmed_2.fastq.gz

		echo "prinseq deduplication for sample ${SRR}"
		prinseq-lite.pl -verbose -fastq ${folder_path}${SRR}_trimmed_1.fastq \
			-fastq2 ${folder_path}${SRR}_trimmed_2.fastq \
			-derep 1 -out_format 3 -out_good ${folder_path}${SRR}_trimmed_dedup \
			-out_bad null

		echo "removing unneeded files for sample $SRR"
		rm ${folder_path}${SRR}_trimmed_1.fastq ${folder_path}${SRR}_trimmed_2.fastq

	elif [ "$seq_type" == "SE" ]; then
		echo "running fastp trimming and filtering for sample ${SRR}"
		fastp --trim_poly_g --trim_poly_x --length_required $length_filter --average_qual 30 \
			-i ${storage_path}${SRR}.fastq.gz \
			-o ${folder_path}${SRR}_trimmed.fastq.gz -V \
			-h ${folder_path}fastp.html -j ${folder_path}fastp.json

		echo "unzipping the trimmed and filtered files for sample ${SRR} for prinseq"
		gunzip ${folder_path}${SRR}_trimmed.fastq.gz

		echo "prinseq deduplication for sample ${SRR}"
		prinseq-lite.pl -verbose -fastq ${folder_path}${SRR}_trimmed.fastq \
			-derep 1 -out_format 3 -out_good ${folder_path}${SRR}_trimmed_dedup \
			-out_bad null

		echo "removing unneeded files for sample $SRR"
		rm ${folder_path}${SRR}_trimmed.fastq
	fi
	rm ${folder_path}fastp*
	echo "-----------------"



done
echo "complete"
echo "-----------------"


