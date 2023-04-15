#!/bin/bash
#input parameters
list_samples=$1
seq_type=$2
experiment_id=$3
analysis=$4



#input samples
samples_total=$(cat $list_samples)
samples_number=$(cat $list_samples | wc -l)



#fixed parameters
echo "making folder for the experiment"
mkdir -p /home/andrewstap/LINEs/Data_Analysis/${experiment_id}/
echo "setting reference indices"
storage_path="/ngs/2023/02_rnaseq_datasets_Stapran/${experiment_id}/ANALYSIS/"
full_size_ftpr_index="/home/andrewstap/LINEs/REFERENCES/GENOME/GCHR38/RMSK/MAPPINGS/LINES_200/LINES_FOOTPRINTS"
full_size_copy_index="/home/andrewstap/LINEs/REFERENCES/GENOME/GCHR38/RMSK/LINES_index"
flanks_ftpr_index="/home/andrewstap/LINEs/REFERENCES/GENOME/GCHR38/RMSK/MAPPINGS/FLANKS_200/FLANKS_FOOTPRINTS"
flanks_full_index="/home/andrewstap/LINEs/REFERENCES/GENOME/GCHR38/RMSK/FLANKS/UNMASKED/FLANKS_FULL"
genome_index="/home/andrewstap/LINEs/REFERENCES/GENOME/STAR_INDEX"




#setting functions
echo "setting functions"
#bowtie2 alignment
bowtie2_aligner () {
	align_pathway=$1
	reference_index=$2
	file_name=$3

	echo "running $file_name alignment for $seq_type file type"
	case "$seq_type" in "SE" ) \
		bowtie2 -q --threads 20 -a -x $reference_index -U ${storage_path}${SRR}_trimmed_dedup/${SRR}_trimmed_dedup.fastq \
			> ${align_pathway}${SRR}_${file_name}.sam 2> ${align_pathway}${SRR}_align_summary.txt;;
	"PE" ) \
		bowtie2 -q --threads 20 --no-mixed -a -x $reference_index -1 ${storage_path}${SRR}_trimmed_dedup/${SRR}_trimmed_dedup_1.fastq \
	        -2 ${storage_path}${SRR}_trimmed_dedup/${SRR}_trimmed_dedup_2.fastq \
	        > ${align_pathway}${SRR}_${file_name}.sam 2> ${align_pathway}${SRR}_align_summary.txt;;
	* ) \
		echo "unknown sequencing type OR sequencing type not specified";;
	esac	
}

#samtools interpreter of the sam files
samtools_subsetter () {
	align_pathway=$1
	file_name=$2

	samtools view -h -@ 20 -F 0x04 ${align_pathway}${SRR}_${file_name}.sam > ${align_pathway}${SRR}_${file_name}_mapped.sam
	touch ${align_pathway}${SRR}_${file_name}_final.sam
	grep '@' ${align_pathway}${SRR}_${file_name}_mapped.sam >> ${align_pathway}${SRR}_${file_name}_final.sam
	grep '\<NM:i:[012]\>' ${align_pathway}${SRR}_${file_name}_mapped.sam >> ${align_pathway}${SRR}_${file_name}_final.sam
	samtools sort ${align_pathway}${SRR}_${file_name}_final.sam > ${align_pathway}${SRR}_${file_name}_final_sorted.sam

	echo "mapped stats"
	samtools stats -@ 20 ${align_pathway}${SRR}_${file_name}_mapped.sam | grep ^SN | cut -f 2- > ${align_pathway}${SRR}_${file_name}_stats_mapped.txt
	result=$(expr $(wc -l ${align_pathway}${SRR}_${file_name}_mapped.sam | cut -d' ' -f1) - 8778)
	echo -e "total number of alignments \t $result" >> ${align_pathway}${SRR}_${file_name}_stats_mapped.txt

	echo "NM012 stats"
	samtools stats -@ 20 ${align_pathway}${SRR}_${file_name}_final_sorted.sam | grep ^SN | cut -f 2- > ${align_pathway}${SRR}_${file_name}_stats.txt
	result=$(expr $(wc -l ${align_pathway}${SRR}_${file_name}_final_sorted.sam | cut -d' ' -f1) - 8778)
	echo -e "total number of alignments \t $result" >> ${align_pathway}${SRR}_${file_name}_stats.txt

	echo "calculating coverage"
	samtools depth ${align_pathway}${SRR}_${file_name}_final_sorted.sam > ${align_pathway}${SRR}_${file_name}_coverage.txt
	rm ${align_pathway}*sam	
}

#genome aligner
genome_aligner () {
	output_folder=$1
	sample_name=$2

	case "$seq_type" in "SE" ) \	
	STAR --genomeDir $genome_index \
	--runThreadN 20 \
	--readFilesIn ${storage_path}${sample_name}_trimmed_dedup/${sample_name}_trimmed_dedup.fastq \
	--outSAMtype BAM Unsorted \
	--quantMode GeneCounts \
	--outFileNamePrefix ${output_folder}${sample_name}_ \
	--outFilterMultimapNmax 100 \
	--winAnchorMultimapNmax 100;;
	
	"PE" ) \
	STAR --genomeDir $genome_index \
	--runThreadN 20 \
	--readFilesIn ${storage_path}${sample_name}_trimmed_dedup/${sample_name}_trimmed_dedup_1.fastq \
		${storage_path}${sample_name}_trimmed_dedup/${sample_name}_trimmed_dedup_2.fastq \
	--outSAMtype BAM Unsorted \
	--quantMode GeneCounts \
	--outFileNamePrefix ${output_folder}${sample_name}_ \
	--outFilterMultimapNmax 100 \
	--winAnchorMultimapNmax 100;;

	* ) \
		echo "unknown sequencing type OR sequencing type not specified";;
	esac

}



#main part
for i in $(seq 1 $samples_number); do
	SRR=$(echo $samples_total | cut -d' ' -f$i)
	echo "alignment for $SRR"
	echo "making directories"
	experiment_path="/home/andrewstap/LINEs/Data_Analysis/${experiment_id}/$SRR/"
	full_align_path="${experiment_path}FULL/regions_to_peptides/"
	full_align_path_ftpr="${experiment_path}FULL/peptides_to_regions/"
	flanks_align_path="${experiment_path}FLANKS/regions_to_peptides/"
	flanks_align_path_ftpr="${experiment_path}FLANKS/peptides_to_regions/"
	align_genome="${experiment_path}STAR_GENOME/"
	mkdir -p $experiment_path
	mkdir -p $align_genome
	mkdir -p $full_align_path
	mkdir -p $full_align_path_ftpr
	mkdir -p $flanks_align_path
	mkdir -p $flanks_align_path_ftpr

	fileinsert1="reg_to_pep"
	fileinsert2="pep_to_reg"

	case "$analysis" in "GENOME_ONLY" )
	genome_aligner $align_genome $SRR;;
	"FULL_ONLY" )
	bowtie2_aligner $full_align_path $full_size_copy_index $fileinsert1
	samtools_subsetter $full_align_path $fileinsert1
	bowtie2_aligner $full_align_path_ftpr $full_size_ftpr_index $fileinsert2
	samtools_subsetter $full_align_path_ftpr $fileinsert2;;
	"FLANKS_ONLY" )
	bowtie2_aligner $flanks_align_path $flanks_full_index $fileinsert1
	samtools_subsetter $flanks_align_path $fileinsert1
	bowtie2_aligner $flanks_align_path_ftpr $flanks_ftpr_index $fileinsert2
	samtools_subsetter $flanks_align_path_ftpr $fileinsert2;;
	"FLANKS_FULL" )
	bowtie2_aligner $full_align_path $full_size_copy_index $fileinsert1
	samtools_subsetter $full_align_path $fileinsert1
	bowtie2_aligner $full_align_path_ftpr $full_size_ftpr_index $fileinsert2
	samtools_subsetter $full_align_path_ftpr $fileinsert2
	bowtie2_aligner $flanks_align_path $flanks_full_index $fileinsert1
	samtools_subsetter $flanks_align_path $fileinsert1
	bowtie2_aligner $flanks_align_path_ftpr $flanks_ftpr_index $fileinsert2
	samtools_subsetter $flanks_align_path_ftpr $fileinsert2;;
	"ALL" )
	genome_aligner $align_genome $SRR
	bowtie2_aligner $full_align_path $full_size_copy_index $fileinsert1
	samtools_subsetter $full_align_path $fileinsert1
	bowtie2_aligner $full_align_path_ftpr $full_size_ftpr_index $fileinsert2
	samtools_subsetter $full_align_path_ftpr $fileinsert2
	bowtie2_aligner $flanks_align_path $flanks_full_index $fileinsert1
	samtools_subsetter $flanks_align_path $fileinsert1
	bowtie2_aligner $flanks_align_path_ftpr $flanks_ftpr_index $fileinsert2
	samtools_subsetter $flanks_align_path_ftpr $fileinsert2;;
	* )
	echo "Unknown task type";;
	esac	

	echo "---------------------"
done

echo "COMPLETE"
echo "---------------------"


