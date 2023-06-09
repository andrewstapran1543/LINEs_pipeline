### LINEs_pipeline
This repository describes the analysis pipeline for finding new retroelement-based cancer specific antigens<br/>
Below is the outline of how the analysis proceeds, and in the folders you can find:
<ul>
  <li>Reference fasta and bed files for LINE retroelements</li>
  <li>Reference MHC-immunopeptidome with description of datasets used for antigen search</li>
  <li>Scripts that were used for the different steps of the analysis</li>
</ul>
Below is the outline of the analysis pipeline:
<img src=https://github.com/andrewstapran1543/LINEs_pipeline/blob/main/analysis_outline.png width="600">
Folder <a href="https://github.com/andrewstapran1543/LINEs_pipeline/tree/main/SCRIPTS">SCRIPTS</a> in this repository has the following scripts:<br/>
<ul>
<li><b>trimmer_filter.sh</b> - this script trims and filter the raw reads. Takes in the following parameters: 1) samples list, 2) seq_type - PE or SE, 3) experiment_id  - how the folder with analysis files will be called, 4) length filter - lower threshold for filtering the reads by length using fastp</li>
<li><b>aligner.sh</b> - this script aligns the processed reads to 1) LINEs, 2) LINE FLANKS, 3) LINE peptide footprints, 4) LINE FLANKS peptide footprints. Takes in the following parameters: 1) samples list, 2) seq_type - PE or SE, 3) experiment_id  - how the folder with analysis files is called, 4) ANALYSIS TYPE - FULL_ONLY (for aligning to LINE copies only), FLANKS_ONLY (for aligning to LINE FLANKS only), FLANKS_FULL (for aligning to both LINE FLANKS and LINE copies), ALL (similar to previous but also aligning to genome for TETranscript). The reference files can be found <a href="https://github.com/andrewstapran1543/LINEs_pipeline/tree/main/REFERENCES">here</a></li>
<li><b>translator.py</b> - this script takes in the fasta file and translates it in 6 possible frames into peptide sequences. Takes in the following parameters: 1) fasta file in nucleotide format, 2) location of the output peptide sequence file</li>
<li><b>integrator.py</b> - this script combines resulting coverage file into a single output final dataframe. Takes in the following parameters: 1) location of the input read number file (for each sample), 2) experiment_id (how the folder with analysis files is called), 3) analysis type (FLANKS or FULL) - for alignment to FULL LINE copies or FLANK regions, 4) mode (REG_PEP or FOOTPR) - for original sequences or peptide mapping footprint sequences</li>
<li><b>calculator.py</b> - this script calculates continuous regions overexpressed in tumor compared to normal tissues based on final dataframe generated by previous script. Takes in the following parameters: 1) experiment_id (how the folder with analysis files is called), 2) analysis type (FLANKS or FULL) - for alignment to FULL LINE copies or FLANK regions, 3) location of the file where tumor-normal samples are decribed, 4) mode (REG_PEP or FOOTPR) - for original sequences or peptide mapping footprint sequences, 5) threshold for how many pairs out of total tumor-normal pairs have considered position overexpressed in tumor compared to normal tissue, 6) threshold for how long an overexpressed region should be (in nucleotides)</li>
<li><b>peptide_mapper.py</b> - this script takes in the continuous regions file, translates them and maps the peptides from <a href="https://github.com/andrewstapran1543/LINEs_pipeline/blob/main/IMMUNOPEPTIDOME_DATASETS/united_peptidome.txt">immunopeptidome list</a> - returns table with mapping details and the list of peptides. Takes in the following parameters: 1) experiment_id (how the folder with analysis files is called), 2) analysis type (FLANKS or FULL) - for alignment to FULL LINE copies or FLANK regions, 3) mode (REG_PEP or FOOTPR) - for original sequences or peptide mapping footprint sequences, 4) MAPPING MODE for blast - RELAXED or STRINGENT</li>
<li><b>cont_assembler.py</b> - this script combines three previous scripts together - integrator, calculator, peptide_mapper. Takes in the following parameters: 1) location of the input read number file (for each sample), 2) location of the file where tumor-normal samples are decribed, 3) experiment_id (how the folder with analysis files is called), 4) threshold for how many pairs out of total tumor-normal pairs have considered position overexpressed in tumor compared to normal tissue, 5) threshold for how long an overexpressed region should be (in nucleotides), 6) MAPPING MODE for blast - RELAXED or STRINGENT</li>
</ul>
The preliminary results for four cancer types (colorectal, renal, kidney, and prostate) can be found <a href="https://github.com/andrewstapran1543/LINEs_pipeline/tree/main/RESULTS">here</a><br/><br/>
Here is the list of RNA-seq datasets used for the analysis:
<ul>
  <li><a href="https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=149603">PROSTATE cancer</a> - PE sequencing</li>
  <li><a href="https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=24446862">KIDNEY cancer</a> - PE sequencing</li>
  <li><a href="https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=39868">HBV-related hepatocellualr carcinoma</a> - SE sequencing</li>
  <li>Colorectal cancer - dataset from the laboratory</li>
</ul>
Trimmed reads for all datasets can be found at the following address on the server: <i>/ngs/2023/02_rnaseq_datasets_Stapran/</i><br/><br/>
<b><i>IMPORTANTLY!!!</i></b><br/>Change the paths to folders in the <a href="https://github.com/andrewstapran1543/LINEs_pipeline/tree/main/SCRIPTS">scripts</a> before running them
