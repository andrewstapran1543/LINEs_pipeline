# LINEs_pipeline
## This is a description of how a new pipeline was applied. The pipeline approaches finding new TE-derived antigens from two sides:
## 1) First approach - align the RNAseq reads to full-size genome LINE copies (RepeatMasker). Then identify continuous overexpressed regions longer than some threshold (in my case 150) - based on pair comparisons between matched norm-tumor samples (or median/mean values). Further translate the overexpressed fasta sequences into 6 peptide sequences using python notebook on local computer:
```
/Users/andrewstapran/Desktop/LINEs_new/python_notebooks/cutter_amino_acids.ipynb
```
## Further on we will map the peptides to our translatome:
```
/Users/andrewstapran/Desktop/LINEs_new/python_notebooks/cutter_amino_acids.ipynb
```
## In this way we get the list of peptides that are expressed on MHC and map to the overexpressed regions of LINE copies based on RNAseq data
##
## 2) Second approach - translate the full size LINE copies into 6 ORFs & map the peptides from MHC to translatome. Go back to the nucleotide positions and fasta sequences of where the peptides mapped. Finally, align the RNAseq reads to these fasta sequences - identify the overexpressed regions. Do they coincide with the footprints of where the peptides mapped? If yes - the protein is indeed expressed as MHC
## 
## Finally, get two lists of peptides - overlap them
