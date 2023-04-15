import sys

list_args = list(sys.argv)

fasta = list_args[1]
output = list_args[2]

def SIXPACK_translate(fasta_file,output_file):
    def translate(seq):    
        table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
        protein =""
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
        return protein

    def cutter(seq):
        if len(seq) %3 == 0:
            return seq
        elif len(seq) %3 == 1:
            return seq[:-1]
        elif len(seq) %3 == 2:
            return seq[:-2]

    def complement_strand(seq):
        table1 = {
            'A':'T',
            'G':'C',
            'C':'G',
            'T':'A'
        }
        complement = ''
        for i in range(0,len(seq)):
            complement += table1[seq[i]]
        return complement[::-1]
    
    list_names = ['FW_ORF1','FW_ORF2','FW_ORF3','RV_ORF1','RV_ORF2','RV_ORF3'] 
    resulting_peptides = open(output_file,'w')

 
    with open(fasta_file,'r') as sequences:
        lines = sequences.readlines()
    i = 0
    while i < len(lines):
        if 'n' in lines[i+1] or 'N' in lines[i+1]:
            pass
            i += 2

        else:
            name_seq = {}

            name = lines[i].strip().replace('>','')
            list_names_individual = [name+'_'+element for element in list_names]

            i += 1
            sequence = lines[i].strip().replace('>','').upper()

            #if analysis == 'LINES':
            #    sequence = sequence_pre
            #elif analysis == 'FLANKS':
            #    sequence = sequence_pre[::-1]

            FW_seq1 = translate(cutter(sequence))
            FW_seq2 = translate(cutter(sequence[1:]))
            FW_seq3 = translate(cutter(sequence[2:]))

            complement1 = complement_strand(sequence)
            RV_seq1 = translate(cutter(complement1))
            RV_seq2 = translate(cutter(complement1[1:]))
            RV_seq3 = translate(cutter(complement1[2:]))
            list_seq = [FW_seq1,FW_seq2,FW_seq3,RV_seq1,RV_seq2,RV_seq3]

            for number,name_seq in enumerate(list_names_individual):
                resulting_peptides.write(f'>{name_seq}\n')
                resulting_peptides.write(f'{list_seq[number]}\n')

            i+=1
    return print('Translation in three frames COMPLETE')

SIXPACK_translate(fasta,output)