#!/usr/bin/python3

import sys
from Bio import SeqIO


'''
create classes for this utils function.. this can be used in different objects.

'''




 def extract_id(path):

    id_list=[]
    for id in sorted(os.listdir(path)):
        id_list.append(id_list.split('.p')[0])

    return id_list

def to_df(model,w):

    d=w//2
    gor_columns=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    ss = ['H', 'E', 'C','R']
    window = range(-d,d+1)
    gor_df = pd.DataFrame(model,pd.MultiIndex.from_product([ss,window],names=['R_SS', 'PW']),gor_columns)

    return gor_df



def print_model(model,dir):
    model.to_csv(path_or_buf=dir+"information_matrix_df.tsv", sep='\t')
    #model.to_csv(path_or_buf=dir+"gor_df.tsv", sep='\t')
    #np.savetxt(dir+"gor_model.tsv", model, delimiter="\t")

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
        sequence=str(seq_record.seq)
        if len(sequence)<50 or len(sequence)>300 :    print(seq_record.id)


for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
        sequence=str(seq_record.seq)
        if sequence.isupper()==False:    print(sequence)
def a function:
    wanted = [line.strip() for line in open(sys.argv[2])]
    seqiter = SeqIO.parse(open(sys.argv[1]), 'fasta')
    for seq in seqiter:
        if seq.id[:6] in wanted:
            print('>'+seq.id[:6])
            print(seq.seq)

for record in SeqIO.parse(sys.argv[1], "fasta"):
    if record.seq.count('X') == 0:
        print(record.format("fasta"))
