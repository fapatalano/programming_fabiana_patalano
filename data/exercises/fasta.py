file=open("./blosum.txt","r")
def couple(bases):
    couples=[]
    for x in range(len(bases)):
        for i in range(len(bases)):
            c=bases[x]+bases[i]
            couples.append(c)
    return couples

def make_list_file(file):
    lines=[]
    for line in file:
        line=line.rstrip()
        line=line.split()
        sline=[]
        for el in line:
            sline.append(el)
        lines.append(sline)
    del lines[0]
    for i in range(len(lines)):
        del lines[i][0]
    lfinal=[]
    for i in range(len(lines)):
        for l in range(len(lines)):
         lfinal.append(lines[i][l])
    return lfinal


def dictionary(file, ammino):
    c=couple(ammino)
    #list=[]
    listel=make_list_file(file)
    diz={}
    i=0
    for s in c:
        diz[s]=listel[i]
        i+=1
    return diz


ammino="ARNDCQEGHILKMFPSTWYV"

blosum50=dictionary(file, ammino)

file.close()


def scoring(seq1, seq2):
    s = 0
    for i in range(len(seq1)):
        blosum_key = seq1[i]+seq2[i]
        if blosum50.key(blosum_key):
            s = s + int(blosum50[blosum_key])

    return s

file1=open("./hba_hbb_human.fasta")

def read_input(fasta_file):
    seqs = []
    for line in fasta_file:
        if not line.startswith('>'):
            seqs.append(line.strip())
    return seqs

seqs=read_input(file1)
print seqs

if len(seqs) == 2:
    print "the scoring of the alignment is: ", scoring(seqs[0], seqs[1])
else:
    print "too many sequences in the array"

file1.close()
