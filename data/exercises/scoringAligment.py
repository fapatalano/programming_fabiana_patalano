
''' write a function that reads the two sequences from the file.fasta into Python strings and put them in two variables seq1,seq2, respectively.
The function should return a touple (seq1,seq2). write a function that reads a table from the text file BLOSUM50.txt and puts its value in a dictionary
B where pairs of aa letters are the keys. the function sholud return a dictionary B.'''

def makeDictionary(file,aminoacids):
    matrix={}
    for aminoacid in aminoacids:
        row=file.readline().split()
        for i in range(len(row)):
            couple=aminoacid+aminoacids[i]
            matrix[couple]= float(row[i])
    return matrix

def calcAlScore(matrix,seq1,seq2,d):
    score=0
    for i in range(len(seq1)):
        print('seq1',seq1[i])
        pair=seq1[i]+seq2[i]
        if "-" in pair:
            score-=d
        elif pair[::-1] in matrix:
            value=matrix[pair[::-1]]
            score+=float(value)
        else:
            value=matrix[pair]
            score+=float(value)
    return score

def scoringAl(matrix1,matrix2,aminoacid,d,file):
    matrix=[]
    seq1,seq2=readLines(file)
    pam=makeDictionary(matrix1,aminoacid)
    blosum=makeDictionary(matrix2,aminoacid)
    matrix.append(pam)
    matrix.append(blosum)
    for i in range(len(seq1)):
        for j in matrix:
            print ("The score is: ", calcAlScore(j,seq1[i], seq2[i],d))

'''write a function that reads the sequence fasta into python string and put them in the vairiables seq1,seq2.
the function should return a tuple (seq1,seq2)'''

# def readFile(fasta_file):
#     sequences = []
#     for line in fasta_file:
#         if not line.startswith('>'):
#             sequences.append(line.rstrip())
#     print(sequences)
#     return sequences

def readLines(file):
    lines=file.readlines()
    seq1=[]
    seq2=[]
    for i in range(1,len(lines)-1,4):
        seq1.append(lines[i].rstrip())
        seq2.append(lines[i+2].rstrip())
    return(seq1,seq2)


aminoacid="ARNDCQEGHILKMFPSTWYV"
pam=open("./PAM250.txt","r")
fasta_file=open("./alignments.fasta","r")
#fasta_file=open("./2seq.fasta","r")
d=2
scoringAl(pam,blosum,aminoacid,d,fasta_file)
#makeDictionary(blosum,aminoacid)
