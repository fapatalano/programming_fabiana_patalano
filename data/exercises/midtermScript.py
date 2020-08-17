''' write a function that reads the TWO sequences from the file sequences.fasta into Python strings and put them in two variables seq1,seq2, respectively.
    The function should return a touple (seq1,seq2).
    write a function that reads a table from the text file BLOSUM50.txt and puts its value in a dictionary B, where pairs of aa letters are the keys.
    The function should return a dictionary B.'''


def makeDictionary(file,aminoacids):
    B={}
    for aminoacid in aminoacids:
        row=file.readline().split()
        for i in range(len(row)):
            couple=aminoacid+aminoacids[i]
            B[couple]= float(row[i])
    return B
def calcAlScore(matrix,seq1,seq2,d):
    score=0
    for i in range(len(seq1)):
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
    for j in matrix:
        print(seq1)
        print(seq2)
        print ("The score is: ", calcAlScore(j,seq1,seq2,d))

def readLines(file):
    lines=file.readlines()
    print(lines)
    for i in range(1,len(lines)/2,2):
        print(i)
    seq1=str(lines[1].rstrip())
    seq2=str(lines[3].rstrip())
    return(seq1,seq2)

aminoacid="ARNDCQEGHILKMFPSTWYV"
pam=open("./PAM250.txt","r")
blosum=open("./blosum62.txt","r")
fasta_file=open("./alignments.fasta","r")
#fasta_file=open("./2seq.fasta","r")
d=2
scoringAl(pam,blosum,aminoacid,d,fasta_file)
