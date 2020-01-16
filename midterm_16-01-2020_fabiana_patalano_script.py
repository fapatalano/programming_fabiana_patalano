def createCouple(aminoacid):
    couples=[]
    for x in range(len(aminoacid)):
        for i in range (len(aminoacid)):
            c=aminoacid[x]+aminoacid[i]
            if c[::-1] not in couples:
                couples.append(c)
    return couples

def collectScores(matrix,aminoacid):
    lines=[]
    for line in matrix:
        line=line.rstrip()
        line=line.split()
        lines.append(line)
    x=20
    scores=[]
    for n in range(len(aminoacid)):
        columns=[]
        for j in range(x):
            columns.append(lines[j][n])
        x-=1
        scores.extend(columns)
        del lines[0]
    return scores

def createDictionary(file,aminoacid):
    coupleList=createCouple(aminoacid)
    scores=collectScores(file,aminoacid)
    dictionary={}
    i=0
    for couple in coupleList:
        dictionary[couple]=scores[i]
        i+=1
    return dictionary

def calcAlScore(matrix,seq1,seq2,d):
    score=0
    for i in range(len(seq1)):
        pair=seq1[i]+seq2[i]
        if "-" in pair:
            score-=d
        elif pair[::-1] in matrix:
            value=matrix[pair[::-1]]
            score+=int(value)
        else:
            value=matrix[pair]
            score+=int(value)
    return score

def scoringAl(matrix1,matrix2,sequences,aminoacid,d):
    matrices=[]
    dict1=createDictionary(matrix1,aminoacid)
    dict2=createDictionary(matrix2,aminoacid)
    matrices.append(dict1)
    matrices.append(dict2)
    for j in range(len(sequences)):
        for matrix in matrices:
            print(seq[j],seq[j+1])
            print ("The score is: ", calcAlScore(matrix,sequences[j], sequences[j+1],d))

aminoacid="ARNDCQEGHILKMFPSTWYV"
pam=open("./PAM250.txt","r")
blosum=open("./blosum62.txt","r")
fasta_file=open("./alignments.fasta","r")
d=2

def readFile(fasta_file):
    sequences = []
    for line in fasta_file:
        if not line.startswith('>'):
            sequences.append(line.strip())
    return sequences

seq=readFile(fasta_file)
scoringAl(pam,blosum,seq,aminoacid,d)
