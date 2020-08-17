def calcZeroMatrix(seq1,seq2,d,subMatrix):
    matrix=[]
    rowScore=0
    for i in range(len(seq1)+1):
        score=d
        row=[]
        for j in range(len(seq2)+1):
            if i== 0:                   # fill row
                score-=d
                row.append(score)
            elif j==0:                 # fill columns
                rowScore-=d
                row.append(rowScore)
            else:
                s=0
                row.append(s)
        matrix.append(row)
    return matrix

def matrixScore(a,b,matrix):
    key=a+b
    if key in matrix or key[::-1] in matrix:
        return matrix[key]
    else:
        print("ERROR",key)

def calcMax(r,c,score,gapPenalty,scoresMatrix):
    scores=[]
    d=float(scoresMatrix[r-1][c-1])+float(score)
    u=float(scoresMatrix[r-1][c])-gapPenalty
    l=float(scoresMatrix[r][c-1])-gapPenalty
    if d >= u and d>=l:
        scores.append(d)
        scores.append('d')
    elif u>l and u>d :
        scores.append(u)
        scores.append('u')
    elif l>u and l>d:
        scores.append(l)
        scores.append('l')
    print(scores)
    return scores

def fillMatrix(seq1,seq2,d,subMatrix):
    tracebackMatrix=[]
    scoreMatrix=calcZeroMatrix(seq1,seq2,d,subMatrix)
    for r in range(1,len(seq1)+1):       #row
        tracebackList=[]
        for c in range(1,len(seq2)+1):   #columns
            score=matrixScore(seq1[r-1],seq2[c-1],subMatrix)
            max=calcMax(r,c,score,d,scoreMatrix)
            scoreMatrix[r][c]=max[0]
            tracebackList.extend(max[1])
        tracebackMatrix.append(tracebackList)
    print("list",tracebackMatrix)
    return tracebackMatrix

def calcAlignment(seq1,seq2,d,scoringMatrix):
    tracebackMatrix=fillMatrix(seq1,seq2,d,scoringMatrix)
    seqAligned1= ''
    seqAligned2 = ''
    ls1=len(seq1)-1
    ls2=len(seq2)-1
#    for i in range(len(tracebackMatrix)-1,-1,-1):
    while ls1 != 0 or ls2 != 0:
        if tracebackMatrix[ls1][ls2]== 'd':
            seqAligned1+= seq1[ls1]
            seqAligned2+= seq2[ls2]
            ls1-=1
            ls2-=1
        elif tracebackMatrix[ls1][ls2]== 'u':
            seqAligned2+= '-'
            seqAligned1+= seq1[ls1]
            ls1-=1
        elif tracebackMatrix[ls1][ls2]== 'l':
            seqAligned1+= '-'
            seqAligned2+= seq2[ls2]
            ls2-=1
    seqAligned1 = seqAligned1[::-1]
    print(seqAligned1)
    seqAligned2 = seqAligned2[::-1]
    print(seqAligned2)
    return(seqAligned1, seqAligned2)

seq2=list("ACAATTTTTTTT")
seq1=list("ACTGA")
#seq1 =list("AGAGTCAGTCAGTCAGTCAGCTAGCAGCACGTA")
#seq2 =list("AAAAAAM")
d=2

matrix={"AA":2,"AT":-1,"AC":-1,"AG":0,
"TA":-1,"TT":2,"TC":0,"TG":-1,
"CA":-1,"CT":0,"CC":2,"CG":-1,
"GA":0,"GT":-1,"GC":-1,"GG":2}
calcAlignment(seq1,seq2,d,matrix)
#fillMatrix(seq1,seq2,d,matrix)
#calcZeroMatrix(seq1,seq2,d)
