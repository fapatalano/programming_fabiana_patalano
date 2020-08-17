def createCouple(aa):
    couples=[]
    for x in range(len(aa)):
        for i in range (len(aa)):
            c=aa[x]+aa[i]
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

def createDictionary(matrix,aminoacid):
    coupleList=createCouple(aminoacid)
    scores=collectScores(matrix,aminoacid)
    dictionary={}
    i=0
    for couple in coupleList:
        dictionary[couple]=scores[i]
        i+=1
    print ("VV",dictionary['VV'])
    return dictionary

matrix=open("./PAM250.txt","r")
aminoacid="ARNDCQEGHILKMFPSTWYV"
createDictionary(matrix,aminoacid)
