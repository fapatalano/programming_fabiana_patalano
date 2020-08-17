def createCouple(aa):
    coupleList=[]
    for x in range(len(aa)):
        couples=[]
        for i in range (len(aa)):
            c=aa[x]+aa[i]
            couples.append(c)
        coupleList.append(couples)
    return coupleList

def collectScores(matrix,aminoacid):
    scores=[]
    dictionary={}
    for line in matrix:
        line=line.rstrip()
        line=line.split()
        for element in line:
            if element not in aminoacid:
                scores.append(element)
    return scores

def createDictionary(matrix, aminoacid):
    coupleList=createCouple(aminoacid)
    scores=collectScores(matrix, aminoacid)
    dictionary={}
    i=0
    for couples in coupleList:
        for couple in couples:
            dictionary[couple]=scores[i]
            i+=1
    print(dictionary)
    return dictionary



matrix=open("./blosum50.txt","r")
aminoacid="ARNDCQEGHILKMFPSTWYV"
createDictionary(matrix,aminoacid)
