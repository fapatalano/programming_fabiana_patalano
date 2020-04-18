def readFile(fasta_file):

    ''' This function that reads two amino acid sequences into Python strings from a file in FASTA format

    This function take as input the file and return as output a tuple (seq1,seq2)
    '''

    next(fasta_file)              #remove the first header
    lines=fasta_file.readlines()  #list containing all the lines
    sequence=''
    sequences=[]
    for line in lines:
        line=line.strip()
        if line.startswith('>'):
            sequences.append(sequence)
            sequence=''
        else:
            sequence+=line
    sequences.append(sequence)
    seq1,seq2=sequences[0],sequences[1]
    return (seq1,seq2)

def makeDictionary(file,aminoacids):

    ''' The function reads a substitution matrix from a file and puts value in a dictionary

    Takes in input the file and the amino acids letters. Pairs of  aminoa acid letters are
    the keys of the dictionary. The functio return the dictionary
    '''

    matrix={}               # initialize an empty dictionary
    for aminoacid in aminoacids:
        row=file.readline().split() # splits each line into a list
        for i in range(len(row)):
            couple=aminoacid+aminoacids[i]
            matrix[couple]= float(row[i])           #add in the dictionary also the reverse couple
            matrix[couple[::-1]]=float(row[i])      # as the matrix is symmetric
    return matrix

def calcAlignment(seq1,seq2,dictionary,gp=0):

    '''The function finds the best scoring alignment implementing an ungapped exhaustive alignment algorithm.

    The function takes as input the sequences and the scoring matrix; gap penalty is equal to 0. The
    function will return the best alignment(s) and the score.
    '''

    if len(seq1)<len(seq2):              #a number of dash equal to the lenght of the longest
        long=list(seq2+ '-'*len(seq1))   #sequence are added at the beginning of the shortest one
        short=list('-'*len(seq2) +seq1)  #viceversa a number of dash equal to the shortest are
    else:                                #added at the end of the longest
        long=list(seq1+ '-'*len(seq2))
        short=list('-'*len(seq1) +seq2)
    score=gp*len(short)
    bestAl=[''.join(short),''.join(long),score] #initialize a list where to store the best alignment(s).

    # Iteration: sequeces slides on each other
    for i in range((len(seq1)+ len(seq2))-1):   #the number of iteration is equal to the sum of the lenght of both
        score=0                                 #sequences because the first and the last alignments are equal
        if long[-1]=='-' and short[0]=='-':
            del short[0]                          #delete dash at the beginning
            del long[-1]                          #delete dash at the end
            for a1, a2 in zip(short,long):        #get the score of each of the amino acid aligned
                if a1 == '-' or a2 == '-': score += gp
                else: score += dictionary[a1 + a2]
        elif long[-1] != '-' and short[0] == '-': #if the long sequence has no more dash at the end but the short has
            short+='-'                            #append dash at the end
            del short[0]                          #delete dash at the beginning
            for a1, a2 in zip(short,long):
                if a1 == '-' or a2 == '-': score += gp
                else: score += dictionary[a1 + a2]
        else:
            long.insert(0, '-')                   #append dash at the beginning
            short+='-'
            for a1, a2 in zip(short,long):
                if a1 == '-' or a2 == '-': score += gp
                else: score += dictionary[a1 + a2]
        if score > bestAl[2]: bestAl=[''.join(short),''.join(long),score]   #if the score is better replace the sequences
        elif score == bestAl[2]: bestAl+=[''.join(short),''.join(long),score]   #if the score is equal append the sequences
    return bestAl

fasta_file=open("./fasta_file.fasta","r")
seq1,seq2=readFile(fasta_file)
# seq1="HAGSGK"
# seq2="AGKSHAAAA"
aminoacids="ARNDCQEGHILKMFPSTWYV"
blosum=open("./blosum62.txt","r")
dictionary= makeDictionary(blosum,aminoacids)
bestAl=calcAlignment(seq1,seq2,dictionary)
print('\n'.join(map(str,bestAl)))   #map applies str function to each item of bestAl
