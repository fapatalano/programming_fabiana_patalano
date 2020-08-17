def readFile(fasta_file):

    ''' This function takes as input a FASTA file and return
    as output a tuple containing the first two sequences '''

    lines=fasta_file.readlines()          #list containing all the lines of file
    sequence=''
    sequences=[]
    for i in range(1,len(lines)):
        lines[i]=lines[i].rstrip()
        if lines[i].startswith('>'):      # if the line starts with '>' sequence is firstly added
            sequences.append(sequence)    # to sequences and then is assigned again an empty string
            sequence=''
        else:
            sequence+=lines[i]            # append the line to sequence
    sequences.append(sequence)
    seq1,seq2=sequences[0],sequences[1]
    return (seq1,seq2)

def makeDictionary(file,aminoacids):

    ''' The function takes in input the file with a substitution matrix and the amino acids letters.
    Pairs of amino acid letters are used as keys of the dictionary. The function will return the dictionary
    '''

    matrix={}                                     # initialize an empty dictionary
    for aa in aminoacids:
        row=file.readline().split()               # splits each line into a list
        for i in range(len(row)):
            couple=aa+aminoacids[i]
            matrix[couple]= float(row[i])         # add in the dictionary also the reverse couple
            matrix[couple[::-1]]=float(row[i])
    return matrix

def getScore(aa1,aa2,gp):
    if aa1 == '-' or aa2 == '-': score = gp
    else: score = dictionary[aa1 + aa2]
    return score

def calcAlignment(seq1,seq2,dictionary,gp=0):

    '''The function finds the best scoring alignment implementing an ungapped exhaustive alignment algorithm.
    Takes as input the sequences and the scoring matrix and will return the best alignment(s) and the score.
    Non-aligned residues at the beginning or at the end will score 0
    '''
    if len(seq1)<len(seq2):              #a number of dash equal to the lenght of the longest
        long=list(seq2+ '-'*len(seq1))   #sequence are added at the beginning of the shortest one
        short=list('-'*len(seq2) +seq1)  #viceversa a number of dash equal to the lenght of the shortest
    else:                                #sequence are added at the end of the longest one
        long=list(seq1+ '-'*len(seq2))
        short=list('-'*len(seq1) +seq2)
    score=gp*len(short)
    result=[''.join(short),''.join(long),score]   #initialize a list where to store the best alignment(s).
    # Iteration: sequeces slides on each other
    for i in range((len(seq1)+len(seq2))-1):      #the number of iteration is equal to the sum of the lenght of both
        score=0                                   #sequences minus 1 because the first and the last alignments are equal
        if long[-1]=='-' and short[0]=='-':
            del short[0]                          #delete dash at the beginning
            del long[-1]                          #delete dash at the end
            for aa1, aa2 in zip(short,long):      #get the score of each of the amino acid aligned
                score+=getScore(aa1,aa2,gp)
        elif long[-1] != '-' and short[0] == '-': #if the long sequence has no more dash at the end but the short has
            short+='-'                            #append dash at the end
            del short[0]                          #delete dash at the beginning
            for aa1, aa2 in zip(short,long):
                score+=getScore(aa1,aa2,gp)
        else:
            long.insert(0, '-')                   #append dash at the beginning
            short+='-'                            #append dash at the end
            for aa1, aa2 in zip(short,long):
                score+=getScore(aa1,aa2,gp)
        if score > result[2]: result=[''.join(short),''.join(long),score]       #if the score is better replace the sequences
        elif score == result[2]: result+=[''.join(short),''.join(long),score]   #if the score is equal append the sequences
    return result

fasta_file=open("./fasta_file.fasta","r")
# fasta_file=open("./2seq.fasta","r")
seq1,seq2=readFile(fasta_file)
aminoacids="ARNDCQEGHILKMFPSTWYV"
blosum=open("./blosum62.txt","r")
dictionary= makeDictionary(blosum,aminoacids)
result=calcAlignment(seq1,seq2,dictionary)
print('\n'.join(map(str,result)))   #map applies str function to each item of result
