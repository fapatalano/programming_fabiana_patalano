'''
Local Alignment: Smith & Waterman algorithm

Given two sequences,  seq1 and seq2, with lenghts m and n, define a (m+1)(n+1) matrix F(i.j)
storing the best local alignment between sequeces.

PSEUDOCODE:
compileMatrices(seq1,seq2,matrix,gp):
    INPUT:
        sequences: seq1,seq2
        substitution matrix: matrix
        gap penalty:gp
    OUTPUT:
        scoring matrix: F
        traceback matrix: P
        max value index(es):index
- initialize F and P zeros matrix
- double for loop for Smith & Waterman algorithm
- find the max score and fill P and F matrices:
                |   F(r-1,c-1) + matrix(Ar,Bc) ; 'd' in P matrix
  F(r,c) = Max  |   F(r-1,c) -gp                 'u' in P matrix
                |   F(r,c-1) -gp                 'l' in P matrix
                |   0
- store the maximal score and its indexe(s)
- return F,P,index

calcAlignment(F,P,index,seq1,seq2):
    INPUT:
        scoring matrix: F
        traceback matrix: P
        max value index(es):index
    OUTPUT:
        a list containing aligned sequence(s) and their score:
        result= seqAligned1,seqAligned2,score
- initialize an empty list where to append aligned sequences and the score
- for loop that iterates the list of indexes
- while loop that iterates until a 0 in F matrix is found
- Backtracking starts from the highest value cell in P and proceed following the path:
    d= the two are aligned
    u= introduce a gap in the top sequence
    l= introduce a gap in the left sequence
- for each for loop iteration append to result the two aligned sequen (seqAligned1,seqAligned2) and the score
'''

from input_data import blosum50

def compileMatrices(seq1,seq2,matrix,gp=-4):

    '''This function calculates the local dynamic programming (scoring) matrix F based
    on  the Smith & Waterman algorithm and the associate traceback matrix P'''

# initialization
    lRow=len(seq1)+1
    lCol=len(seq2)+1
    F=[[0]*lRow for x in range(lCol)]   # generates the empty scoring matrix F
    P=[['0']*lRow for x in range(lCol)] # generates the empty traceback matrix P
    maxScore=0
    index=[]
#iteration
    for r in range(1,lCol):
        for c in range(1,lRow):
            d=float(F[r-1][c-1])+matrix[seq1[c-1]+seq2[r-1]]             # d = diagonal
            u=float(F[r-1][c])+gp                                        # l = left
            l=float(F[r][c-1])+gp                                        # u = up
            F[r][c],P[r][c]=max(list(zip((0,d,u,l),('0','d','u','l'))))  # #Fill the F matrix with the maximal value
            if F[r][c]==maxScore:                                        # and the P matrix  with its direction
                maxScore=F[r][c]
                index+=[r,c]                # store the position of the highest score(s) in the matrix
            elif F[r][c]>maxScore:
                maxScore=F[r][c]
                index=[r,c]
    return (F,P,index)

def calcAlignment(F,P,index,seq1,seq2):

    '''The function finds the best local alignment(s) from the traceback matrix P.

    Backtracking process begins from the cell with the highest value in F and
    then traceback to first 0.
    '''

    result=[]                       # initialize an empty list where to store
    score=F[index[0]][index[1]]     # the aligned sequences and the score
    for i in range(0,len(index)-1,2):
        r=index[i]                  # highest value coordinates
        c=index[i+1]
        seqAligned1= ''
        seqAligned2 = ''
# Backtracking
        while F[r][c] != 0:          #  it terminates when it reaches the first 0 in F
            if P[r][c]== 'd':        # 'd' - the letters from two sequences are aligned
                seqAligned1= seq1[c-1]+seqAligned1
                seqAligned2=seq2[r-1]+seqAligned2
                r-=1
                c-=1
            elif P[r][c]== 'u':      # 'u' - a gap is introduced in seq1
                seqAligned1= '-'+seqAligned1
                seqAligned2= seq2[r-1]+seqAligned2
                r-=1
            elif P[r][c]== 'l':     # 'l' - a gap is introduced in seq2
                seqAligned2= '-'+seqAligned2
                seqAligned1= seq1[c-1]+seqAligned1
                c-=1
        result+=[seqAligned1]       # append to a list aligned sequences and their score
        result+=[seqAligned2]
        result+=[str(score)]
    return(result)

seq1="PALQPTQGAM"
seq2="QQMEELGMAP"

#from input_data import seq1,seq2
F,P,index=compileMatrices(seq1,seq2,blosum50)
result=calcAlignment(F,P,index,seq1,seq2)
print("\n".join(result)) #join list items using newline as separator
