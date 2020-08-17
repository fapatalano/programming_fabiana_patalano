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
    d= the two letters are aligned
    u= introduce a gap in the top sequence
    l= introduce a gap in the left sequence
- for each for loop iteration append to result the two aligned sequences (seqAligned1,seqAligned2) and the score
'''

def compileMatrices(seq1,seq2,matrix,gp=-4):

    '''This function calculates the local dynamic programming (scoring) matrix F based
    on  the Smith & Waterman algorithm and the associate traceback matrix P'''

# initializatione
    nCol=len(seq1)+1                    # seq1 is the top sequence -> define the number of columns
    nRow=len(seq2)+1                    # seq2 is the left sequence -> define the number of rows
    F=[[0]*nCol for x in range(nRow)]   # populates scoring matrix F with 0
    P=[['0']*nCol for x in range(nRow)] # populates traceback matrix P with 0
    maxScore=0
    index=[]
#iteration
    for r in range(1,nRow):
        for c in range(1,nCol):
            d=F[r-1][c-1]+matrix[seq1[c-1]+seq2[r-1]]              # d = diagonal
            u=F[r-1][c]+gp                                         # l = left
            l=F[r][c-1]+gp                                         # u = up
            F[r][c],P[r][c]=max(zip((0,d,u,l),('0','d','u','l')))  # Fill the F matrix with the max value
            if F[r][c]==maxScore:                                  # and the P matrix  with its direction
                maxScore=F[r][c]
                index+=[[r,c]]                # index stores the position of the highest score(s) in the matrix
            elif F[r][c]>maxScore:
                maxScore=F[r][c]
                index=[[r,c]]
    return (F,P,index)

def calcAlignment(F,P,index,seq1,seq2):

    '''The function finds the best LOCAL alignment(s) from the traceback matrix P.

    Backtracking process begins from the cell with the highest value in F and
    then traceback to first 0.
    '''
    result=[]                              # result will store the aligned sequences and the score
    for i in range(len(index)):
        score=F[index[i][0]][index[i][1]]  # final score of the alignment
        r=index[i][0]                      # highest value coordinates
        c=index[i][1]
        seqAligned1= ''                    # initialization of two variables
        seqAligned2 = ''                   # that will be the  final aligned sequences
# Backtracking
        while F[r][c] != 0:                # it stops when it reaches the first 0 in F matrix
            if P[r][c]== 'd':              # 'd' - the letters from two sequences are aligned
                seqAligned1= seq1[c-1]+seqAligned1
                seqAligned2= seq2[r-1]+seqAligned2
                r-=1
                c-=1
            elif P[r][c]== 'u':            # 'u' - a gap is introduced in seq1
                seqAligned1= '-'+seqAligned1
                seqAligned2= seq2[r-1]+seqAligned2
                r-=1
            elif P[r][c]== 'l':            # 'l' - a gap is introduced in seq2
                seqAligned2= '-'+seqAligned2
                seqAligned1= seq1[c-1]+seqAligned1
                c-=1
        result.append(seqAligned1)         # append to result aligned sequences and their score
        result.append(seqAligned2)
        result.append(str(score))
    return result

from input_data import blosum50,seq1,seq2
F,P,index=compileMatrices(seq1,seq2,blosum50)
result=calcAlignment(F,P,index,seq1,seq2)
print("\n".join(result)) #join each element of result using newline as separator
