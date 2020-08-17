def seq_align(seq1,seq2):
    if len(seq1) != len(seq2):
        return "Two sequences have different lenght"
    score=0
    for i in range(len(seq1)):
        if seq1[i]==seq2[i]:
            score+=1
        else:
             score-=1
    return score

seq1="AAAAGCTAT"
seq2="GTCTCCTC"
print(seq_align(seq1,seq2))
