import math

def couple(bases):
    couples=[]
    for x in range(len(bases)):
        for i in range (x,len(bases)):
            c=bases[x]+bases[i]
            couples.append(c)
    print (couples)
    return couples

def freq(count,leng):
    r=float(count)/leng
    return r

def calc_score(seq1,seq2,bases):
    frequences={}
    for i in bases:
        count=seq1.count(i)+seq2.count(i)
        frequence=freq(count,(len(s1)*2))
        frequences[i]=frequence
    couples=couple(bases)
    counter=1
    for coup in couples:
        for l in range(len(s1)):
            if coup == (seq1[l]+seq2[l]):
                counter+=1
            result=freq(counter,len(seq1))
            frequences[coup]=result
    return frequences

def calc_matrix(bases):
    frequences=calc_score(s1,s2,bases)
    #print(frequences)
    for i in range(len(bases)):
        matrix=[]
        for x in range (i,len(bases)):
            d1=frequences[bases[i]]
            d2=frequences[bases[x]]
            z=bases[i]+bases[x]
            n=frequences[str(z)]
            log_odd=math.log((float(n)/d1*d2),10)
            #print(str(z),log_odd)
            matrix.append(log_odd)
s1="ATCG"
s2="AGGT"
bases="ATCG"
#calc_score(s1,s2)
calc_matrix(bases)
