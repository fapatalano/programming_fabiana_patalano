import numpy as np

def doMatrixP(m1,m2):
    x=np.dot(m1,m2)
    print(x)
    return x

'''def molt_mat(m1,m2):
    if len(m1)==len(m2[0]):
        matrix=[]
        for i in range(len(m1)):    #scorre le righe
            row=[]
            for j in range(len(m2[0])):  #scorre le colonne
                sum=0
                for k in range(len(m2)):
                    sum+=m1[i][k]*m2[k][j]
                row.append(sum)
            matrix.append(row)
        return matrix
    else:
        return False'''

m1=[[1,2],[3,4],[5,6]]
m2=[[1,2,3],[3,4,5]]
doMatrixP(m1,m2)
#print(molt_mat(m1,m2))
