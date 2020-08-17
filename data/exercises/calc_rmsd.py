GYC=open("./1GYC.pdb","r")
MOD=open("./MODELLO.pdb","r")
def CA (pdb,x):
    ca=[]
    for line in pdb:
        line=line.rstrip()
        coord=[]
        if "CA" in line and "ATOM" in line:
            line=line.split()
            print(line)
            for i in range (x,x+3):
                coord.append(float(line[i]))
            ca.append(coord)
    return ca

l1=CA(GYC,6)
l2=CA(MOD,5)
GYC.close()
MOD.close()

def calc_rmsd(list1,list2,n_coord):
    s=0
    for i in range(len(list1)):
        for j in range(n_coord):
            x=(list1[i][j]-list2[i][j])**2
            s+=x
    rmsd=(s/len(list1))**(1/2)
    return rmsd

print(calc_rmsd(l1,l2,3))
