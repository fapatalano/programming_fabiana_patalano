#!/bin/python3

class svc():

    class_code = {'H': '0', 'E': '1', 'C': '2','-':'2'}
    k=5

    def __init__(self,profile_path,cv_path,dssp_path,w):
        self.w=w
        self.profile_path = profile_path
        self.cv_path=cv_path
        self.dssp_path=dssp_path
        self.profile_id=svc.extract_id(profile_path)
        self.cv_set=svc.k_fold_ids(self)
        self.X,self.Y,self.test_fold=svc.encode(self)


    def extract_id(path):
        id_list=[]
        for id in sorted(os.listdir(path)):
            id_list.append(id.split('.p')[0])
        return id_list

    def k_fold_ids(self):
        k_fold={}
        k=0
        for cv_set in sorted(os.listdir(self.cv_path)):
            tmp=[]
            with open (self.cv_path+cv_set) as file:
                [tmp.append(i.strip()) for i in file.readlines() if i.strip() in self.profile_id]
                k_fold[k]=tmp
                k+=1
        return k_fold
#
    def encode(self):
        '''
        The encoding required by scikitlearn module is in the format:
        '''
        Y=[]
        test_fold=[]
        for id in  self.profile_id:
            [test_fold.append(key) for key, value in self.cv_set.items() if id in value]
            with open(self.profile_path+id+'.prf') as file:
                profile=pd.read_csv(file,sep='\s+',header=None).iloc[:,1:].values
                padding= np.zeros((self.w//2,20))
                profile=np.vstack((padding,profile,padding))
                dssp=open(self.dssp_path+id+'.dssp').readlines()[1].strip()
                X=profile[0:self.w].flatten()
                for i in range(1,len(dssp)):
                    Y+=svc.class_code[dssp[i]]
                    row=profile[i:i+self.w].flatten()
                    X=np.vstack((X,row))
        return X,Y,test_fold

    def fit_model(self):
        ps = PredefinedSplit(self.test_fold)
        PredefinedSplit(self.test_fold)
        for train_index, test_index in ps.split():
            print("TRAIN:", train_index, "TEST:", test_index)

        # X_train, X_test = X[train_index], X[test_index]
        # y_train, y_test = y[train_index], y[test_index]

if __name__ == '__main__':

    import argparse
    import os
    import pandas as pd
    import numpy as np
    from sklearn.model_selection import GridSearchCV,PredefinedSplit
    from sklearn.metrics import classification_report
    from sklearn.svm import SVC

    parser= argparse.ArgumentParser()
    parser.add_argument("--profile",help="path to the profile directory")
    parser.add_argument("--dssp",help="path to the dssp directory")
    parser.add_argument("--w",help="window size", type=int)
    parser.add_argument("--cv",help="path to the cv set directory")
    args=parser.parse_args()
    SVM=svc(args.profile,args.cv,args.dssp,args.w)
    SVM.fit_model()
