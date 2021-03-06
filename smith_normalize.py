#!/usr/bin/env python
# coding: utf-8

# reference: https://qiita.com/yuji0001/items/64dc97cd4dcebf83d0a8


import numpy as np
import numpy.linalg as LA

def Eij(i,j,n):
    E = np.eye(n)
    E[i,i] = 0
    E[j,j] = 0
    E[i,j] = 1
    E[j,i] = 1
    return E

def Ei(i,n):
    E = np.eye(n)
    E[i,i] = -1
    return E

def Ec(i,j,c,n):
    E = np.eye(n)
    E[i,j] = c
    return E    

# A[k:,k:]上の0以外の最小の絶対値をA[k,k]に移動する変換
def moveMN(A,k):
    tmp_A = A[k:,k:]
    a = np.abs(tmp_A[tmp_A != 0]).min()
    i = np.where(np.abs(tmp_A) == a)[0][0] + k
    j = np.where(np.abs(tmp_A) == a)[1][0]+ k
    P = Eij(k,j,A.shape[1])
    invQ = Eij(i,k,A.shape[0])
    B = invQ.dot(A).dot(P)
    if B[k,k]<0:
        Pi =Ei(k,A.shape[1])
        B = B.dot(Pi)
        P = P.dot(Pi)
    return invQ.astype(int),B.astype(int),P.astype(int)

#A[k,k]を使って，A[k+1:,k]を0に整える変換．(整いきれなかったら要素はA[k,k]より小さくなる)
def rowR(A,k):
    B = A.copy()
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    for i in range(k+1,A.shape[0]):
        q = A[i,k]//A[k,k]
        #残渣
        r = A[i,k]%A[k,k]
        invQi = Ec(i,k,-q,A.shape[0])
        B = invQi.dot(B)
        invQ = invQi.dot(invQ)
    return invQ.astype(int),B.astype(int),P.astype(int)

#A[k,k]を使って，A[k,k+1]を0に整える変換．(整いきれなかったら要素はA[k,k]より小さくなる)
def colR(A,k):
    B = A.copy()
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    for i in range(k+1,A.shape[1]):
        q = A[k,i]//A[k,k]
        #残渣
        r = A[k,i]%A[k,k]
        Pi = Ec(k,i,-q,A.shape[1])
        B = B.dot(Pi)
        P = P.dot(Pi)
    return invQ.astype(int),B.astype(int),P.astype(int)

# A[k+1:,k+1:]においてA[k,k]で割り切れない要素A[i,j]をA[k,k]の残差に変換する変換
def remR(A,k):
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    #  Find i,j
    i = np.where(A[k+1:,k+1:]%A[k,k] !=0)[0][0] +k+1
    j = np.where(A[k+1:,k+1:]%A[k,k] !=0)[1][0] +k+1
    q = A[i,j]//A[k,k]
    r = A[i,j]%A[k,k]
    invQi = Ec(i,k,q,A.shape[0])
    Pi = Ec(k,j,-1,A.shape[1])
    B = invQi.dot(A).dot(Pi)
    P = P.dot(Pi)
    invQ = invQi.dot(invQ)
    return invQ.astype(int),B.astype(int),P.astype(int)


# Main Function
def Smith_Normalization(A):
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    A0 = A.copy()
    # limit of optimization
    N = 1000
    for k in range(min(A0.shape)):
        # If A0[k:,k:] is zero matrix, then stop calculation
        if np.sum(np.abs(A0[k:,k:]))==0:
            break
        for i in range(N):
            if i == N-1 : 
                print("Error: Time Out")
            # minimize A[k,k]
            invQi,A1,Pi = moveMN(A0,k)
            invQ = invQi.dot(invQ)
            P = P.dot(Pi)
            # make zero row A[k+1:,k]
            invQi,A2,Pi = rowR(A1,k)
            invQ = invQi.dot(invQ)
            P = P.dot(Pi)
            # if row A2[k+1:,k] is zero vector ?
            if np.abs(A2[k+1:,k]).sum() ==0:
                # make zero col A[k,k+1:]
                invQi,A3,Pi = colR(A2,k)
                invQ = invQi.dot(invQ)
                P = P.dot(Pi)
                # if col A3[k+1:,k] is zero vector ?
                if np.abs(A3[k,k+1:]).sum() ==0:
                    # A[k,k]|A[k+1:,k+1:]?
                    if np.sum(A3[k+1:,k+1:]%A3[k,k]) == 0:                    
                        A0 = A3.copy()            
                        break;
                    else:
                        # reduce A[k+1:,k+1:]
                        invQi,A0,Pi = remR(A3,k)
                        invQ = invQi.dot(invQ)
                        P = P.dot(Pi)
                else:
                    A0 = A3.copy()            
            else:
                A0 = A2.copy()

    B = A0.copy().astype(int)
    P = P.astype(int)
    invQ = invQ.astype(int)
    return invQ,B,P


# Check result
#invU,B,V = Smith_Normalization(A)
#A=U*B*inv(V)


