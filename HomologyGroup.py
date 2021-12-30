#!/usr/bin/env python
# coding: utf-8

# In[348]:


from math import gcd
import numpy as np
import numpy.linalg as LA
from sympy.combinatorics import Permutation
import smith_normalize as smith


def perm_sgn(l1,l2):
    return Permutation([ l2.index(i) for i in l1 ]).signature()

def simplex_boudary(s,C1):
    chain = np.zeros(len(C1))
    for i in range(len(s)):
        sb = s.copy()
        sb.pop(i)
        isb = [ set(s) for s in C1 ].index(set(sb))
        chain[isb] = pow(-1,i)*perm_sgn(C1[isb],sb)
    return chain
    
def calc_ith_boundary(C1,C2):
    D = np.zeros( (len(C2),len(C1)) ).astype(int)
    for j,s in enumerate(C2):
        D[j] = simplex_boudary(s,C1)
    return D.T

def calc_boundary_operators(SimplicicalComplex):
    size = list(set([ len(s) for s in SimplicicalComplex ]))
    ChainGroupList = [ [ s for s in SimplicicalComplex if len(s) == size[i] ] for i in range(len(size)) ]
    BoundaryList = [np.zeros( (1,len(ChainGroupList[0])) ).astype(int)]
    for i in range(len(size)-1):
        BoundaryList.append( calc_ith_boundary(ChainGroupList[i],ChainGroupList[i+1]) )
    BoundaryList.append(np.array([]))
    return BoundaryList


def colR(M):
    if(M.shape[0]*M.shape[1]==0): return M
    for j in range(M.shape[1]):
        i = list(map(lambda a:np.where(a!=0)[0][0],M.T))[j]
        for j2 in range(M.shape[1]):
            if(j2==j): continue
            c = np.sign(M[i,j2])*(abs(M[i,j2])//M[i,j])
            M[:,j2] -= c*M[:,j]
    return M

def colR_same_tor(M,tortion):
    tor = tortion.copy()
    if(M.shape[0]*M.shape[1]==0): return M
    for j in range(M.shape[1]):
        i = list(map(lambda a:np.where(a!=0)[0][0],M.T))[j]
        for j2 in range(M.shape[1]):
            gt = gcd(tor[j],tor[j2])
            # avoid mixing non-trivial tortion
            if(j2==j or (tor[j2]>tor[j]) or (tor[j]>1 and tor[j2]>1 and gt==1) ): continue
            c = np.sign(M[i,j2])*(abs(M[i,j2])//M[i,j])
            tor[j2] = gt
            M[:,j2] -= c*M[:,j]
    return M,tor

def calc_cohomology(D):
    r = LA.matrix_rank(D)
    invU,S,V = smith.Smith_Normalization(D)
    U = LA.inv(invU)
    Z = V[:,r:]
    tortion = np.diag(S[:r,:r])
    B = U[:,:r]
    Z = colR(Z)
    B, tortion = colR_same_tor(B,tortion)
    return B,tortion,Z

def calc_ith_homology(D1,D2):
    if( (D1==0).all() ): Z1 = np.identity(D1.shape[1]).astype(int)
    else: B0, tortion0, Z1 = calc_cohomology(D1)
    if( len(D2)==0 ):
        B1 = np.array([[]])
        tortion1 = np.array([])
    else: B1,tortion1, Z2 = calc_cohomology(D2)
    Z_nonzero = list(map(lambda z:np.where(z!=0)[0][0],Z1.T))
    B_nonzero = list(map(lambda z:np.where(z!=0)[0][0],B1.T))
    non_trivial = [ [i,1] for i,z in enumerate(Z_nonzero) if not(z in B_nonzero) ]
    non_trivial += [ [np.where(Z_nonzero == ib )[0][0],tortion1[i]] for i,ib in enumerate(B_nonzero) if tortion1[i]>1 ]
    return Z1[:,[ x[0] for x in non_trivial ]],[ x[1] for x in non_trivial ]

def calc_HomologyGroupList(SimplicicalComplex):
    BoundaryList = calc_boundary_operators(SimplicicalComplex)
    HomologyGroupList = []
    for i in range( len(BoundaryList)-1 ):
        HomologyGroupList.append( calc_ith_homology(BoundaryList[i],BoundaryList[i+1]) )
    return HomologyGroupList

