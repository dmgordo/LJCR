import pandas as pd
import numpy as np
import matplotlib as mpl
import json

f = open('ds.json','r')
diffsets = json.load(f)
f.close()

print(f'read {len(diffsets.keys())} data items\n')

# diffsets is a python dictionary; each entry is the "name" of a set of parameters, e.g. "DS(11,5,2,[11])"
# the dictionary entries contain 
#          "status": either "All", "Yes", "Open" or "No", (all known, exist, open, or known not to exist),
#          "comment": information about how the status is known,
#          "sets": a list of each known difference set
#          "G_rep": the representation of the group elements in the sets, if not given by the invariant factors

# these functions help to access the dictionary

#status of parameters
def get_status(D):
    if 'status' in diffsets[D]:
        return diffsets[D]['status']
    return None

def get_comment(D):
    if 'comment' in diffsets[D]:
        return diffsets[D]['comment']

# number of sets known for these parameters
def num_sets(D):
    if "sets" not in diffsets[D]:
        return 0

    return len(diffsets[D]["sets"])

# pull parameters out from name D
def get_v(D):
    S = D.split(',')
    return int(S[0].split('(')[1])

def get_k(D):
    S = D.split(',')
    return int(S[1])

def get_lam(D):
    S = D.split(',')
    return int(S[2])

def get_G(D):
    S = D.split('[')[1].split(']')[0].split(',')
    G = []
    for i in range(len(S)):
        component = S[i]
        G += [int(component)]
    return G

#get the ith set as a list
def get_set(D,i):
    if 'sets' not in diffsets[D] or (len(diffsets[D]['sets']) <= i):
        print('error: no such set')
        return
    return diffsets[D]['sets'][i]




def get_ds(v,k,lam,G,i):
    dsname = f'DS({v},{k},{lam},{G})'.replace(' ','')
    if dsname not in diffsets:
        print(f'{dsname} not in database')
        return

    D = diffsets[dsname]
    if 'sets' not in D:
        print(f'no {dsname} difference sets in database')
        return

    if len(D['sets'])<=i:
        print(f'only {len(D["sets"])} {dsname} difference sets in database')
        return

    return [v,k,lam,G,D['sets'][i]]

# get number of DS(v,k,lam) sets, if any
def num_sets(ds):
    if 'sets' not in ds:
        return 0
    return len(ds['sets'])

# print out information about a given DS
def get_ds_data(v,k,lam,G):
    dsname = f'DS({v},{k},{lam},{G})'.replace(' ','')
    if dsname not in diffsets:
        print(f'{dsname} not in database')
        return

    D = diffsets[dsname]
    if D['status']=="All":
        if num_sets(D)>1:
            print(f'There are exactly {num_sets(D)} DS({v},{k},{lam}) in group {G}')
        else:
            print(f'There is exactly {num_sets(D)} DS({v},{k},{lam}) in group {G}')

    if D['status']=="Yes":
        if num_sets(D)>1:
            print(f'There are at least {num_sets(D)} DS({v},{k},{lam}) in group {G}')
        else:
            if num_sets(D)>0:
                print(f'There is at least {num_sets(D)} DS({v},{k},{lam}) in group {G}')
            else:
                print(f'There is at least one DS({v},{k},{lam}) in group {G}, but it is not in this dataset')

    if D['status']=="No":
            print(f'No DS({v},{k},{lam}) exists in group {G}')

    if 'comment' in D:
        print(f'Reference: {D["comment"]}\n')

    if 'G_rep' in D:
        print(f'DS given as elements of {D["G_rep"]}')

    for i in range(num_sets(D)):
        ds = D['sets'][i]
        if num_sets(D) > 1:
            print(f'{i}:\t',end='')
        for j in range(len(ds)):
            print(f'{ds[j]} ',end='')
        print('')


# find all (v,k,lambda) difference sets for any group
def get_ds_allgroups(v,k,lam):
    dsname = f'DS({v},{k},{lam}'
    dslen = len(dsname)
    for d in diffsets:
        if d[:dslen]==dsname:
            G = d[dslen+1:-1]
            get_ds_data(v,k,lam,G)
            print('')


# convert an int or vector to an Additive Abelian Group element
def toGroup(G,g):
    # for a cyclic group elements are integers
    if G.is_cyclic():
        return G(vector([g]))

    # for noncyclic they're lists
    return G(g)

# convert an DS to its representation as an element of the group ring of G
def ds_as_gp_ring_elt(D):

    print(D)
    G = AdditiveAbelianGroup(D[3])
    R = GroupAlgebra(G,ZZ)

    A = R.zero()
    S = D[4]

    A = 0
    for g in S:
        A = A +  R(toGroup(G,g))

    return A


# get the ith (v,k,lambda) signed difference set in G, as a group ring element
def get_ds_in_groupring(v,k,lam,G,i):
    dsname = f'DS({v},{k},{lam},{G})'.replace(' ','')
    if dsname not in diffsets:
        print(f'{dsname} not in database')
        return

    D = diffsets[dsname]
    if ('sets' not in D) or (len(D['sets']) <= i):
        print('no such set')
        return
    if 'G_rep' in D:
        R = GroupAlgebra(AdditiveAbelianGroup(D['G_rep']),ZZ)
    else:
        R = GroupAlgebra(AdditiveAbelianGroup(G),ZZ)

    A = R.zero()
    P = D['sets'][i][0]
    N = D['sets'][i][1]

    for d in P:
        A = A + toGroupRing(R,d)
    for d in N:
        A = A - toGroupRing(R,d)

    return A

# multiply group elts in group ring element by t
def gp_ring_elt_map(A,t):
    R = A.parent()
    G = R.group()
    
    B = R.zero()
    for a in A:
        B = B + a[1]*R(t*a[0])
        
    return B
    

# test by getting the group ring version
def is_ds(D):

    v = D[0]
    k = D[1]
    lam = D[2]
    G = AdditiveAbelianGroup(D[3])
    R = GroupAlgebra(G,ZZ)

    A = ds_as_gp_ring_elt(D)

    B = R.zero()
    for a in A:
        B = B + a[1]*R(-a[0])
        
    C = A*B

    for g in G:

        if A.coefficient(g)!=0 and  A.coefficient(g)!=1:
            return False

        if g==G.zero():
            if C.coefficient(g) != k:
                return False
        else:
            if C.coefficient(g) != lam:
                return False
            
    return True

# code to create tables for showing a list of difference sets
def init_tab():
    T = {}
    T['v'] = []
    T['k'] = []
    T['lambda'] = []
    T['n'] = []
    T['status'] = []
    T['comment'] = []
    return T

def add_tab_entry(T,D):
    v = int(D.split(',')[0].split('(')[1])
    k = int(D.split(',')[1])
    lam = int(D.split(',')[2].split(')')[0])
    n = k-lam
    T['v'] += [v]
    T['k'] += [k]
    T['lambda'] += [lam]
    T['n'] += [n]
    T['status'] += [diffsets[D]['status']]
    T['comment'] += [diffsets[D]['comment']]

def show_tab(T):
    df = pd.DataFrame(T)
    df = df.style.hide(axis='index')
    return df


