import pandas as pd
import numpy as np
import matplotlib as mpl
import json

f = open('sds.json','r')
signed_diffsets = json.load(f)
f.close()

print(f'read {len(signed_diffsets.keys())} data items\n')

# signed_diffsets is a python dictionary; each entry is the "name" of a set of parameters, e.g. "SDS(89,12,1,[89])"
# the dictionary entries contain 
#          "status": either "All", "Yes", "Open" or "No", (all known, exist, open, or known not to exist),
#          "comment": information about how the status is known,
#          "sets": a list of [P,N] sets for each known signed difference set
#          "G_rep": the representation of the group elements in P and N, if not given by the invariant factors

# these functions help to access the dictionary

#status of parameters
def get_status(D):
    if 'status' in signed_diffsets[D]:
        return signed_diffsets[D]['status']
    return None

def get_comment(D):
    if 'comment' in signed_diffsets[D]:
        return signed_diffsets[D]['comment']
    return None


# number of sets known for these parameters
def num_sets(D):
    if "sets" not in signed_diffsets[D]:
        return 0

    return len(signed_diffsets[D]["sets"])

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
def get_P(D,i):
    if 'sets' not in signed_diffsets[D] or (len(signed_diffsets[D]['sets']) <= i):
        print('error: no such set')
        return
    return signed_diffsets[D]['sets'][i][0]

def get_N(D,i):
    if 'sets' not in signed_diffsets[D] or (len(signed_diffsets[D]['sets']) <= i):
        print('error: no such set')
        return
    return signed_diffsets[D]['sets'][i][1]


# print one set
def print_set(D,i):
    S = D.split(',')
    v = int(S[0].split('(')[1])
    k = int(S[1])
    lam = int(S[2])
    G = S[3].split(')')[0]
    print(f'\n{D}',end='')
    if "G_rep" in signed_diffsets[D]:
        print(f' as elements of {signed_diffsets[D]["G_rep"]}',end='')
    print(f'\n\tP={signed_diffsets[D]["sets"][i][0]}, N={signed_diffsets[D]["sets"][i][1]}')
    

# get the ith (v,k,lambda) signed difference set in G, as a list
def get_sds(v,k,lam,G,i):
    sdsname = f'SDS({v},{k},{lam},{G})'.replace(' ','')
    if sdsname not in signed_diffsets:
        print(f'{sdsname} not in database')
        return

    D = signed_diffsets[sdsname]
    if ('sets' not in D) or (len(D['sets']) <= i):
        print('no such set')
        return

    if 'G_rep' in D:
        G = D['G_rep']

    return [v,k,lam,G,D['sets'][i][0],D['sets'][i][1]]

# convert an int or vector to an Additive Abelian Group element
def toGroup(G,g):
    # for a cyclic group elements are integers
    if G.is_cyclic():
        return G(vector([g]))

    # for noncyclic they're lists
    return G(g)

# convert an SDS to its representation as an element of the group ring of G
def sds_as_gp_ring_elt(D):

    G = AdditiveAbelianGroup(D[3])
    R = GroupAlgebra(G,ZZ)

    A = R.zero()
    P = D[4]
    N = D[5]

    A = 0
    for g in P:
        A = A +  R(toGroup(G,g))

    for g in N:
        A = A -  R(toGroup(G,g))

    return A


# get the ith (v,k,lambda) signed difference set in G, as a group ring element
def get_sds_in_groupring(v,k,lam,G,i):
    sdsname = f'SDS({v},{k},{lam},{G})'.replace(' ','')
    if sdsname not in signed_diffsets:
        print(f'{sdsname} not in database')
        return

    D = signed_diffsets[sdsname]
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
        print(A)
    for d in N:
        A = A - toGroupRing(R,d)

    return A

# get number of SDS(v,k,lam) sets, if any
# NOTE: this is different than num_sets() above in that its argument is the dictionary entry, not the SDS "name"
def set_count(sds):
    if 'sets' not in sds:
        return 0
    return len(sds['sets'])

# print out information about a given SDS
def get_sds_data(v,k,lam,G):
    sdsname = f'SDS({v},{k},{lam},{G})'.replace(' ','')
    if sdsname not in signed_diffsets:
        print(f'{sdsname} not in database')
        return

    D = signed_diffsets[sdsname]
    if D['status']=="All":
        if set_count(D)>1:
            print(f'There are exactly {set_count(D)} SDS({v},{k},{lam}) in group {G}')
        else:
            print(f'There is exactly {set_count(D)} SDS({v},{k},{lam}) in group {G}')

    if D['status']=="Yes":
        if set_count(D)>1:
            print(f'There are at least {set_count(D)} SDS({v},{k},{lam}) in group {G}')
        else:
            if set_count(D)>0:
                print(f'There is at least {set_count(D)} SDS({v},{k},{lam}) in group {G}')
            else:
                print(f'There is at least one SDS({v},{k},{lam}) in group {G}, but it is not in this dataset')

    if D['status']=="No":
            print(f'No SDS({v},{k},{lam}) exists in group {G}')

    if 'comment' in D:
        print(f'Reference: {D["comment"]}\n')

    if 'G_rep' in D:
        print(f'SDS given as elements of {D["G_rep"]}')

    for i in range(set_count(D)):
        sds = D['sets'][i]
        if set_count(D) > 1:
            print(f'{i}:\t',end='')
        print(f'P = {sds[0]}, N = {sds[1]}')

def get_cyclic_sds_data(v,k,lam):
    get_sds_data(v,k,lam,[v])

# find all (v,k,lambda) difference sets for any group
def get_sds_allgroups(v,k,lam):
    sdsname = f'SDS({v},{k},{lam}'
    sdslen = len(sdsname)
    for d in signed_diffsets:
        if d[:sdslen]==sdsname:
            G = d[sdslen+1:-1]
            get_sds_data(v,k,lam,G)
            print('')

# multiply group elts in group ring element by t
def gp_ring_elt_map(A,t):
    R = A.parent()
    G = R.group()
    
    B = R.zero()
    for a in A:
        B = B + a[1]*R(t*a[0])
        
    return B
    


# test by getting the group ring version
def is_sds(D):

    v = D[0]
    k = D[1]
    lam = D[2]
    G = AdditiveAbelianGroup(D[3])
    R = GroupAlgebra(G,ZZ)

    A = sds_as_gp_ring_elt(D)

    B = R.zero()
    for a in A:
        B = B + a[1]*R(-a[0])
        
    C = A*B

    for g in G:

        if abs(A.coefficient(g))>1:
            return False

        if g==G.zero():
            if C.coefficient(g) != k:
                return False
        else:
            if C.coefficient(g) != lam:
                return False
            
    return True


# code to create tables for showing a list of signed difference sets
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
    T['status'] += [signed_diffsets[D]['status']]
    T['comment'] += [signed_diffsets[D]['comment']]

def show_tab(T):
    df = pd.DataFrame(T)
    df = df.style.hide(axis='index')
    return df


