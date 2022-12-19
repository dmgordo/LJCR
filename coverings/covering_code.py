import pandas as pd
import numpy as np
import matplotlib as mpl

# coverdata is a python dictionary; each entry is the "name" of a set of parameters, e.g. "C(7,3,2)"
# the dictionary entries contain 
#          "size": the number of blocks in the best known such covering design
#          "low_bd": the best lower bound known
#          "imps": improvements submitted over the years
#                       each improvement is a list containing:
#                            size of the covering
#                            method used to produce it
#                            submitter
#                            a timestamp of when it was submitted

# pull parameters out from name C
def get_v(C):
    S = C.split(',')
    return int(S[0].split('(')[1])

def get_k(C):
    S = C.split(',')
    return int(S[1])

def get_t(C):
    S = C.split(',')
    return int(S[2].split(')')[0])


def print_improvement(imp):
    sz = imp[0]
    method = imp[1]
    creator = imp[2]
    timestamp = imp[3]
    if len(method)>0:
        print(f'{imp[0]}:  {imp[3]} by {imp[2]} using {imp[1]}')
    else:
        print(f'{imp[0]}:  {imp[3]} by {imp[2]}')


def get_cover_data(v,k,t):
    covname = f'C({v},{k},{t})'
    if covname not in coverdata:
        print(f'{covname} not in database')
        return

    CD = coverdata[covname]
    sz = CD['size']
    lb = CD['low_bd']
    imp = CD['imps'][0]
    if sz == lb:
        print(f'{covname} = {sz}')
    else:
        print(f'{lb} <= {covname} <= {sz}')
    print(f'last update ',end='')
    print_improvement(imp)

def get_cover(v,k,t):
    covname = f'C({v},{k},{t})'
    if covname not in covers:
        print(f'{covname} not in database')
        return

    C = covers[covname]
    print(f'{covname} has {len(C)} blocks')
    return [v,k,t,C]

def show_cover(v,k,t):
    covname = f'C({v},{k},{t})'
    if covname not in covers:
        print(f'{covname} not in database')
        return

    C = covers[covname]
    print(f'{covname} has {len(C)} blocks')
    for i in range(len(C)):
        print(f'{C[i]}')


# code to create tables for showing a list of covering designs
def init_tab():
    T = {}
    T['v'] = []
    T['k'] = []
    T['t'] = []
    T['size'] = []
    T['lower bd'] = []
    T['creator'] = []
    T['method'] = []
    T['timestamp'] = []
    return T

def add_tab_entry(T,C):
    v = C.split(',')[0].split('(')[1]
    k = C.split(',')[1]
    t = C.split(',')[2].split(')')[0]
    T['v'] += [v]
    T['k'] += [k]
    T['t'] += [t]
    T['size'] += [coverdata[C]['size']]
    T['lower bd'] += [coverdata[C]['low_bd']]
    T['creator'] += [coverdata[C]['imps'][0][2]]
    T['method'] += [coverdata[C]['imps'][0][1]]
    T['timestamp'] += [coverdata[C]['imps'][0][3]]

def show_tab(T):
    df = pd.DataFrame(T)
    df = df.style.hide(axis='index')
    return df

# function to show a table of improvements for a given (v,k,t)
def show_history(v,k,t):
    T = {}
    T['size'] = []
    T['creator'] = []
    T['method'] = []
    T['timestamp'] = []
    indx = []
    covname = f'C({v},{k},{t})'
    if covname not in coverdata:
        print(f'{covname} not in database')
        return

    CD = coverdata[covname]
    sz = CD['size']
    lb = CD['low_bd']
    imp = CD['imps']
    if sz == lb:
        print(f'{covname} = {sz}')
    else:
        print(f'{lb} <= {covname} <= {sz}')
    
    print('update history:')
    for i in range(len(imp)):
        T['size'] += [int(imp[i][0])]
        indx += [int(imp[i][0])]
        T['method'] += [imp[i][1]]
        T['creator'] += [imp[i][2]]
        T['timestamp'] += [imp[i][3]]

    table = pd.DataFrame(T)
    table = table.style.hide(axis='index')
    return table

# check that a list of blocks is a covering
# this will take a long time for large parameters
def is_cover(C):
    v = C[0]
    k = C[1]
    t = C[2]
    L = C[3][:-1]

    for c in Combinations(range(1,v+1),t):
        found_it = False
        for b in L:
            if Set(c).issubset(Set(b)):
                found_it = True
                break

        if found_it == False:
            return False

    return True

