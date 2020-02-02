import networkx as nx
from networkx.algorithms import isomorphism
from itertools import permutations
from functools import reduce
import operator
class MNI:
    def __init__(self,pattern):
        self.nodes = pattern.nodes
        self.table = dict()
        for i in self.nodes:
            self.table[i]=dict()
        self.supp = dict()
    def add(self, dic):
        for i in self.nodes:
            if dic[i] in self.table[i].keys():
                self.table[i][dic[i]] += 1
            else:
                self.table[i][dic[i]] = 1
            self.supp[i] = len(self.table[i].values())
    def support(self):
        return min(self.supp.values())
    
def neighbor(a,G):
    # neighbor of a in G
    return reduce(operator.add,[[(i,j) for j in G.neighbors(i) if j not in a.nodes] for i in a.nodes],[])

class FRINGE:
    def __init__(self):
        self.MFS = set()
        self.MIFS = set()
        
def EVALUATE(G, tau, pattern):
    sgs_nodes = permutations(G.nodes, k)
    sgs =  [G.subgraph(i) for i in sgs_nodes]
    mni = MNI(pattern)
    for i in sgs:
        embedding = isomorphism.GraphMatcher(pattern, i)
        if embedding.is_isomorphic():
            mni.add(embedding.mapping)
    print(mni.table)
    return mni.support()>=tau
