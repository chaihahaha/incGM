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
        v = self.supp.values()
        return min(v) if v else 0
    
def neighbor(a,G):
    # neighbor of a in G
    return reduce(operator.add,[[(i,j) for j in G.neighbors(i) if j not in a.nodes] for i in a.nodes],[])

def same(g1,g2):
    return (g1.nodes == g2.nodes) and (g1.edges == g2.edges)
def union(g1,g2):
    g=g1.copy()
    g.add_nodes_from(g2.nodes)
    g.add_edges_from(g2.edges)
    return g

class FRINGE:
    def __init__(self):
        self.MFS = []
        self.MIFS = []
        
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

def UPDATEFRINGE(fringe, S, isFreq, tau, G):
    if isFreq:
        exist = False
        for i in fringe.MFS:
            if same(i,S):
                exist=True
        if not exist:
            print('1', S.edges)
            fringe.MFS.append(S)

        for i in fringe.MIFS:
            if same(i,S):
                print('2',i.edges)
                fringe.MIFS.remove(i)
                break
        for i in fringe.MFS:
            if len(i) == len(S) and not same(i,S):
                u = union(i,S)
                if not EVALUATE(G, tau,u):
                    exist=False
                    for j in fringe.MIFS:
                        if same(j,u):
                            exist=True
                    if not exist:
                        print('3',u.edges)
                        fringe.MIFS.append(u)
    return
