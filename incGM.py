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
    k = len(pattern)
    sgs_nodes = []
    components = nx.connected_components(G)
    for i in components:
        sgs_nodes.append(permutations(list(i),k))
    sgs_nodes = permutations(G.nodes, k)
    sgs =  [G.subgraph(i) for i in sgs_nodes]
    mni = MNI(pattern)
    for i in sgs:
        embedding = isomorphism.GraphMatcher(pattern, i)
        if embedding.is_isomorphic():
            mni.add(embedding.mapping)
    return mni.support()>=tau
def exist(u, arr):
    for j in arr:
        if same(j,u):
            return True
    return False
def PRUNE(fringe, S):
    fringe.MFS = [i for i in fringe.MFS if not subset(i,S)]
    return
def UPDATEFRINGE(fringe, S, isFreq, tau, G):
    # return # of deleted in MIFS
    count = 0
    if isFreq:
        if not exist(S,fringe.MFS):
            fringe.MFS.append(S)

        for i in fringe.MIFS:
            if same(i,S):
                fringe.MIFS.remove(i)
                count += 1
                break
        for i in fringe.MFS:
            if len(i) == len(S) and not same(i,S):
                u = union(i,S)
                if not exist(u,fringe.MIFS):
                    fringe.MIFS.append(u)
    return count
def subset(g1,g2):
    # return if g1 <= g2
    for i in g1.nodes:
        if not g2.has_node(i):
            return False
    return True

def incGM(G, fringe, tau, newedge):
    newgraph = nx.Graph()
    newgraph.add_edge(*newedge)
    if not subset(newgraph, G):
        fringe.MIFS.append(newgraph)
    
    G.add_edge(*e)
    i = 0
    while 0 <= i <len(fringe.MIFS):
        S = fringe.MIFS[i]
        isFreq = EVALUATE(G,tau,S)
        delete = UPDATEFRINGE(fringe, S, isFreq, tau, G)
        i = i + 1 - delete
    return fringe.MFS

G = nx.Graph()
fringe = FRINGE()
tau=1
edges = [(0,1),(0,2),(0,3),(0,4)]

for e in edges:
    print()
    print("MFS",[list(i.edges) for i in incGM(G,fringe, tau, e)])
