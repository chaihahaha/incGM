import networkx as nx
from networkx.algorithms import isomorphism
import time

def neighbor(a,G):
    # neighbor of a in G
    return G.edges(a.nodes)

def neighbor_graph(arr,G):
    sgss = []
    for a in arr:
        edges = neighbor(a,G)
        union_edges = set(a.edges)
        sgs = []
        for e in edges:
            union_edges.add(e)
            sgs.append(G.edge_subgraph(union_edges))
            union_edges.remove(e)
        sgss += sgs
    return sgss

def neighbor_graphs(a, G, n):
    seed = [a]
    for i in range(n):
        seed = neighbor_graph(seed,G)
    return [G.subgraph(i.nodes) for i in seed]

def same(g1,g2):
    return (g1.nodes == g2.nodes) and (g1.edges == g2.edges)

def union(g1,g2,G):
    u_nodes = set(g1.nodes) | set(g2.nodes)
    return G.subgraph(u_nodes)

class FRINGE:
    def __init__(self):
        self.MFS = []
        self.MIFS = []
        
class MNI_table:
    def __init__(self,S):
        self.nodes = S.nodes
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
    def frequent(self, tau):
        return self.support()>=tau

class FELS:
    def __init__(self,S):
        self.S = S
        self.mni = MNI_table(self.S)
        self.inverted_index = dict()
        self.embeddings = []
    def add(self, embedding):
        self.mni.add({v: k for k, v in embedding.items()})
            
            
class FELS_dict:
    def __init__(self):
        self.elements = dict()
    def add(self,S,embedding):
        exists = False
        for j in fels_dict.keys():
            if same(j,S):
                exists = True
        if not exists:
            self.elements[S] = FELS(S)
        self.elem(S).add(embedding)
    def keys(self):
        return self.elements.keys()
    def elem(self, S):
        for i in self.elements.keys():
            if same(i,S):
                return self.elements[i]
        return self.elements[S]

    def is_frequent(self, S, tau):
        return self.elem(S).mni.frequent(tau)


def SEARCHLIMITED(S,newgraph,G):
    n_v = len(S)
    n_e = S.size()
    search_region = neighbor_graphs(newgraph,G, n_v - 2)
    embeddings = []
    for graph in search_region:
        gm = isomorphism.GraphMatcher(graph, S)
        for i in gm.subgraph_isomorphisms_iter():
            embeddings.append(i)
    return embeddings


def FELSUpdate(embeds, S,tau):
    for embedding in embeds:
        fels_dict.add(S, embedding)
        if fels_dict.is_frequent(S, tau):
            break
        
        
def EVALUATE(G, tau, S):
    gm = isomorphism.GraphMatcher(G,S)
    for i in gm.subgraph_isomorphisms_iter():
        fels_dict.add(S, i)
        if fels_dict.is_frequent(S, tau):
            return True
    return False

def exist(u, arr):
    for j in arr:
        if same(j,u):
            return True
    return False

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
        for i in range(len(fringe.MFS)):
            MFSi = fringe.MFS[i]
            if len(MFSi) == len(S) and not same(MFSi,S):
                u = union(MFSi,S,G)
                if not EVALUATE(G,tau,u):
                    fringe.MIFS.append(u)
    return count

fels_dict = FELS_dict()

def incGM_plus(G, fringe, tau, newedge):
    newgraph = nx.Graph()
    newgraph.add_edge(*newedge)
    if not G.has_edge(*newedge):
        fringe.MIFS.append(newgraph)
        
    G.add_edge(*newedge)
    i = 0
    while 0 <= i <len(fringe.MIFS):
        S = fringe.MIFS[i]
        embeds = SEARCHLIMITED(S, newgraph,G)
        if not embeds:
            isFreq = EVALUATE(G,tau,S)
        else:
            FELSUpdate(embeds, S, tau)
            isFreq = fels_dict.is_frequent(S, tau)
        delete = UPDATEFRINGE(fringe, S, isFreq, tau, G)
        i = i + 1 - delete
    return fringe.MFS

G = nx.Graph()
fringe = FRINGE()
tau=200
with open("citeseerInt.cites","r") as f:
    txt = f.read()
edges = [[int(i) for i in j.split(",")] for j in txt.split("\n") if j]

for e in edges:
    print()
    print("MFS",[list(i.edges) for i in incGM_plus(G,fringe, tau, e)])
