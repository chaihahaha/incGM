import networkx as nx
from networkx.algorithms import isomorphism
from itertools import permutations, combinations
from functools import reduce
import operator

def neighbor(a,G):
    # neighbor of a in G
    return reduce(operator.add,[[(i,j) for j in G.neighbors(i) if j not in a.nodes] for i in a.nodes],[])

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

def same(g1,g2):
    return (g1.nodes == g2.nodes) and (g1.edges == g2.edges)

def union(g1,g2,G):
    u_edges = set(g1.edges) + set(g2.edges)
    return G.edge_subgraph(u_edges)

class FRINGE:
    def __init__(self):
        self.MFS = []
        self.MIFS = []
        
class MNI_table:
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
    def frequent(self, tau):
        return self.support()>=tau

class FELS:
    def __init__(self,S):
        self.S = S
        self.mni = MNI_table(self.S)
        self.inverted_index = dict()
        self.invalid = set()
    def add(self, embedding):
        exists = False
        for i in embedding.nodes:
            if i in self.inverted_index.keys():
                if not exist(embedding, self.inverted_index[i]):
                    self.inverted_index[i].append(embedding)
                else:
                    exists = True
            else:
                self.inverted_index[i] = [embedding]
        if not exists:
            self.mni.add(dict(zip(list(self.S.nodes), list(embedding.nodes))))
            
            
class FELS_dict:
    def __init__(self):
        self.elem = dict()
    def add(self,S, fels):
        self.elem[S] = fels

def SEARCHLIMITED(S,newgraph,G):
    n_v = len(S)
    n_e = S.size()
    neighbor_graphs = []
    seed = [newgraph]
    for i in range(n_e-1):
        seed = neighbor_graph(seed,G)
    embeddings = []
    for i in seed:
        if len(i)==n_v:
            gm = isomorphism.GraphMatcher(i, S)
            is_iso = gm.is_isomorphic()
            if is_iso:
                embeddings.append(i)
    return embeddings

def has_nodes(g, s):
    for i in s:
        if g.has_node(i):
            return True
    return False

def FELSUpdate(embeds, S,tau):
    if S in fels_dict.elem.keys():
        valid_nodes = fels_dict.elem[S].inverted_index.keys()
        embeds.sort(key=lambda i: has_nodes(i,valid_nodes))
    for embedding in embeds:
        if S not in fels_dict.elem.keys():
            fels_dict.add(S, FELS(S))
        fels_dict.elem[S].add(embedding)
        if MNI(S,tau):
            break
        
        
def EVALUATE(G, tau, S):
    n_v = len(S)
    n_e = S.size()
    S_nodes = list(S.nodes)
    sgs_nodes = []
    components = nx.connected_components(G)
    for i in components:
        sgs_nodes += list(combinations(list(i),n_v))
    sgs = []
    for i in sgs_nodes:
        n_v_subsets = list(combinations(G.subgraph(i).edges,n_e))
        for j in n_v_subsets:
            n_e_subgraph = G.edge_subgraph(j)
            if nx.is_connected(n_e_subgraph):
                sgs.append(n_e_subgraph)
    sgs_nodes = [i.nodes for i in sgs]
    if S in fels_dict.elem.keys():
        valid_nodes = fels_dict.elem[S].inverted_index.keys()
        sgs.sort(key=lambda i: has_nodes(i,valid_nodes))
    mni = MNI(S)
    count_iso = 0
    for i in range(len(sgs)):
        gm = isomorphism.GraphMatcher(S, sgs[i])
        is_iso = gm.is_isomorphic()
        if is_iso:
            embedding = dict(zip(S_nodes,list(sgs_nodes[i])))
            mni.add(embedding)
            fels_dict.elem[S].add(sgs[i])
            fels_dict.elem[S].mni = mni
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
        for i in range(len(fringe.MFS)):
            MFSi = fringe.MFS[i]
            if len(MFSi) == len(S) and not same(MFSi,S):
                u = union(MFSi,S,G)
                if not EVALUATE(G,tau,u) and not exist(u,fringe.MIFS):
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

def MNI(S, tau):
    return fels_dict.elem[S].mni.frequent(tau)

fels_dict = FELS_dict()

def incGM_plus(G, fringe, tau, newedge):
    newgraph = nx.Graph()
    newgraph.add_edge(*newedge)
    if not subset(newgraph, G):
        fringe.MIFS.append(newgraph)
        
    G.add_edge(*e)
    for i in fels_dict.elem.keys():
        fels_dict.elem[i].invalid |= set(newedge)
    i = 0
    while 0 <= i <len(fringe.MIFS):
        S = fringe.MIFS[i]
        embeds = SEARCHLIMITED(S, newgraph,G)
        if not embeds:
            isFreq = EVALUATE(G,tau,S)
        else:
            FELSUpdate(embeds, S, tau)
            isFreq = MNI(S, tau)
        delete = UPDATEFRINGE(fringe, S, isFreq, tau, G)
        i = i + 1 - delete
    return fringe.MFS

G = nx.Graph()
fringe = FRINGE()
tau=5
with open("citeseerInt.cites","r") as f:
    txt = f.read()
edges = [[int(i) for i in j.split(",")] for j in txt.split("\n") if j]

for e in edges:
    print()
    print("MFS",[list(i.edges) for i in incGM(G,fringe, tau, e)])
