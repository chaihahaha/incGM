import networkx as nx
from networkx.algorithms import isomorphism
import time

def neighbor(a,G):
    # a is a node set
    # neighbor of a in G
    return frozenset(G.edge_subgraph(G.edges(a)))

def neighbor_graphn(a, G, n):
    # a is a node set, G is a graph, n is # of iteration
    for i in range(n):
        a = neighbor(a,G)
    return a

class FRINGE:
    def __init__(self):
        self.MFS = []
        self.MIFS = []
    def addMFS(self, mfs):
        if mfs not in self.MFS:
            self.MFS.append(mfs)
            return True
        return False
    def addMIFS(self, mifs):
        if mifs not in self.MIFS:
            self.MIFS.append(mifs)
            return True
        return False
    def removeMIFS(self, mifs):
        if mifs in self.MIFS:
            self.MIFS.remove(mifs)
            return True
        return False

class MNI_table:
    def __init__(self,S_nodes):
        self.nodes = S_nodes
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
    def __init__(self,S_nodes):
        self.S_nodes = S_nodes
        self.mni = MNI_table(self.S_nodes)
        self.inverted_index = dict()
        self.embeddings = []
    def add(self, embedding):
        self.mni.add(embedding)


class FELS_dict:
    def __init__(self):
        self.elements = dict()

    def add(self,G2S_embedding):
        # index only by node is enough
        S2G_embedding = {v: k for k, v in G2S_embedding.items()}
        S_nodes = frozenset(S2G_embedding.keys())
        subG = frozenset(G2S_embedding.keys())
        if S_nodes not in fels_dict.keys():
            self.elements[S_nodes] = FELS(S_nodes)
        if subG not in fels_dict.keys():
            self.elements[subG] = FELS(subG)
        self.elem(S_nodes).add(S2G_embedding)
        self.elem(subG).add(G2S_embedding)
        return

    def keys(self):
        return self.elements.keys()

    def elem(self, S_nodes):
        return self.elements[S_nodes]

    def is_frequent(self, S_nodes, tau):
        return self.elem(S_nodes).mni.frequent(tau)


def SEARCHLIMITED(S_nodes,newnodes,G):
    n_v = len(S_nodes)
    search_region = neighbor_graphn(newnodes,G, n_v - len(newnodes))
    embeddings = []
    gm = isomorphism.GraphMatcher(G.subgraph(search_region), G.subgraph(S_nodes))
    for i in gm.subgraph_isomorphisms_iter():
        embeddings.append(i)
    return embeddings


def FELSUpdate(embeds, S_nodes,tau):
    for embedding in embeds:
        fels_dict.add(embedding)
        if fels_dict.is_frequent(S_nodes, tau):
            break


def EVALUATE(G, tau, S_nodes):
    isFreq = False
    gm = isomorphism.GraphMatcher(G,G.subgraph(S_nodes))
    for i in gm.subgraph_isomorphisms_iter():
        fels_dict.add(i)
        if fels_dict.is_frequent(S_nodes, tau):
            isFreq = True
            break
    return isFreq

def UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G):
    deleted = False
    if isFreq:
        added = fringe.addMFS(S_nodes)
        deleted = fringe.removeMIFS(S_nodes)

        for i in range(len(fringe.MFS)):
            MFSi = fringe.MFS[i]
            if len(MFSi) == len(S_nodes) and MFSi != S_nodes:
                u = MFSi | S_nodes
                if u not in fringe.MIFS:
                    if not EVALUATE(G,tau,u):
                        joined = fringe.addMIFS(u)
    return deleted

fels_dict = FELS_dict()

def incGM_plus(G, fringe, tau, newgraph):
    G.add_edges_from(newgraph.edges)
    newnodes = frozenset(newgraph.nodes)
    fringe.MIFS.append(newnodes)
    i = 0
    while 0 <= i <len(fringe.MIFS):

        S_nodes = fringe.MIFS[i]
        embeds = SEARCHLIMITED(S_nodes, newnodes,G)
        if not embeds:
            isFreq = EVALUATE(G,tau,S_nodes)
            print("EVALUATE!")
        else:
            FELSUpdate(embeds, S_nodes, tau)
            isFreq = fels_dict.is_frequent(S_nodes, tau)
        delete = UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G)
        i = i + 1 - int(delete)
    return fringe.MFS

G = nx.Graph()
fringe = FRINGE()
tau=200
with open("citeseerInt.cites","r") as f:
    txt = f.read()
edges = [[int(i) for i in j.split(",")] for j in txt.split("\n") if j]

count = 0
batchsize = 10
for i in range(100):
    print()
    tik = time.time()
    newgraph = nx.Graph([edges[i]])
    print("MFS",[list(i.edges) for i in incGM_plus(G,fringe, tau, newgraph)])
    tok = time.time()
    count += tok-tik
    print(tok-tik, count/G.size())

for i in range(100,len(edges)//batchsize):
    print()
    tik = time.time()
    newgraph = nx.Graph(edges[i:i+batchsize])
    print("MFS",[list(i.edges) for i in incGM_plus(G,fringe, tau, newgraph)])
    tok = time.time()
    count += tok-tik
    print(tok-tik, count/G.size())
