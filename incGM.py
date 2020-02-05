import networkx as nx
from networkx.algorithms import isomorphism
import matplotlib.pyplot as plt
import time

def neighbor(a,G):
    # a is a node set
    # neighbor of a in G
    return frozenset(G.edge_subgraph(G.edges(a)))

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
        self.embeddings = []
        self.inverted_index = dict() # dict of set of node sets(subgraph)
    def add(self, embedding):
        self.mni.add(embedding)
    def add_inverted(self, embedding):
        emb_nodes = frozenset(embedding.keys())
        for i in embedding.keys():
            if i not in self.inverted_index.keys():
                self.inverted_index[i] = set()
            self.inverted_index[i].add(emb_nodes)


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
        self.elem(S_nodes).add_inverted(G2S_embedding)
        self.elem(subG).add(G2S_embedding)
        self.elem(subG).add_inverted(S2G_embedding)
        return

    def keys(self):
        return self.elements.keys()

    def elem(self, S_nodes):
        return self.elements[S_nodes]

    def is_frequent(self, S_nodes, tau,G):
        if not nx.is_connected(G.subgraph(S_nodes)):
            return False
        return self.elem(S_nodes).mni.frequent(tau)

    def iso_graphs(self, S_nodes):
        gs = set()
        inv = self.elem(S_nodes).inverted_index
        for i in inv.keys():
            for j in inv[i]:
                gs.add(j)
        return gs


def SEARCHLIMITED(S_nodes,search_region,G):
    embeddings = False
    if not nx.is_connected(G.subgraph(S_nodes)):
        return False
    while len(search_region) < len(S_nodes):
        old_region = search_region
        search_region = neighbor(search_region,G)
        if old_region == search_region:
            break
    gm = isomorphism.GraphMatcher(G.subgraph(search_region), G.subgraph(S_nodes))
    for i in gm.subgraph_isomorphisms_iter():
        embeddings = True
        fels_dict.add(i)
        if fels_dict.is_frequent(S_nodes, tau,G):
            break
    return embeddings

def EVALUATE(G, tau, S_nodes):
    if not nx.is_connected(G.subgraph(S_nodes)):
        return False
    components = nx.connected_components(G)
    for cnodes in components:
        gm = isomorphism.GraphMatcher(G.subgraph(cnodes),G.subgraph(S_nodes))
        for i in gm.subgraph_isomorphisms_iter():
            fels_dict.add(i)
            if fels_dict.is_frequent(S_nodes, tau, G):
                return True
    return False

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
    fringe.addMIFS(newnodes)
    i = 0
    while 0 <= i <len(fringe.MIFS):
        S_nodes = fringe.MIFS[i]
        embeds = SEARCHLIMITED(S_nodes, newnodes,G)
        if not embeds:
            isFreq = EVALUATE(G,tau,S_nodes)
        else:
            isFreq = fels_dict.is_frequent(S_nodes, tau, G)
        delete = UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G)
        i = i + 1 - int(delete)
    return fringe.MFS

base = nx.gnm_random_graph(15,25,1)
pos = nx.spring_layout(base)
nx.draw(base,pos=pos)
plt.savefig("base.png")
plt.clf()
G = nx.Graph()
tau = 7
fringe = FRINGE()
for e in base.edges:
    tik = time.time()
    incGM_plus(G,fringe,tau,base.subgraph(e))
    tok = time.time()
    print(tok-tik)
distinct = [i for i in fringe.MFS]

for i in range(len(distinct)-1):
    j = i+1
    while 0<=j<len(distinct):
        gi,gj = (G.subgraph(distinct[i]), G.subgraph(distinct[j]))
        if nx.is_isomorphic(gi,gj):
            distinct.pop(j)
            j -= 1
        j += 1
for i in distinct:
    nx.draw(G.subgraph(i),pos=pos)
    plt.savefig(str(i) + ".png")
    plt.clf()
