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
        self.table = dict() # lower bound of solution (marked solutions)
        self.domain = dict() # upper bound of solution
        for i in self.nodes:
            self.table[i]=set()
    def add(self, dic):
        for i in self.nodes:
            self.table[i].add(dic[i])
    def frequent(self, tau):
        minsupp = min([len(self.table[i]) for i in self.nodes])
        return minsupp>=tau
    def infrequent(self, tau):
        maxsupp = min([len(self.domain[i]) for i in self.nodes])
        return maxsupp<tau


class FELS:
    def __init__(self,S_nodes):
        self.S_nodes = S_nodes
        self.mni = MNI_table(self.S_nodes)
        self.embeddings = set()
        self.inverted_index = dict() # dict of set of node sets(subgraph)
    def add(self, embedding):
        self.embeddings.add(frozenset(embedding.values()))
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

    def remove(self, S_nodes, v, u):
        self.domain(S_nodes)[v] -= {u}
        return

    def keys(self):
        return self.elements.keys()

    def elem(self, S_nodes):
        return self.elements[S_nodes]

    def mni(self, S_nodes):
        return self.elements[S_nodes].mni.table

    def domain(self, S_nodes):
        return self.elements[S_nodes].mni.domain

    def intersection(self, S_nodes, ds, G):
        if S_nodes not in self.keys():
            self.elements[S_nodes] = FELS(S_nodes)
        for i in S_nodes:
            inter = frozenset(G)
            for j in ds:
                if i in j and j in self.keys():
                    inter &= frozenset(self.mni(j)[i])
            self.domain(S_nodes)[i] = inter - frozenset(self.mni(S_nodes)[i])
        return

    def is_frequent(self, S_nodes, tau,G):
        if not nx.is_connected(G.subgraph(S_nodes)):
            return False
        return self.elem(S_nodes).mni.frequent(tau)

    def is_infrequent(self, S_nodes, tau):
        if not nx.is_connected(G.subgraph(S_nodes)):
            return True
        return self.elem(S_nodes).mni.infrequent(tau)

    def iso_graphs(self, S_nodes):
        gs = set()
        inv = self.elem(S_nodes).inverted_index
        for i in inv.keys():
            for j in inv[i]:
                gs.add(j)
        return gs

def direct_subgraphs(G, S_nodes):
    ds = []
    for i in S_nodes:
        if nx.is_connected(G.subgraph(S_nodes - {i})):
            ds.append(S_nodes - {i})
    return ds

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

def candidates(G, S_nodes, v, u):
        # candidate embeddings to be checked which is generated from domains of MNI
        embeddings = {frozenset({u})}
        for i in (S_nodes-{v}):
            newembeddings = set()
            for j in fels_dict.domain(S_nodes)[i]:
                for k in embeddings:
                    union = k | {j}
                    if nx.is_connected(G.subgraph(union)):
                        newembeddings.add(union)
            embeddings = newembeddings
        return embeddings

def exists_embeddings(G, S_nodes, v, u):
    exists = False
    frequent = False
    cand = candidates(G, S_nodes, v, u)
    for c in cand:
        gm =  isomorphism.GraphMatcher(G.subgraph(S_nodes), G.subgraph(c))
        for e in gm.isomorphisms_iter():
            fels_dict.add(e)
            if e[v] == u:
                exists = True
            if fels_dict.is_frequent(S_nodes, tau, G):
                frequent = True
                break
        else:
            continue
        break
    return exists, frequent

def EVALUATE(G, tau, S_nodes):
    if not nx.is_connected(G.subgraph(S_nodes)):
        return False
    ds = direct_subgraphs(G, S_nodes)
    fels_dict.intersection(S_nodes, ds, G)
    for v in S_nodes:
        count = 0
        if fels_dict.is_infrequent(S_nodes, tau):
            return False
        for u in fels_dict.domain(S_nodes)[v]:
            if u in fels_dict.mni(S_nodes)[v]:
                count += 1
            else:
                exists, frequent = exists_embeddings(G, S_nodes, v, u)
                if frequent:
                    return True
                if exists:
                    count += 1
                else:
                    fels_dict.remove(S_nodes, v, u)
            if count == tau:
                break
        else:
            return False
    if fels_dict.is_frequent(S_nodes, tau, G) == False:
        print(fels_dict.mni(S_nodes))
    assert fels_dict.is_frequent(S_nodes, tau, G) == True
    return True

def UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G):
    deleted = False
    if isFreq:
        added = fringe.addMFS(S_nodes)
        deleted = fringe.removeMIFS(S_nodes)

        for i in range(len(fringe.MFS)):
            MFSi = fringe.MFS[i]
            if len(MFSi) == len(S_nodes) and MFSi != S_nodes:
                u = MFSi | S_nodes
                if u not in fringe.MIFS and u not in fringe.MFS:
                    if not EVALUATE(G,tau,u):
                        joined = fringe.addMIFS(u)
    return deleted

fels_dict = FELS_dict()

def incGM_plus(G, fringe, tau, newgraph):
    G.add_edges_from(newgraph.edges)
    newnodes = frozenset(newgraph.nodes)
    fringe.addMIFS(newnodes)
    i = 0
    sl = 0
    ev = 0
    while 0 <= i <len(fringe.MIFS):
        S_nodes = fringe.MIFS[i]
        sl_tik = time.time()
        embeds = SEARCHLIMITED(S_nodes, newnodes,G)
        sl_tok = time.time()
        sl += sl_tok - sl_tik
        print("search limited:", sl)
        if not embeds:
            ev_tik = time.time()
            isFreq = EVALUATE(G,tau,S_nodes)
            ev_tok = time.time()
            ev += ev_tok - ev_tik
            print("evaluate:",ev)
        else:
            isFreq = fels_dict.is_frequent(S_nodes, tau, G)

        ev_tik = time.time()
        delete = UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G)
        ev_tok = time.time()
        ev += ev_tok - ev_tik
        print("evaluate:",ev)
        i = i + 1 - int(delete)
    return fringe.MFS

base = nx.gnm_random_graph(15,25,1)
pos = nx.spring_layout(base)
nx.draw_networkx_nodes(base,pos=pos,node_color='#000000')
nx.draw_networkx_edges(base,pos=pos,edge_color='#000000')
plt.savefig("base.png")
plt.clf()
G = nx.Graph()
tau = 7
fringe = FRINGE()
for e in base.edges:
    tik = time.time()
    incGM_plus(G,fringe,tau,base.subgraph(e))
    tok = time.time()
    print("TOTAL:",tok-tik)
    print("=================")
    print()
    print()
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
    nx.draw_networkx_nodes(G.subgraph(i),pos=pos,node_color='#000000')
    nx.draw_networkx_edges(G.subgraph(i),pos=pos,edge_color='#000000')
    plt.savefig(str(i) + ".png")
    plt.clf()
