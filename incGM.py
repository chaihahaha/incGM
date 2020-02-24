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


class FELS:
    def __init__(self,S_nodes, G):
        self.inverted_index = dict() # dict of set of node sets(subgraph)
        self.nodes = tuple(S_nodes)
        self.domain = dict() # upper bound of solution
        self.embeddings = set() # set of tuples
        self.G = G
        self.blacklistnodes = {i:0 for i in G.nodes}
        self.invalid_col = {i:0 for i in self.nodes}

    def add(self, dic):
        emb_nodes =tuple(dic[i] for i in self.nodes)
        self.embeddings.add(emb_nodes)

        # add inverted
        for i in dic.values():
            if i not in self.inverted_index.keys():
                self.inverted_index[i] = set()
            self.inverted_index[i].add(emb_nodes)
        return

    def remove(self, v):
        # remove the outofdate infos of old nodes
        for i in set(G.nodes)-set(self.blacklistnodes.keys()):
            self.blacklistnodes[i] = 0
        if v in self.inverted_index.keys():
            remove_list = list(self.inverted_index[v])
            for emb_nodes in remove_list:
                self.embeddings.remove(emb_nodes)
                for k in self.inverted_index.keys():
                    self.inverted_index[k] -= {emb_nodes}
            self.inverted_index.pop(v)
        return

    def mni(self):
        # generate mni table according to stored embeddings
        table = {i:set() for i in self.nodes} # MNI table lower bound
        for i in range(len(self.nodes)):
            for tup in self.embeddings:
                table[self.nodes[i]].add(tup[i])
        return table

    def frequent(self, tau):
        table = self.mni()
        minsupp = min([len(table[i]) for i in self.nodes])
        return minsupp>=tau

    def infrequent(self, tau):
        maxsupp = min([len(self.domain[i]) for i in self.nodes])
        return maxsupp<tau

    def blacklist(self, dic):
        # add invalid node to blacklist
        for i in dic.values():
            self.blacklistnodes[i] += 1


class FELS_dict:
    def __init__(self):
        self.elements = dict()

    def add(self,G2S_embedding, G):
        # add an embedding
        S2G_embedding = {v: k for k, v in G2S_embedding.items()}
        S_nodes = frozenset(S2G_embedding.keys())
        subG = frozenset(G2S_embedding.keys())
        if S_nodes not in fels_dict.keys():
            self.elements[S_nodes] = FELS(S_nodes, G)
        if subG not in fels_dict.keys():
            self.elements[subG] = FELS(subG, G)
        self.elem(S_nodes).add(S2G_embedding)
        self.elem(subG).add(G2S_embedding)
        return

    def remove(self, S_nodes, v, u):
        # remove from domain to decrease upper bound, push down prunning
        self.domain(S_nodes)[v] -= {u}
        return

    def update_mni(self, v):
        # remove embeddings in case inverted index or mni is outofdate
        for i in self.keys():
            self.elem(i).remove(v)

    def keys(self):
        return self.elements.keys()

    def elem(self, S_nodes):
        return self.elements[S_nodes]

    def mni(self, S_nodes):
        return self.elements[S_nodes].mni()

    def domain(self, S_nodes):
        return self.elements[S_nodes].domain

    def intersection(self, S_nodes, ds, G, tau):
        # intersections of direct subgraphs
        if S_nodes not in self.keys():
            self.elements[S_nodes] = FELS(S_nodes, G)
        for i in S_nodes:
            inter = frozenset(G)
            for j in ds:
                if i in j and j in self.keys():
                    inter &= frozenset(self.mni(j)[i])
            self.domain(S_nodes)[i] = inter
        return

    def is_frequent(self, S_nodes, tau,G):
        # frequentness of lower bound
        if not nx.is_connected(G.subgraph(S_nodes)):
            return False
        return self.elem(S_nodes).frequent(tau)

    def is_infrequent(self, S_nodes, tau):
        # infrequentness of upper bound
        if not nx.is_connected(G.subgraph(S_nodes)):
            return True
        return self.elem(S_nodes).infrequent(tau)

    def iso_graphs(self, S_nodes):
        gs = set()
        inv = self.elem(S_nodes).inverted_index
        for i in inv.keys():
            for j in inv[i]:
                gs.add(j)
        return gs

    def score(self, S_nodes, dic):
        # invalidness of invalid nodes which cannot match with S
        s = 0
        for v in dic.values():
            s += self.elem(S_nodes).blacklistnodes[v]
        return s

    def invalid_col(self, S_nodes, v):
        # add invalid column of mni
        self.elem(S_nodes).invalid_col[v] += 1
    def invalid_score(self, S_nodes, v):
        # invalidness of a column
        return self.elem(S_nodes).invalid_col[v]

def direct_subgraphs(G, S_nodes):
    # direct subgraph missing one node
    ds = []
    for i in S_nodes:
        if S_nodes - {i} and nx.is_connected(G.subgraph(S_nodes - {i})):
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
        fels_dict.add(i, G)
        if fels_dict.is_frequent(S_nodes, tau,G):
            break
    return embeddings

def perm(domain, i, dic, L, mni_table, key_list, G):
    coincide1, coincide2 = False, False
    if i == len(domain.keys()) - 1:
        coincide1 = True
        for k in range(i):
            v = key_list[k]
            if dic[v] not in mni_table[v] :
                coincide1 = False
        coincide2 = True
        for k in range(i):
            v = key_list[k]
            if dic[v] not in key_list:
                coincide2 = False
    if i < len(domain.keys()):
        v = key_list[i]
        # not all mapped nodes in lower bound(added)
        filter1 = set() if not coincide1 else mni_table[v]
        # not all mapped nodes in automorphism(added)
        filter2 = set() if not coincide2 else set(key_list)
        # cannot have duplicate nodes in one embedding, filter old out
        filter3 = set(map(lambda k:dic[k], set(key_list[:i])))
        # new node should not be disconneted to old ones
        filter4 = set(filter(lambda k:not G.subgraph(filter3).edges([k]), domain[v]))
        filtered = domain[v] - filter1 - filter2 - filter3 - filter4
        for j in filtered:
            dic[v] = j
            perm(domain, i+1, dic, L, mni_table, key_list, G)
    else:
        L.append(dic.copy())

def candidates(G, S_nodes, v, u):
    # candidate embeddings to be checked which is generated from domains of MNI
    domain = fels_dict.domain(S_nodes).copy()
    domain[v] = frozenset({u})

    mni_table = fels_dict.mni(S_nodes)

    embeddings = []
    dic = {i:0  for i in S_nodes}
    key_list = sorted(domain.keys())
    perm(domain, 0, dic, embeddings, mni_table, key_list, G)

    # graph node reordering
    embeddings.sort(key=lambda d: fels_dict.score(S_nodes, d))
    return embeddings

def exists_embeddings(G, S_nodes, v, u):
    # check candidate isomorphisms in domain with v->u mapping
    exists = False
    cand = candidates(G, S_nodes, v, u)
    for c in cand:
        iso = True
        S_edges = G.subgraph(S_nodes).edges
        C_edges = G.subgraph(c.values()).edges
        for e in S_edges:
            if len(S_edges) != len(C_edges) or (c[e[0]],c[e[1]]) not in C_edges:
                iso = False
                fels_dict.elem(S_nodes).blacklist(c)
                break
        if iso:
            exists = True
            fels_dict.add(c, G)
    return exists

def EVALUATE(G, tau, S_nodes):
    if not nx.is_connected(G.subgraph(S_nodes)):
        return False
    ds = direct_subgraphs(G, S_nodes)
    fels_dict.intersection(S_nodes, ds, G, tau) # push down prunning
    # Automorphisms
    gm = isomorphism.GraphMatcher(G.subgraph(S_nodes), G.subgraph(S_nodes))
    for automorphism in gm.isomorphisms_iter():
        fels_dict.add(automorphism, G)
    if fels_dict.is_frequent(S_nodes, tau, G):
        return True

    # MNI col reordering
    s_node_list = list(S_nodes)
    s_node_list.sort(key=lambda v: fels_dict.invalid_score(S_nodes, v), reverse=True)

    for v in s_node_list:
        count = 0
        if fels_dict.is_infrequent(S_nodes, tau):
            return False
        if fels_dict.is_frequent(S_nodes, tau, G):
            return True

        # Lazy search ??
        for u in fels_dict.domain(S_nodes)[v]:
            if u in fels_dict.mni(S_nodes)[v]:
                count += 1
            else:
                exists = exists_embeddings(G, S_nodes, v, u)
                if exists:
                    count += 1
                else:
                    fels_dict.remove(S_nodes, v, u)
            if count == tau:
                break
        else:
            fels_dict.invalid_col(S_nodes, v)
            return False
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
    for v in newnodes:
        fels_dict.update_mni(v)
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
nx.draw_networkx_nodes(base,pos=pos,node_color='#000000')
nx.draw_networkx_edges(base,pos=pos,edge_color='#000000')
plt.savefig("base.png")
plt.clf()
G = nx.Graph()
tau = 7
fringe = FRINGE()
cnt = 0
for e in base.edges:
    cnt += 1
    print(cnt, "adding:",e)
    tik = time.time()
    incGM_plus(G,fringe,tau,base.subgraph(e))
    tok = time.time()
    print("TOTAL:",tok-tik)
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
