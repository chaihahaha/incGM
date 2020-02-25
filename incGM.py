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
        self.domain = {i:frozenset(G) for i in S_nodes} # upper bound of solution
        self.mni = {i:dict() for i in S_nodes} # upper bound of solution
        self.embeddings = set() # set of tuples
        self.G = G
        self.blacklistnodes = {i:0 for i in G.nodes}
        self.blacklistedges = {i:0 for i in G.edges}
        self.invalid_col = {i:0 for i in self.nodes}

    def add(self, dic):
        emb_nodes =tuple(dic[i] for i in self.nodes)

        # add to embeddings
        self.embeddings.add(emb_nodes)

        # add to mni
        for i in self.nodes:
            if dic[i] not in self.mni[i]:
                self.mni[i][dic[i]] = 1
            else:
                self.mni[i][dic[i]] += 1

        # add inverted
        for i in dic.values():
            if i not in self.inverted_index.keys():
                self.inverted_index[i] = set()
            self.inverted_index[i].add(emb_nodes)

        # add to domain
        for i in self.nodes:
            mni_col_i = self.mni[i]
            self.domain[i] |= set([j for j in mni_col_i.keys() if mni_col_i[j]>0])
        return

    def remove(self, v):
        # remove the outofdate infos of old nodes
        if v in self.inverted_index.keys():
            remove_list = list(self.inverted_index[v])
            for emb_nodes in remove_list:
                # remove from embeddings
                self.embeddings.remove(emb_nodes)

                # remove from mni
                emb = dict(zip(self.nodes, emb_nodes))
                for i in self.nodes:
                    mni_col_i = self.mni[i]
                    if emb[i] in mni_col_i.keys():
                        if mni_col_i[emb[i]] > 0:
                            mni_col_i[emb[i]] -= 1

                # remove from inverted index
                for k in self.inverted_index.keys():
                    self.inverted_index[k] -= {emb_nodes}
            self.inverted_index.pop(v)
        return

    def update_domain(self, v):
        for i in self.nodes:
            self.domain[i] |= {v}

    def mni_valid(self):
        mni = self.mni
        mni_val = {k:{i for i in mni[k].keys() if mni[k][i]>0} for k in self.nodes}
        return mni_val

    def frequent(self, tau):
        mni = self.mni
        minsupp = min([len([i for i in mni[k].keys() if mni[k][i]>0]) for k in self.nodes])
        #print(self.mni,minsupp>=tau)
        return minsupp>=tau

    def infrequent(self, tau):
        maxsupp = min([len(self.domain[i]) for i in self.nodes])
        return maxsupp<tau


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
            self.elem(i).update_domain(v)

    def keys(self):
        return self.elements.keys()

    def elem(self, S_nodes):
        return self.elements[S_nodes]

    def mni(self, S_nodes):
        return self.elements[S_nodes].mni_valid()

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
                    inter &= self.domain(j)[i]
            self.domain(S_nodes)[i] &= inter
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

    def valid_node_edge(self, S_nodes, dic, e):
        # add invalid node to blacklist
        invalid_nodes = self.elem(S_nodes).blacklistnodes
        invalid_edges = self.elem(S_nodes).blacklistedges
        for i in dic.values():
            if i not in invalid_nodes.keys():
                invalid_nodes[i] = 0
            invalid_nodes[i] -= 1
        if e not in invalid_edges.keys():
            invalid_edges[e] = 0
        invalid_edges[e] -= 1

    def invalid_node_edge(self, S_nodes, dic, e):
        # add invalid node to blacklist
        invalid_nodes = self.elem(S_nodes).blacklistnodes
        invalid_edges = self.elem(S_nodes).blacklistedges
        for i in dic.values():
            if i not in invalid_nodes.keys():
                invalid_nodes[i] = 0
            invalid_nodes[i] += 1
        if e not in invalid_edges.keys():
            invalid_edges[e] = 0
        invalid_edges[e] += 1

    def invalid_node_score(self, S_nodes, v):
        # invalidness of invalid nodes which cannot match with S
        invalid_nodes = self.elem(S_nodes).blacklistnodes
        if v not in invalid_nodes.keys():
            invalid_nodes[v] = 0
        return invalid_nodes[v]

    def invalid_edge_score(self, S_nodes, e):
        # invalidness of invalid nodes which cannot match with S
        invalid_edges = self.elem(S_nodes).blacklistedges
        if e not in invalid_edges.keys():
            invalid_edges[e] = 0
        return invalid_edges[e]

    def invalid_col(self, S_nodes, v):
        # add invalid column of mni
        self.elem(S_nodes).invalid_col[v] += 1
    def invalid_col_score(self, S_nodes, v):
        # invalidness of a column
        return self.elem(S_nodes).invalid_col[v]

def direct_subgraphs(G, S_nodes):
    # direct subgraph missing one node
    ds = []
    for i in S_nodes:
        if S_nodes - {i} and nx.is_connected(G.subgraph(S_nodes - {i})):
            ds.append(S_nodes - {i})
    return ds

def candidates(G, S_nodes, v, u):
    # return candidate embeddings generator
    # candidate embeddings to be checked which is generated from domains of MNI
    domain = fels_dict.domain(S_nodes).copy()
    domain[v] = frozenset({u})

    mni_table = fels_dict.mni(S_nodes)

    dic = {i:0  for i in S_nodes}

    # MNI column reordering
    key_list = sorted(list(domain.keys()), key=lambda v: fels_dict.invalid_col_score(S_nodes, v), reverse=True)
    return perm_limited(domain, 0, dic, mni_table, key_list, G, set(G.nodes))

def exists_embeddings(G, S_nodes, v, u, cand, exhastive):
    # check whether where exists an isomorphism in domain with v->u mapping
    # return ifexists, False if not timeout
    # return generator, True if timeout
    if not cand:
        cand = candidates(G, S_nodes, v, u)

    tik = time.time()
    for c in cand:
        iso = True
        S_edges = list(G.subgraph(S_nodes).edges)
        C_edges = G.subgraph(c.values()).edges
        for e in S_edges:
            if len(S_edges) != len(C_edges) or (c[e[0]],c[e[1]]) not in C_edges:
                iso = False
                fels_dict.invalid_node_edge(S_nodes, c, e)
                break
        if iso:
            fels_dict.add(c, G)
            return True, False
        if time.time()-tik>0.5 and not exhastive:
            return cand, True
    return False, False

def perm_limited(domain, i, dic, mni_table, key_list, G, limitation):
    # yield all the possible unevaluated embeddings in domain
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
    lim = set(G.nodes)
    if i>len(key_list)-len(limitation)-1 and not (set(key_list[:i]) & limitation):
        lim = limitation
    if i < len(domain.keys()):
        v = key_list[i]
        # not all mapped nodes in lower bound(added)
        filter1 = set() if not coincide1 else mni_table[v]
        # not all mapped nodes in automorphism(added)
        filter2 = set() if not coincide2 else set(key_list)
        # cannot have duplicate nodes in one embedding, filter old out
        filter3 = set(map(lambda k:dic[k], set(key_list[:i])))
        # last node should not be disconneted to old ones
        filter4 = set() if not i == len(domain.keys())-1 else set(filter(lambda k:not nx.is_connected(G.subgraph(filter3|{k})), domain[v]))
        filtered = (domain[v] - filter1 - filter2 - filter3 - filter4)&lim
        for j in filtered:
            dic[v] = j
            yield from perm_limited(domain, i+1, dic, mni_table, key_list, G, limitation)
    else:
        yield dic

def candidates_limited(G, S_nodes, limitation):
    # return candidate embeddings generator
    # candidate embeddings to be checked which is generated from domains of MNI
    domain = fels_dict.domain(S_nodes).copy()

    mni_table = fels_dict.mni(S_nodes)

    dic = {i:0  for i in S_nodes}

    key_list = list(domain.keys())
    return perm_limited(domain, 0, dic, mni_table, key_list, G, limitation)

def exists_embeddings_limited(G, S_nodes, limitation):
    # check whether where exists an isomorphism in domain with limitaiton in mapped graph
    # return ifexists
    cand = candidates_limited(G, S_nodes, limitation)

    for c in cand:
        iso = True
        S_edges = list(G.subgraph(S_nodes).edges)
        C_edges = G.subgraph(c.values()).edges
        for e in S_edges:
            if len(S_edges) != len(C_edges) or (c[e[0]],c[e[1]]) not in C_edges:
                iso = False
                fels_dict.invalid_node_edge(S_nodes, c, e)
                break
        if iso:
            fels_dict.add(c, G)
            return True
    return False

def EVALUATE(G, tau, S_nodes):
    if not nx.is_connected(G.subgraph(S_nodes)):
        return False
    ds = direct_subgraphs(G, S_nodes)
    fels_dict.intersection(S_nodes, ds, G, tau) # push down prunning
    # Automorphisms
    gm = isomorphism.GraphMatcher(G.subgraph(S_nodes), G.subgraph(S_nodes))
    for automorphism in gm.isomorphisms_iter():
        fels_dict.add(automorphism, G)


    # MNI column reordering
    s_node_list = list(S_nodes)
    s_node_list.sort(key=lambda v: fels_dict.invalid_col_score(S_nodes, v), reverse=True)

    for v in s_node_list:
        # count # of the valid rows in column v
        count = 0

        # store search state when timed out
        timedout_search = []

        if fels_dict.is_infrequent(S_nodes, tau):
            return False
        if fels_dict.is_frequent(S_nodes, tau, G):
            return True

        # Lazy search
        # graph node reordering
        d_v = sorted(list(fels_dict.domain(S_nodes)[v]), key=lambda d: fels_dict.invalid_node_score(S_nodes, d))
        for u in d_v:
            if u in fels_dict.mni(S_nodes)[v]:
                count += 1
            else:
                exists, timeout = exists_embeddings(G, S_nodes, v, u, None, False)
                if timeout:
                    timedout_search.append((exists, u))
                else:
                    if exists:
                        count += 1
                    else:
                        fels_dict.remove(S_nodes, v, u)
            if count == tau:
                break
        else:
            # Resume timed-out search if needed
            if len(timedout_search) + count >=tau:
                for (cand, u) in timedout_search:
                    exists, timeout = exists_embeddings(G, S_nodes, v, u, cand, True)
                    assert timeout==False
                    if exists:
                        count += 1
                    else:
                        fels_dict.remove(S_nodes, v, u)
                    if count == tau:
                        break
            fels_dict.invalid_col(S_nodes, v)
            return False
    assert fels_dict.is_frequent(S_nodes, tau, G) == True
    return True

def SEARCHLIMITED(S_nodes,newnodes,G):
    if not nx.is_connected(G.subgraph(S_nodes)):
        return False
    ds = direct_subgraphs(G, S_nodes)
    fels_dict.intersection(S_nodes, ds, G, tau) # push down prunning

    if fels_dict.is_infrequent(S_nodes, tau):
        return True
    if fels_dict.is_frequent(S_nodes, tau, G):
        return True

    # limited search
    exists = exists_embeddings_limited(G, S_nodes, newnodes)
    return exists

def UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G, add, delete):
    if isFreq:
        fringe.addMFS(S_nodes)
        delete.append(S_nodes)

        for i in range(len(fringe.MFS)):
            MFSi = fringe.MFS[i]
            if len(MFSi) == len(S_nodes) and MFSi != S_nodes:
                u = MFSi | S_nodes
                if u not in fringe.MIFS + add and u not in fringe.MFS:
                    if not EVALUATE(G,tau,u):
                        add.append(u)
                    else:
                        fringe.addMFS(u)
                    fels_dict.intersection(u, [MFSi, S_nodes], G, tau)

    return

fels_dict = FELS_dict()

def incGM_plus(G, fringe, tau, newgraph):
    G.add_edges_from(newgraph.edges)
    newnodes = frozenset(newgraph.nodes)
    for v in newnodes:
        fels_dict.update_mni(v)
    fringe.addMIFS(newnodes)
    delete, add = [], []
    for S_nodes in fringe.MIFS:
        embeds = SEARCHLIMITED(S_nodes, newnodes,G)
        if not embeds:
            isFreq = EVALUATE(G,tau,S_nodes)
        else:
            isFreq = fels_dict.is_frequent(S_nodes, tau, G)
        UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G, add, delete)
    for i in delete:
        fringe.removeMIFS(i)
    for i in add:
        fringe.addMIFS(i)
    return fringe.MFS

base = nx.gnm_random_graph(5,9,1)
pos = nx.spring_layout(base)
nx.draw_networkx_nodes(base,pos=pos,node_color='#000000')
nx.draw_networkx_edges(base,pos=pos,edge_color='#000000')
plt.savefig("base.png")
plt.clf()
G = nx.Graph()
tau = 3
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
