import networkx as nx
from networkx.algorithms import isomorphism
import matplotlib.pyplot as plt
import time
from ortools.sat.python import cp_model
from itertools import combinations
class VarArraySolutionFinder(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    def __init__(self, X, S_nodes, fels_dict, G):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__solution_count = 0
        self.X = X
        self.S_nodes = S_nodes
        self.fels_dict = fels_dict
        self.G = G

    def on_solution_callback(self):
        self.__solution_count += 1
        self.fels_dict.add({i: self.Value(self.X[i]) for i in self.S_nodes},self.G)

    def solution_count(self):
        return self.__solution_count

def var_from_domain(model, name, domain):
    "initialize a variable with integer domain defined by domain"
    domain = cp_model.Domain.FromIntervals([[i] for i in domain])
    val = model.NewIntVarFromDomain(domain, name)
    return val

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

    def is_infrequent(self, S_nodes, tau, G):
        # infrequentness of upper bound
        if not nx.is_connected(G.subgraph(S_nodes)):
            return True
        return self.elem(S_nodes).infrequent(tau)

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

def exists_embeddings(G, S_nodes, v, u):
    model = cp_model.CpModel()
    S = G.subgraph(S_nodes)
    X = {i:var_from_domain(model, "node"+str(i),fels_dict.domain(S_nodes) if i!=v else [u]) for i in S_nodes}

    model.AddAllDifferent(X.values())

    # S -> subgraph of G
    GEdges = list(G.edges) + [e[::-1] for e in G.edges]
    for nodei, nodej in combinations(S_nodes,2):
        if (nodei,nodej) in S.edges:
            model.AddAllowedAssignments((X[nodei],X[nodej]), GEdges)
        else:
            model.AddForbiddenAssignments((X[nodei],X[nodej]), GEdges)

    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    if status == cp_model.FEASIBLE:
        fels_dict.add({i: solver.Value(X[i]) for i in S_nodes},G)
        return True
    else:
        return False

def exists_embeddings_limited(G, S_nodes, subset):
    model = cp_model.CpModel()
    S = G.subgraph(S_nodes)
    X = {i:var_from_domain(model, "node"+str(i),fels_dict.domain(S_nodes)) for i in S_nodes}
    model.AddAllDifferent(X.values())
    # S -> subgraph of G
    GEdges = list(G.edges) + [e[::-1] for e in G.edges]
    for nodei, nodej in combinations(S_nodes,2):
        if (nodei,nodej) in S.edges:
            model.AddAllowedAssignments((X[nodei],X[nodej]), GEdges)
        else:
            model.AddForbiddenAssignments((X[nodei],X[nodej]), GEdges)
    # subset <= X
    for element in subset:
        x_eq_subset = [model.NewBoolVar("node"+str(k)+"=="+str(element)) for k in S_nodes]
        for x,e in zip(X.values(),x_eq_subset):
            # e <=> x == element
            model.Add(x==element).OnlyEnforceIf(e)
            model.Add(x!=element).OnlyEnforceIf(e.Not())
        # element is in X <=> exists x which equals to element
        model.AddBoolOr(x_eq_subset)
    solver = cp_model.CpSolver()
    solution_finder = VarArraySolutionFinder(X, S_nodes, fels_dict, G)
    status = solver.SearchForAllSolutions(model, solution_finder)
    if solution_finder.solution_count()>0:
        return True
    else:
        return False


def EVALUATE(G, tau, S_nodes):
    S = G.subgraph(S_nodes)
    if not S or not nx.is_connected(S):
        return False
    ds = direct_subgraphs(G, S_nodes)
    fels_dict.intersection(S_nodes, ds, G, tau) # push down prunning
    # Automorphisms
    gm = isomorphism.GraphMatcher(S, S)
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

        if fels_dict.is_infrequent(S_nodes, tau, G):
            return False
        if fels_dict.is_frequent(S_nodes, tau, G):
            return True

        # Lazy search
        # graph node reordering
        d_v = list(fels_dict.domain(S_nodes)[v])
        for u in d_v:
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

def SEARCHLIMITED(S_nodes,newnodes, tau, G):
    S = G.subgraph(S_nodes)
    if not S or not nx.is_connected(S):
        return False
    ds = direct_subgraphs(G, S_nodes)
    fels_dict.intersection(S_nodes, ds, G, tau) # push down prunning

    if fels_dict.is_infrequent(S_nodes, tau, G):
        return True
    if fels_dict.is_frequent(S_nodes, tau, G):
        return True

    # limited search
    exists = exists_embeddings_limited(G, S_nodes, newnodes)
    return exists

def UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G):
    deleted = False
    if isFreq:
        fringe.addMFS(S_nodes)
        deleted = fringe.removeMIFS(S_nodes)

        for i in range(len(fringe.MFS)):
            MFSi = fringe.MFS[i]
            if len(MFSi) == len(S_nodes) and MFSi != S_nodes:
                u = MFSi | S_nodes
                if u not in fringe.MIFS and u not in fringe.MFS:
                    fringe.addMIFS(u)
                    fels_dict.intersection(u, [MFSi, S_nodes], G, tau)
    return deleted

fels_dict = FELS_dict()

def incGM_plus(G, fringe, tau, newgraph):
    G.add_edges_from(newgraph.edges)
    newnodes = frozenset(newgraph.nodes)
    for v in newnodes:
        fels_dict.update_mni(v)
    fringe.addMIFS(newnodes)
    i = 0
    while 0 <= i < len(fringe.MIFS):
        S_nodes = fringe.MIFS[i]
        embeds = SEARCHLIMITED(S_nodes, newnodes,tau,G)
        if not embeds:
            isFreq = EVALUATE(G,tau,S_nodes)
        else:
            isFreq = fels_dict.is_frequent(S_nodes, tau, G)
        delete = UPDATEFRINGE(fringe, S_nodes, isFreq, tau, G)
        i = i + 1 - int(delete)
    return fringe.MFS

if __name__=="__main__":
    base = nx.gnm_random_graph(15,25,1)
    pos = nx.spring_layout(base)
    nx.draw_networkx_nodes(base,pos=pos)
    nx.draw_networkx_labels(base,pos=pos,labels=dict(zip(base.nodes,base.nodes)))
    nx.draw_networkx_edges(base,pos=pos,edge_color='#000000')
    plt.savefig("base.png")
    plt.clf()
    G = nx.Graph()
    tau = 4
    fringe = FRINGE()
    cnt = 0
    for e in base.edges:
        cnt += 1
        print(cnt, "adding:",e)
        tik = time.time()
        incGM_plus(G,fringe,tau,base.subgraph(e))
        tok = time.time()
        print("TOTAL:",tok-tik)
    print("Num of MFS:",len(fringe.MFS))
    #if fringe.MFS:
    #    n1,n2 = len(fringe.MFS)-4, len(fringe.MFS)
    #    for i in range(n1,n2):
    #        print("Embeddings of ",fels_dict.elem(fringe.MFS[i]).nodes,":")
    #        print(fels_dict.elem(fringe.MFS[i]).embeddings)
    #        print()
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
        nx.draw_networkx_nodes(G.subgraph(i),pos=pos)
        nx.draw_networkx_labels(G.subgraph(i),pos=pos,labels=dict(zip(i,i)))
        nx.draw_networkx_edges(G.subgraph(i),pos=pos,edge_color='#000000')
        plt.savefig(str(i) + ".png")
        plt.clf()
