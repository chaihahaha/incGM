from incGM import FRINGE, incGM_plus
import networkx as nx
import time
import matplotlib.pyplot as plt

def add_edge(f, G, cnt):
    line = f.readline()
    e = tuple(int(i) for i in line.split(","))

    newgraph = nx.Graph()
    if e[0] != e[1]:
        print(cnt, "Adding edge:",e)
        newgraph.add_edge(*e)
    else:
        print(cnt, "Adding node:",e[0])
        newgraph.add_node(e[0])
    print(newgraph.nodes)

    tik = time.time()
    incGM_plus(G, fringe, tau, newgraph)
    tok = time.time()

    print("TOTAL:", tok-tik)
    print("Num of MFS:", len(fringe.MFS))
    print("Num of MIFS:", len(fringe.MIFS))
    print()
    return line

tau = 160
fringe = FRINGE()
G = nx.Graph()
n_edges = 694
with open("citeseerInt.cites","r") as f:
    cnt = 1
    has_next = add_edge(f,G, cnt)
    while has_next:
        cnt += 1
        line = add_edge(f,G, cnt)
        if cnt>n_edges:
            break

distinct = [i for i in fringe.MFS]
for i in range(len(distinct)-1):
    j = i+1
    while 0<=j<len(distinct):
        gi,gj = (G.subgraph(distinct[i]), G.subgraph(distinct[j]))
        if nx.is_isomorphic(gi,gj):
            distinct.pop(j)
            j -= 1
        j += 1
pos = nx.spring_layout(G)
for i in distinct:
    nx.draw_networkx_nodes(G.subgraph(i),pos=pos)
    nx.draw_networkx_labels(G.subgraph(i),pos=pos,labels=dict(zip(i,i)))
    nx.draw_networkx_edges(G.subgraph(i),pos=pos,edge_color='#000000')
    plt.savefig(str(i) + ".png")
    plt.clf()

