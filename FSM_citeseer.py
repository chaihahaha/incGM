from incGM import FRINGE, incGM_plus
import networkx as nx
import time

def add_edge(f, G, cnt):
    line = f.readline()
    e = tuple(int(i) for i in line.split(","))

    print(cnt, "Adding edge:",e)
    tik = time.time()
    incGM_plus(G, fringe, tau, nx.Graph([e]))
    tok = time.time()

    print("TOTAL:", tok-tik)
    print("Num of MFS:", len(fringe.MFS))
    print()
    return line

tau = 20
fringe = FRINGE()
G = nx.Graph()
with open("citeseerInt.cites","r") as f:
    cnt = 1
    has_next = add_edge(f,G, cnt)
    while has_next:
        cnt += 1
        line = add_edge(f,G, cnt)


