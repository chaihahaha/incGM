with open("citeseer/citeseer.cites","r") as f:
    cites = f.read()
edges = cites.split("\n")
edges = [i.split("\t") for i in edges if i != ""]
useddigits = set()
usedstring = set()
for i in range(len(edges)):
    edges[i] = [int(i) if i.isdigit() else i for i in edges[i] ]
for i in range(len(edges)):
    for j in edges[i]:
        if isinstance(j,int):
            useddigits.add(j)
        else:
            usedstring.add(j)
used = list(usedstring|useddigits)
label = [i for i in range(len(used))]
map2int = dict(zip(used,label))
for i in range(len(edges)):
    edges[i] = [map2int[j] for j in edges[i] ]

for i in range(len(edges)):
    edges[i] = "{},{}".format(edges[i][0],edges[i][1])


with open("citeseer/citeseerInt.cites","w") as f:
    f.write("\n".join(edges))
