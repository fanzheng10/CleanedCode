__author__ = 'fanzheng'

from General import *
import networkx as nx

def conGraph(conf, cut = 0.02):
    G = nx.DiGraph()
    with open(conf) as cf:
        for l in cf:
            if not re.match('contact', l):
                continue
            info = l.strip().split()
            res1, res2 = info[1].replace(',', ''), info[2].replace(',', '')
            cond = float(info[3])
            if cond > cut:
                G.add_edge(res1, res2, weight = cond)
                G.node[res1]['aa'] = info[4]
                G.node[res2]['aa'] = info[5]
    return G