__author__ = 'fanzheng'

from General import *
import PDB, conGraph

par = argparse.ArgumentParser()
par.add_argument('--p', nargs= '+', help = 'input structures')
args = par.parse_args()

for p in args.p:

    residues = PDB.ConRes(p)
    G = conGraph.conGraph(changeExt(p, 'conf'), 0.02)
    path = getPath(p)
    selflist = path + '/residue.txt'
    pairlist = path + '/contact.txt'
    with open(selflist, 'w') as sf:
        for r in residues:
            sf.write(r.getChid() + str(r.getResnum()) + '\n')
    with open(pairlist, 'w') as pf:
        edges = G.edges()
        edges = sorted(edges, key = lambda x:(x[0][0], int(x[0][1:]), x[1][0], int(x[1][1:])))
        for e in edges:
            pf.write(e[0] + ' ' + e[1] + '\n')