__author__ = 'fanzheng'

from General import *
import PDB

def createPepFile(shortfile, longfile):
    '''read file with 1-letter peptide sequences, and create a file with 3-letter sequences,
    with the last element be a peptide name; Good for typical L-sequences'''
    with open(shortfile) as sf, open(longfile, 'w') as lf:
        for l in sf:
            seq1 = l.strip().split()[0]
            info = list(seq1)
            seq3 = []
            for res1 in info:
                res3 = PDB.s2t(res1)
                seq3.append(res3)
            outstr = ' '.join(seq3 + ['|', seq1])
            lf.write(outstr + '\n')


