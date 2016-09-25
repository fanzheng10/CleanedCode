__author__ = 'fanzheng'

from General import *
import PDB, Analyze
import operator

def createPepFile(shortfile, longfile):
    '''
    :param shortfile: input file with 1-letter peptide sequence
    :param longfile: output file with 3-letter peptide sequence, the last element of each row being a peptide name
    :return:
    '''
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


def clustFPD(inp='score.sc', radius=1):
    '''use the cluster script provided by FlexPePDock, to cluster and generate PDB files
    '''
    os.system('rm c.*.pdb clog')
    nmodels = len(open(inp).readlines())
    cmd = [PATH_fpddemo + '/scripts/clustering/cluster.sh',
           str(nmodels), str(radius), inp, 'native.pdb', 'decoys.silent', 'total_score']
    sub.call(cmd)

def extractRosettaSilent(silent, model):
    cmd = [PATH_rosetta + '/main/source/bin/extract_pdbs.linuxgccrelease',
           '-in:file:silent', silent,
           '-in:file:tags', model]
    sub.call(cmd)

def selectBestModel(name = removePath(os.getcwd())+'.pdb', col = (1, -1)):
    dirs = listDir('.')
    bestsc = 1000.0
    bestdir = ''
    bestname = ''
    for d in dirs:
        scf = d + '/score.sc'
        scores = Analyze.readColumn(scf, col[0], skiprow=1)
        names = Analyze.readColumn(scf, col[1], skiprow=1)
        scores, names = map(float, scores[1:]), names[1:]
        min_ind, min_sc = min(enumerate(scores), key=operator.itemgetter(1))
        min_name = names[min_ind]
        if min_sc < bestsc:
            bestsc, bestdir, bestname = min_sc, d, min_name
    extractRosettaSilent(bestdir+'/decoys.silent', bestname)
    os.system('mv '+bestname +'.pdb '+ name)
    return bestsc, bestdir, bestname