__author__ = 'fanzheng'

from General import *
import Master, conGraph, Terms, Stability, Analyze, Cluster
import pickle

SB = selfbin(sys.argv[0])
par = argparse.ArgumentParser()
par.add_argument('--p', required=True, help='input PDB file')
par.add_argument('--homof', help='the file with homologous information')
args = par.parse_args()

# save input argument
argdict = vars(args)
os.system('rm *.pkl')
pklname = str(os.getpid()) + '.shared.pkl'
pklpath = os.path.realpath(pklname)
pkl = open(pklname, 'w')
pickle.dump(argdict, pkl)
pkl.close()

# get a contact graph
confind_out = changeExt(args.p, 'conf')
if not os.path.isfile(confind_out):
    Master.confind(p=args.p, o=confind_out, rLib=PATH_rotLib)

G = conGraph.conGraph(confind_out)

Homo =  Analyze.findHomo(args.homof)
homos = Homo[getBase(args.p)]

Jobs = []

selfts = []
# search all local fragments
for n in G.nodes():
    t = Terms.Term(parent=args.p, seed=n)
    t.makeFragment(flank=2)
    selfts.append(t)
    crind = t.findResidue(n[0], n[1:])
    rmsd_eff = Stability.rmsdEff(t.getSegLen())
    cmd = ['python', SB + '/EnsemblePreparation.py',
           '--p', t.pdbf,
		   '--homo', ' '.join(homos),
		   '--rmsd', str(rmsd_eff),
		   '--crind', str(crind),
           '--ncon', '0']
    cmd = ' '.join(cmd)
    jobid = n
    job = Cluster.jobOnCluster([cmd], jobid, args.he + '_' + changeExt(t.pdbf, 'seq'))
    Jobs.append(job)
    job.submit(3)
    time.sleep(0.5)

# search all pair fragments
pairts = []
for e in G.edges():
    r1, r2 = e[0], e[1]
    t = Terms.Term(parent=args.p, seed=r1, contact=[r2])
    t.makeFragment(flank=1)
    pairts.append(t)
    crind =  t.findResidue(r1[0], r1[1:]) # should consider two positions in redundancy removal; but worry about that later
    rmsd_eff = Stability.rmsdEff(t.getSegLen())
    cmd = ['python', SB + '/EnsemblePreparation.py',
           '--p', t.pdbf,
		   '--homo', ' '.join(homos),
		   '--rmsd', str(rmsd_eff),
		   '--crind', str(crind),
           '--ncon', '1']
    cmd = ' '.join(cmd)
    jobid = r1 + '_' + r2
    job = Cluster.jobOnCluster([cmd], jobid, args.he + '_' + changeExt(t.pdbf, 'seq'))
    Jobs.append(job)
    job.submit(3)
    time.sleep(0.5)

Cluster.waitJobs(Jobs, giveup_time=1)


# for all possible triple fragments, the criteria of search is at least two sub-pair fragments have more than 100 hits
tripts = []
for n in G.nodes():
    for nn1 in G.neighbors(n):
        for nn2 in G.neighbors(n):
            seqf1 = '_'.join([args.he, getBase(args.p), n, nn1]) + '.seq'
            seqf2 = '_'.join([args.he, getBase(args.p), n, nn2]) + '.seq'
        if (not os.path.isfile(seqf1)) or (not os.path.isfile(seqf2)):
            continue
        else:
            nhit1, nhit12 = 0, 0
            with open(seqf1) as sm1:
                nhit1 = 0
                for lm in sm1:
                    nhit1 +=1
            with open(seqf2) as sm2:
                nhit2 = 0
                for lm in sm2:
                    nhit2 +=1
            if min(nhit1, nhit2) >= args.c2:
                t = Terms.Term(parent=args.p, seed=n, contact=[nn1, nn2])
                t.makeFragment(flank=1)
                tripts.append(t)
                crind =  t.findResidue(r[0], r[1:]) # should consider two positions in redundancy removal; but worry about that later
                rmsd_eff = Stability.rmsdEff(t.getSegLen())
                cmd = ['python', SB + '/EnsemblePreparation.py',
                       '--p', t.pdbf,
                       '--homo', ' '.join(homos),
                       '--rmsd', str(rmsd_eff),
                       '--crind', str(crind),
                       '--ncon', '2']
                cmd = ' '.join(cmd)
                jobid = ' '.join([n, nn1, nn2])
                job = Cluster.jobOnCluster([cmd], jobid, args.he + '_' + changeExt(t.pdbf, 'seq'))
                Jobs.append(job)
                job.submit(3)
                time.sleep(0.5)

Cluster.waitJobs(Jobs, giveup_time=1)

# in this specific code does not consider higher order terms
