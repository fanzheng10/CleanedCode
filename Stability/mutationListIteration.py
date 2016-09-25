__author__ = 'fanzheng'

from General import *
import Terms, Master, Stability, Cluster
import pickle, itertools

SB = selfbin(sys.argv[0])

par = argparse.ArgumentParser()
par.add_argument('--l', required = True, help = 'a list of mutations')
par.add_argument('--i', default = os.getcwd(), help = 'the directory to find the input structures')
par.add_argument('--db', default = '/home/anthill/fzheng/ironfs/searchDB/bc-30-sc-20141022-newpds', help = 'the directory of the searching database')
par.add_argument('--dbl', help = 'list of the searching database')
par.add_argument('--he', help = 'the header of the output files')
par.add_argument('--homof', help='the file with homologous information')
par.add_argument('--c1', default = 20000, type = int, help = 'cutoff for top N matches')
par.add_argument('--c2', default = 100, type = int, help = 'the criteria of increasing the order of sub-TERMs')
par.add_argument('--c3', default = 3, type = int, help = 'the maximum order of sub-TERMs to consider')
par.add_argument('--nr', default=0.4, type=float, help='the level of redundancy removal')
par.add_argument('--maxcon', default=2, type=int, help='the maximum contact an ensemble can have')
args = par.parse_args()

# save input argument
argdict = vars(args)
os.system('rm *.pkl')
pklname = str(os.getpid()) + '.shared.pkl'
pklpath = os.path.realpath(pklname)
pkl = open(pklname, 'w')
pickle.dump(argdict, pkl)
pkl.close()

# a few things should be keep the same:
# the name of the pdb file (without .pdb extension)
# first column of the mutation file
# first column of the homologous file
# there should not be '_' in these

# read the input file and create a mutation list
pdbs = []
positions =[]
with open(args.l) as mutf:
	for l in mutf:
		mut = Stability.getMutation(l.strip())
		pdb, pos = mut.p, mut.dir
		if not pdb in pdbs:
			pdbs.append(pdb)
		if not pos in positions:
			positions.append(pos)


# run confind for all structures
for pdb in pdbs:
	pdbf = args.i + '/' + pdb + '.pdb'
	assert os.path.isfile(pdbf), 'the pdb file '+pdbf +' does not exist; quit...'
	confind_out = changeExt(pdbf, 'conf')
	if os.path.isfile(confind_out):
		continue
	else: # run confind
		Master.confind(p=pdbf, o=confind_out, rLib=PATH_rotLib)


# find homologs
Homo = {}
if args.homof != None:
	with open(args.homof) as hf:
		for hl in hf:
			items = hl.strip().split()
			Homo[items[0]] = items[1:]

# contact identification
pos2cons = {}
pos2pdb = {}

odir = os.getcwd()
for pos in positions:
	os.chdir(odir)
	if not os.path.isdir(pos):
		os.mkdir(pos)
	os.chdir(pos)
	pid, ipos = pos.split('_')
	icid, iresnum = ipos[0], ''.join(ipos[1:])
	pcid = pid + '_' + icid
	conListf = pos + '.conlist'
	cons, _ = Terms.contactList(profile=args.i+'/'+pid+'.conf',
								resnum=iresnum, cid=icid,
								outFile=conListf,
								dmin=0.02, monomer=True)
	pos2pdb[pos] = args.i+'/'+pid+'.pdb'
	pos2cons[pos] = cons
os.chdir(odir)


# recursive fragmentation and master search
jobs = {}

for pos in positions:
	os.chdir(odir)
	os.chdir(pos)
	iTerms = []
	parentpdb = pos2pdb[pos]
	cons = pos2cons[pos]
	seed = pos.split('_')[1]
	pcid = pos[0:5]
	if pcid in Homo:
		homos = Homo[pcid]
	else:
		homos = []

	# create self term without contact
	selfterm = Terms.Term(parent=parentpdb,
					   seed=seed)
	selfterm.makeFragment(flank=2)
	iTerms.append(selfterm)

	# create pair term with 1 contact
	for con in cons:
		c1term =Terms.Term(parent=parentpdb,
						   seed=seed,
						   contact=[con])
		c1term.makeFragment(flank=1)
		iTerms.append(c1term)

	# submit MASTER search
	pdbs, rmsds, crinds = [], [], []

	for tm in iTerms: #
		crind = tm.findResidue(seed[0], seed[1:]) # index of the seed residue in the corresponding term
		rmsd_eff = Stability.rmsdEff(tm.getSegLen(), 18, 1.0)
		pdbs.append(tm.pdbf)
		rmsds.append(str(rmsd_eff))
		crinds.append(str(crind))

	cmd = ['python', SB + '/EnsemblePreparation.py',
		   '--pkl', pklpath,
		   '--p', ' '.join(pdbs),
		   '--homo', ' '.join(homos),
		   '--rmsd', ' '.join(rmsds),
		   '--crind', ' '.join(crinds),
		   '--ncon', '1']
	cmd = ' '.join(cmd)

	# create job object
	jobid = pos + '.1'
	job = Cluster.jobOnCluster([cmd], jobid, odir+'/'+pos + '/.finished.1')
	jobs[jobid] = job
	job.submit(3)
	time.sleep(0.5)

os.chdir(odir)

# now considering higher-order fragments
Ncon = 2

while Ncon <= args.maxcon:
	positions_copy = [x for x in positions]
	while len(positions_copy) > 0:
		for pos in positions_copy:
			time.sleep(1)
			os.chdir(odir)
			os.chdir(pos)
			parentpdb = pos2pdb[pos]
			seed = pos.split('_')[1]

			jobid = pos+'.'+str(Ncon)
			# check if the lower order job is still running
			low_jobid = pos+'.'+str(Ncon-1)
			if not low_jobid in jobs:
				positions_copy.remove(pos)
				continue
			lowj = jobs[low_jobid]
			lowj.checkjob()
			if lowj.running == 0:
				jobs.pop(low_jobid, None)
				if lowj.checkfinish() == 0:
					print('Job about ' + lowj.myid + ' may have died ...')
					if lowj.tried >1:
						print(lowj.myid + ' has failed once, give up ...')
						positions_copy.remove(pos)
						continue
					print('resubmitting ... ')
					os.chdir(pos)
					lowj.submit(12)
					jobs[lowj.myid] = lowj
				else: # if the low job has finished, prepare the current level
					positions_copy.remove(pos)
					# make some new fragments
					cons = pos2cons[pos]
					contuples = itertools.combinations(cons, Ncon)
					morepdbs, morermsds, morecrinds = [], [], []
					seed = pos.split('_')[1]
					for set in contuples:
						sortset = sorted(set, key=lambda x:(x[0], int(x[1:])))
						if Stability.whetherIncludeHighOrder(pos, args.he, sortset, minmatch=args.c2):
							Hterm =Terms.Term(parent=parentpdb,
								   seed=seed,
								   contact=sortset)
							Hterm.makeFragment(flank=1)
							morepdbs.append(Hterm.pdbf)
							morermsds.append(str(Stability.rmsdEff(Hterm.getSegLen(), 18, 1.0)))
							morecrinds.append(str(Hterm.findResidue(seed[0], seed[1:])))

					cmd = ['python', SB + '/EnsemblePreparation.py',
		   					'--pkl', pklpath,
		   					'--p', ' '.join(morepdbs),
							'--homo', ' '.join(homos),
							'--rmsd', ' '.join(morermsds),
		   					'--crind', ' '.join(morecrinds),
		   					'--ncon', str(Ncon)]
					cmd = ' '.join(cmd)

					new_jobid = pos+'.'+str(Ncon)
					job = Cluster.jobOnCluster([cmd], new_jobid, odir+'/'+pos + '/.finished.'+str(Ncon))
					job.submit(3)
					time.sleep(0.5)
	Ncon +=1