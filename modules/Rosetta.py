from General import *
from rosetta import *

# refer to PyRosetta tutorial to understand the function called

CMD_cluster = '/home/ironfs/grigoryanlab/home/fzheng/PDZ.2014/FlexPepDock_AbInitio/scripts/clustering/cluster.sh\
		 500 1 score.sc native.pdb decoys.silent total_score'

def flags_Daa():
	pyrflags = ['-in:file:extra_res_fa']
	pparams = '/home/anthill/fzheng/home/software/PyRosetta.ScientificLinux-r56316.64Bit/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/d-caa/'
	params = os.walk(pparams).next()[2]
	for param in params:
		pyrflags.append(pparams + param)
	pyrflags = ' '.join(pyrflags)
	return pyrflags

def packSingleSite(pose, resns):
	scorefxn = get_fa_scorefxn()
	task = TaskFactory.create_packer_task(pose)
	task.temporarily_fix_everything()
	task.restrict_to_repacking()
	task.temporarily_set_pack_residue(resns, True)
	mover = PackRotamersMover(scorefxn, task)
	mover.apply(pose)


def normalMut(pdbf, resns, resnames, outf, repack=True, daa=False):
	assert len(resns) == len(resnames)
	if daa:
		init(options = flags_Daa())
	else:
		init()
	# init()
	pose = pose_from_pdb(pdbf)
	for i in range(len(resns)):
		if i == 0:
			mutator = MutateResidue(resns[i], resnames[i]+'_p:NtermProteinFull')
		elif i == len(resns)-1:
			mutator = MutateResidue(resns[i], resnames[i]+'_p:CtermProteinFull')
		else:
			mutator = MutateResidue(resns[i], resnames[i]) # how terminal is treated here?
		mutator.apply(pose)

		if repack == True:
			packSingleSite(pose, resns[i])
		elif isinstance(repack, list):
			if repack[i] == 1:
				packSingleSite(pose, resns[i])
	dump_pdb(pose, outf)

def addPepTer(pdbf, resnN, resnameN, outf, plen = 90, resnC = None, resnameC = None, repack = [0,0]):
	if resnC != None:
		assert resnameC != None
	pyrflags = '-in:Ntermini B -in:Ctermini B'
	# processing, remove hydrogens
	lines = open(pdbf).readlines()
	tmpf = 'tmp.pdb'
	tmpo = open(tmpf, 'w')
	for l in lines:
		if l[0:4] == 'ATOM':
			(cid, atomname, seqnum) = (l[21], l[12:16], l[22:26])
			if cid == 'B':
				if seqnum.strip() == str(resnN):
					if atomname.strip() == '1H':
						l = l.replace('1H', ' H')
					elif atomname.strip() in ['2H', '3H']:
						continue
				if resnC != None:
					if (seqnum.strip() == str(resnC)) and (atomname.strip() == 'OXT'):
						continue
			tmpo.write(l)
	tmpo.close()
	init(options=pyrflags)
	pose = pose_from_pdb(tmpf)
	mutator = MutateResidue(resnN + plen, resnameN + '_p:N_acetylated') # not sure if work if N-ter is not a standard aa
	mutator.apply(pose)
	if repack[0] == 1:
		packSingleSite(pose, resnN + plen)
	if resnC != None:
		mutator = MutateResidue(resnC + plen, resnameC + '_p:C_methylamidated')
		mutator.apply(pose)
		if repack[1] == 1:
			packSingleSite(pose, resnN + plen)
	os.remove('tmp.pdb')
	dump_pdb(pose, outf)

