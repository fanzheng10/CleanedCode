import General, Terms, Master, Stability, Cluster
import os, argparse

par = argparse.ArgumentParser()
# mutation list
par.add_argument('--l', required = True, help = 'a list of mutations')
par.add_argument('--i', default = os.getcwd(), required = True, help = 'the directory to find the input structures')
par.add_argument('--sd', help = 'the directory of the searching database')
par.add_argument('--cp', default = '', help = 'the list of contact potential')
par.add_argument('--he', help = 'the header of the output files')
par.add_argument('--pb', default = General._blast, help = 'path to BLAST')
par.add_argument('--pu', default = General._usearch, help = 'path to USearch')
par.add_argument('--c1', default = 20000, type = int, help = 'cutoff for top N matches')
par.add_argument('--c2', default = 100, type = int, help = 'the criteria of increasing the order of sub-TERMs')
par.add_argument('--c3', default = 3, type = int, help = 'the maximum order of sub-TERMs to consider')
args = par.parse_args()


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
	pdbf = args.i + General.changeExt(pdb, 'pdb')
	assert os.path.isfile(pdbf), 'the pdb file '+pdbf +' does not exist; quit...'
	confind_out = General.changeExt(pdbf, '.conf')
	if os.path.isfile(confind_out):
		continue
	else: # run confind
		Master.confind('--pp',
					   p=pdbf, o=confind_out,
					   rLib=General._rotLib)


# contact identification
pos2cons = {}
pos2pdb = {}

for pos in positions:
	odir = os.getcwd()
	os.mkdir(pos)
	os.chdir(pos)
	pid, ipos = pos.split('_')
	icid, iresnum = ipos[0], ''.join(ipos[1:])
	conListf = pos + '.conlist'
	cons, _ = Terms.contactList(profile=args.i+'/'+pid+'.conf',
								resnum=iresnum, cid=icid,
								outFile=conListf,
								dmin=0.02, monomer=True)
	pos2pdb[pos] = args.i+'/'+pid+'.pdb'
	pos2cons[pos] = cons


# recursive fragmentation and master search
for pos in positions:
	iTerms = []
	cons = pos2cons[pos]
	# create self term without contact
	aterm = Terms.Term(parent=pos2pdb[pos],
					   seed=, contact=)




# create homologs list with blast 

# homologs removal 

# redundancy removal

# 3 letter to digital


