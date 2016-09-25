from General import *
import Master, PDB, Cluster
from prody import *

SELFBIN = selfbin(sys.argv[0])

par = argparse.ArgumentParser()
par.add_argument('--p', required = True, help ='PDB file, N2P2 or M3P6')
par.add_argument('--range', required = True, nargs = '+', help = 'range of residues defining the binding site')
par.add_argument('--pname', required = True, help = 'name for this peptide')
par.add_argument('--pseq', required = True, nargs = '+', help = 'peptide sequence, a list with 3-letter names (or more for non-standard)')
par.add_argument('--m', required = True, help = 'template complex structure for generating starting conformation')
par.add_argument('--dname', required = True, help = 'directory name for this starting conformation, like model0.out')
par.add_argument('--ab', choices = [0, 1, 2], default = 0, type = int,
                 help = '0: refinement only; 1: use ab-initio, use fragment picker; 2: use ab-initio, with already made fragment')
par.add_argument('--low', action = 'store_true', help = 'if true, apply low resolution optimization')
par.add_argument('--v', action = 'store_true', help = 'if true, copy everything back, more volume, otherwise, just the score file')
par.add_argument('--capC', action = 'store_true', help = 'if true, also cap C-terminus of the peptide')
args = par.parse_args()

mdir = args.pname + '/' + args.dname
if not os.path.isdir(mdir):
    os.makedirs(mdir)
odir = os.getcwd()
os.chdir(mdir)

# prepare starting conformation
p = parsePDB(args.p, chain = 'A') # assume the receptor is chain A

# parse the range
rangeS = []
for r in args.range:
    r = r.split('-')
    (start, end) = (r[0], r[-1])
    end = str(int(end)+1)
    rangeS.append(start +':'+end)
rangeS = ' '.join(rangeS)
pocket = p.select('resnum '+ rangeS)
writePDB('_pocket.pdb', pocket.copy()) # binding pocket

# use master to create pds files
Master.createPDS(type='query', pdb='_pocket.pdb')
tem = '_'+removePath(args.m)
os.system('cp '+args.m +' '+tem)
Master.createPDS(type='target', pdb=tem)
# search pocket in the template
Master.masterSearch(rmsdcut = 3.0, bbrmsd = False, query='_pocket.pds', target=tem.replace('pdb', 'pds'), topN =1, matchOut='_pocket.match')
# generate full match file for the template
Master.matchInFile(query='_pocket.pds', matchIn='_pocket.match', structOut='_match', outType='full')

ptem = parsePDB('_match/full1.pdb') # this is the superimposed templated structure
assert ptem.numChains() == 2 
chains = []
for chain in ptem.iterChains():
    chains.append(chain)
if chains[0].numResidues() < chains[1].numResidues(): # identify the peptide chain
    pepchain = chains[0]
else:
    pepchain = chains[1]
pepchain.setChid('B') # for convenience always set peptide chain id to 'B'

# start and end residue for the peptide chain
pepstart = p.numResidues() + 1
pepend = pepstart + pepchain.numResidues() -1

# add peptide backbone from template to receptor
rec = p + pepchain.copy()
writePDB('_start.pdb', rec)

pepseq = args.pseq

# mutate the side chains to the specified peptide sequence
import Rosetta # import late because this step is slow
Rosetta.normalMut('_start.pdb', range(pepstart, pepend+1), pepseq, 'start.pdb')

# clean up the temporary files
os.system('rm -r _*')
os.chdir(odir)

# the fpd part
pepfrags = os.path.realpath(odir+'/pepfrags/')
if not os.path.isdir(pepfrags):
    os.makedirs(pepfrags)
flag = SELFBIN + '/flags'
pflag = SELFBIN + '/prepack_flags'
rosetta_db = PATH_rosetta + '/main/database'
rosetta_exe = PATH_rosetta + '/main/source/bin'

ldir = Cluster.createLocalSpace()
os.system('cp '+ mdir + '/start.pdb ' + ldir + '/start.pdb')
os.system('cp '+ mdir + '/start.pdb ' + ldir + '/native.pdb')
os.chdir(ldir)

if args.ab != 0: # if need to change the number of output model, modify the flag files
    # make fragments
    fragf = '/'.join([odir, args.pname, '.fragready'])
    if args.ab == 2:
        while os.path.isfile(fragf) == False:
            time.sleep(10)
    else:
        os.system(PATH_fpddemo + '/scripts/prep_abinitio.sh 0000')
        os.system('cp -dr frags '+pepfrags+'/'+args.pname)
        os.system('touch '+fragf)
    # copy stuff
    os.system('cp '+flag+' flags')
    os.system('cp '+pflag+' prepack_flags')
    if not os.path.isdir('frags'):
        os.system('cp -r '+pepfrags+'/'+args.pname+' frags')
    # run FPD
    os.system(PATH_fpddemo +'/prepack_example')
    os.system(PATH_fpddemo +'/run_example')
else:
    # prepack_flags = ['-database', rosetta_db, '-s start.pdb', '-native native.pdb', '-ex1', '-ex2aro', '-use_input_sc', '-unboundrot native.pdb', '-flexpep_prepack', '-nstruct 1']
    prepack_flags = ['-database', rosetta_db,
                     '-s start.pdb',
                     '-ex1',
                     '-ex2aro',
                     '-flexpep_prepack',
                     '-nstruct 1']
    prepack_flags = ' '.join(prepack_flags)
    os.system(rosetta_exe + '/FlexPepDocking.linuxgccrelease '+ prepack_flags)
    os.system('mv start_0001.pdb start.ppk.pdb; mv score.sc ppk.score.sc')
    
    if args.capC:
        Rosetta.addPepTer('start.ppk.pdb', 1, pepseq[0], 'start.pdb', plen = p.numResidues(), resnC = 6,  resnameC = pepseq[-1])
    # if only add N-ter cap
    else:
        Rosetta.addPepTer('start.ppk.pdb', 1, pepseq[0], 'start.pdb', plen = p.numResidues())

    os.system('cp start.pdb native.pdb')

    rosetta_flags = ['-database', rosetta_db,
                     '-s', 'start.pdb',
                     '-out:file:silent', 'decoys.silent',
                     '-out:file:silent_struct_type', 'binary',
                     '-pep_refine',
                     '-ex1',
                     '-ex2aro',
                     '-nstruct', '200',
                     # '-nstruct', '5', # 5 is for test, 500 is for production
                     '-scorefile', 'score.sc',
                     '-chemical:patch_selectors', 'PEPTIDE_CAP', '-in:Ntermini B',
                     '-score:weights talaris2013_no_rama_paapp_fz']# custom scoring function
    # remember to remove '-use_input_sc', otherwise cannot mutate to D-AA

    if args.capC:
        rosetta_flags.append('-in:Ctermini B')
    if args.low:
        rosetta_flags.append('-lowres_preoptimize')
    rosetta_flags = ' '.join(rosetta_flags)
    os.system(rosetta_exe + '/FlexPepDocking.linuxgccrelease '+ rosetta_flags)
        
# copy stuff back
if args.v:
    os.system('cp -r * '+odir+'/'+args.pname+'/'+args.dname)
else:
    os.system('cp score.sc '+odir+'/'+args.pname+'/'+args.dname)
os.chdir(odir)
os.system('touch ' + args.pname+'/'+args.dname+'/.finished')
Cluster.destroyLocalSpace(ldir)
