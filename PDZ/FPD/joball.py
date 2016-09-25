from General import *
import PDB, Cluster
SELFBIN = selfbin(sys.argv[0])

# some constant to declare
bp_range = {'N2P2':'16-23 65-73', 'M3P6': '9-16 60-68', 'P95P3': '23-30 72-80'}
pdbbin = PATH_thesisData + '/PDZ/templatePDBs/'
templateList = PATH_thesisData + '/PDZ/' + 'domains_withPep_sim.txt'

par = argparse.ArgumentParser(description = 'execute Rosetta FlexPepDock for a list of peptide (Do you remember pyr?)')
par.add_argument('--p', choices = bp_range.keys(), help ='receptor structure')
par.add_argument('--l', required = True, help = 'a list of all peptide sequences')
par.add_argument('--n', help = 'the first and the last line to execute in the list file, like 1-10, if not defined, run all')
par.add_argument('--v', action = 'store_true', help = 'if true, return everything, otherwise just score file')
par.add_argument('--pepl', default = 6, type = int, help = 'length of peptide')
par.add_argument('--abinitio', action = 'store_true', help = 'it automatically include low resolution optimization')
par.add_argument('--low', action = 'store_true', help = 'if true, apply low resolution optimization')
par.add_argument('--capC', action = 'store_true', help = 'if true, also cap C-terminus of the peptide')
args = par.parse_args()

tems = open(templateList).readlines()
ntem = len(tems)
lines = open(args.l).readlines()
if args.n != None:
    (start, end) = map(int, args.n.split('-'))
else:
    (start, end) = (1, len(lines))
jobs = []

for i in range(len(lines)):
    if i+1 < start:
        continue
    if i+1 > end:
        break
    l = lines[i].strip()
    # parse peptide sequence, a little tricky, notice the format of peptide file
    arr = l.split()
    pepid = arr[-1]
    pep = arr[:-2]

    first = True
    for m in range(ntem):
        if tems[m][0] == '#':
            continue
        tem = tems[m].split()
        (temp, peplen) = (tem[0], tem[5])
        tempdb = pdbbin + '/' + temp + '.pdb'
        if int(peplen) < args.pepl: # peptide length in template is not long enough
            continue
        if not os.path.isfile(pepid +'/model'+str(m)+'.out/score.sc'):
            cmd = ['python', SELFBIN + '/fpd.py', 
            '--p', absPath(args.p + '.pdb'),
            '--range', bp_range[args.p], 
            '--pname', pepid, 
            '--pseq', ' '.join(pep), 
            '--m', tempdb, 
            '--dname', 'model'+str(m)+'.out']

            if not args.abinitio:
                cmd.extend(['--ab', '0']) # not use ab_initio
            elif first == True:
                cmd.extend(['--ab', '1']) # use ab_initio and create new fragment database
                first = False
            else:
                cmd.extend(['--ab', '2']) # use ab_initio but use availabel fragment database

            if args.v:
                cmd.append('--v')
            if args.low:
                cmd.append('--low')
            if args.capC:
                cmd.append('--capC')  

            cmd = ' '.join(cmd)

            job = Cluster.jobOnCluster(['pyr', cmd], pepid + ' model'+str(m), pepid +'/model'+str(m)+'.out/.finished')
            jobs.append(job)
            job.submit(12)
            time.sleep(1)

Cluster.waitJobs(jobs)

        
    
    
    