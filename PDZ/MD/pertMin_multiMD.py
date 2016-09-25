__author__ = 'fanzheng'

from General import *
from prody import *
import Cluster

par = argparse.ArgumentParser()
par.add_argument('--p', required = True, help = 'provide a template file')
par.add_argument('--n', type = int, default = 64, help = 'number of perturbed structures generated')
par.add_argument('--e', type = float, default = 0.01, help = 'the maximum perturbation for each heavy atom, in angstrom')
par.add_argument('--change', nargs = '+', help = 'see instruction of multiMDsingle.py')
par.add_argument('--opts', nargs = 2, help = 'see instruction of multiMDsingle.py')
args = par.parse_args()
SB = selfbin(sys.argv[0])
odir = os.getcwd()
hrs = 12

### pert
mol = parsePDB(args.p)
mol2 = mol.copy()
chB = mol.select('chain B')

for n in range(args.n):
    os.chdir(odir)
    dirname = 'pert' + str(n)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    os.chdir(dirname)
    # need to clean the current directory
    if len(os.listdir('.')) != 0:
        os.system('rm -r *')
    chB = mol.select('chain B')
    coords = chB.getCoords()
    # chain B are subjected to pertubation
    for i in range(len(coords)):
        # generate random number between -0.001 and 0.001
        for j in range(3):
            ran = random.uniform(-args.e, args.e)
            coords[i][j] += ran
    chB2 = mol2.select('chain B')
    chB2.setCoords(coords)
    npdb = args.p + '.pert' + str(n)
    writePDB(npdb, mol2)

    cmd = ['python', SB + '/multiMDsingle.py', '--p', npdb, '--hrs', str(hrs)]
    if args.change != None:
        cmd.append('--change')
        cmd.extend(args.change)
    cmd = ' '.join(cmd)
    Cluster.qsub([cmd], hrs = hrs, opts = ['#$ -pe smp 4', '#$ -q medium'])