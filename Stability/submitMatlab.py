__author__ = 'fanzheng'

from General import *
import Cluster, Stability

SB = selfbin(sys.argv[0])
par = argparse.ArgumentParser()
par.add_argument('--l', required = True, help = 'a list of mutations')
par.add_argument('--path', default= SB +'/MATLAB', help = 'the path for MATLAB scripts')
par.add_argument('--head', required = True, help = 'header of the .seq files')
par.add_argument('--ext', required = True, help = 'the extension of the output .mat file')
par.add_argument('--lamb', default = '[]', help = 'control the stength of constraint')
par.add_argument('--noprior', action = 'store_true', help = 'whether use prior')
par.add_argument('--puse', action = 'store_true', help = 'if true, only load the data and output parameter usage')
par.add_argument('--conpot', required = True, help = 'file for contact potential')
args = par.parse_args()

dirs = os.walk('.').next()[1]
dirs.sort()

dirs =[]
with open(args.l) as mutf:
    for l in mutf:
        mut = Stability.getMutation(l.strip())
        if not mut.dir in dirs:
            dirs.append(mut.dir)

noprior, puse = '[]', '[]'
if args.noprior:
    noprior = '1'
if args.puse:
    puse = '1'

for d in dirs:
    if os.path.isfile(d + '/' + d + '.' + args.ext + '.mat'):
        continue
    cmds = []
    cmds.append(' '.join(['matlab', '-nodisplay', '-nojvm', '-r', '"addpath(\'' + args.path + '\');main(\'' + d + '\',\'' + args.head + '\',\'' + d + '/' + d + '.' + args.ext + '\',' + args.lamb + ',' + noprior + ',' + puse + ',\'' + args.conpot + '\');\"']))

    Cluster.qsub(cmds, hrs = 3)