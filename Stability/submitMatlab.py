__author__ = 'fanzheng'

from General import *
import Stability, Cluster

SB = selfbin(sys.argv[0])
par = argparse.ArgumentParser()
par.add_argument('--path', default= SB +'/MATLAB', help = 'the path for MATLAB scripts')
par.add_argument('--head', required = True, help = 'header of the .seq files')
par.add_argument('--ext', required = True, help = 'the extension of the output .mat file')
par.add_argument('--lamb', default = '[]', help = 'control the stength of constraint')
par.add_argument('--puse', action = 'store_true', help = 'if true, only load the data and output parameter usage')
args = par.parse_args()

dirs = os.walk('.').next()[1]
dirs.sort()
with open(args.tab) as tf:
    for tl in tf:
        mut = Stability.getMutation(tl.strip())
        if not mut.dir in dirs:
            dirs.append(mut.dir)

for d in dirs:
    if os.path.isfile(d + '/' + d + '.' + args.ext + '.mat'):
        continue
    cmds = []
    if not args.puse:
        cmds.append(' '.join(['matlab', '-nodisplay', '-nojvm', '-r', '"addpath(\'' + args.path + '\');main(\'' + d + '\',\'' + args.head + '\',\'' + d + '/' + d + '.' + args.ext + '\',' + args.lamb + ', []);\"']))
    else:
        cmds.append(' '.join(['matlab', '-nodisplay', '-nojvm', '-r', '"addpath(\'' + args.path + '\');main(\'' + d + '\',\'' + args.head + '\',\'' + d + '/' + d + '.' + args.ext + '\',' + args.lamb + ', 1);\"']))
    Cluster.qsub(cmds, hrs = 3)