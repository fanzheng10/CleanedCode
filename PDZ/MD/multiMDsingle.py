__author__ = 'fanzheng'

from General import *
import threading
import Cluster

sb = selfbin(sys.argv[0])

par = argparse.ArgumentParser()
par.add_argument('--p', required = True, help = 'input pdb')
# par.add_argument('--v', action = 'store_true', help = 'if true, copy everything back from local directory')
par.add_argument('--change', nargs ='+', help = 'if want to change the settings in setup command, provide strings in pairs, before and after change respectively; ')
par.add_argument('--nc', default = 4, type = int, help = 'the number of cores to use')
par.add_argument('--hrs', default = 12, type = float, help = 'time allocated to run this job, in hours')
args = par.parse_args()

ldir1 = Cluster.createLocalSpace()
ldir2 = Cluster.createLocalSpace() # data files are here
default_top = PATH_thesisData + '/PDZ/Dcmap/top_all22_prot_namd.inp'
default_par = PATH_thesisData + '/PDZ/Dcmap/par_all22_prot.inp'

namd_args = ['-p', args.p,
             '-n', '1000000',
             '--ns', '1000',
             '--ts', '2',
             '--pad', '9',
             '--Pext', '1.01325',
             '--mini', '5000',
             '--top', default_top,
             '--par', default_par,
             # '--psfopts', '"first nter; last cter|first ace; last ct3"',
             '--psfopts', '"first nter; last cter|first ace; last cter"',
             '--np', '6', '--pme',
             '--odir', ldir2,
             '--tmp', ldir1,
             '--dry']
cmd_setup = ['perl -w', sb + '/mdNAMD.pl'] + namd_args
cmd_setup = ' '.join(cmd_setup)

if args.change != None:
    assert len(args.change) % 2 == 0
    for i in range(0, len(args.change), 2):
        cmd_setup = cmd_setup.replace(args.change[i], args.change[i+1])

cmd_run = ['/home/grigoryanlab/library/bin/charmrun', '+p'+str(args.nc), '/home/grigoryanlab/library/bin/namd2', 'namdrun_run.ctl >& namdrun_run.out']
cmd_run = ' '.join(cmd_run)

odir = os.getcwd()

# run setup codes
os.system(cmd_setup)
Cluster.destroyLocalSpace(ldir1)
os.chdir(ldir2)
os.system('rm -r vmd*')

# need to have a subprocess with timeout, because wants to stop MD and allow certain time before job terminated
class Command:
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            print('Thread started')
            self.process = sub.Popen(self.cmd, shell=True)
            self.process.communicate()
            print('Thread finished')

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print('Terminating process')
            self.process.terminate()
            thread.join()
        print(self.process.returncode)

# run MD
command = Command(cmd_run)
timeout = args.hrs  * 3600 - 600 # 10 min to do the data analysis and clean up
command.run(timeout)

# clean up
os.chdir(odir)
os.system('mv '+ldir2 + '/* .')