from General import *
import threading

namd_args = ['-p', 'dummy.pdb',
			 '-n', '1000000',
			 '--ns', '1000',
			 '--ts', '2',
			 '--pad', '9',
			 '--Pext', '1.01325',
			 '--mini', '5000',
			 '--top', '/home/grigoryanlab/home/fzheng/PDZ.2014/MD/Dcmap/top_all22_prot_namd.inp',
			 '--par', '/home/grigoryanlab/home/fzheng/PDZ.2014/MD/Dcmap/par_all22_prot.inp',
			 '--psfopts', '"first nter; last cter|first ace; last ct3"',
			 '--np', '6', '--pme',
			 '--odir', '.',
			 '--tmp', '.',
			 '--dry']


def MDconstructor(**kwargs):
	namd_args_mod = namd_args
	for (key, value) in kwargs:
		ind = namd_args.index(key)
		namd_args_mod[ind+1] = value
		return namd_args_mod


# use pertMin method to generate slightly different ensembles
def pertMin(coords, N=64, e=0.01):
	coordsset = []
	for n in range(N):
		for i in range(len(coords)):
			for j in range(3):
				ran = random.uniform(-e, e)
		coords[i][j] += ran
		coordsset.append(coords)
	return coordsset


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


### interact with VMDscript
def vmd_addSelection():



def vmd_parseFunction():



def vmd_executeFunction():