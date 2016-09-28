__author__ = 'fanzheng'

from General import *
import Master

par = argparse.ArgumentParser()
par.add_argument('--p', required=True, help='input PDB file')
par.add_argument('--v', action='store_true', help='if true keep intermediate data')
par.add_argument('')