__author__ = 'fanzheng'

from General import *
# run vmd script
par = argparse.ArgumentParser()
par.add_argument('--l', required = True, help = 'a list of peptides')
par.add_argument('--Dc', nargs = '*', help = 'the index of peptides with terminus from ASP side chain')
par.add_argument('--o', required = True, help = 'output file')
args = par.parse_args()

vmdbin = '~/home/Thesis/Codes/VMDscript'

lines = open(args.l).readlines()
pepids = []
for l in lines:
    pepids.append(l.strip().split()[-1])

odir = os.getcwd()
outf = open(args.o, 'w')
for pi in range(len(pepids)):
    os.chdir(odir)
    if not os.path.isdir(pepids[pi]):
        continue
    var1=1
    if args.Dc != None and pi in args.Dc:
        var1=0
    if pepids[pi][3] == 'S': # amino acid at -2
        var2=0
    elif pepids[pi][3] == 'T':
        var2=1
    else:
        continue
    for i in range(64):
        os.chdir(odir +'/'+pepids[pi])
        os.chdir('pert' + str(i))
        if os.path.isfile('mdRmsf.dat') and os.path.isfile('mdDist.dat'):
            continue
        os.system(PATH_vmd + ' -dispdev text -e ' + vmdbin + '/analysis.vmdscript -args ' + str(var1) + ' ' + str(var2))

    # extract features
    features = []

    # process all rmsf (11 features)
    rmsf_pi = []
    for i in range(64):
        os.chdir(odir + '/'+ pepids[pi])
        os.chdir('pert' + str(i))
        rmsf_ti = []
        with open('mdRmsf.dat') as mdr:
            for l in mdr:
                values = [float(x) for x in l.strip().split(':')[1].split()]
                rmsf_ti.append(round(sum(values)/len(values), 3))
        rmsf_pi.append(rmsf_ti)
    rmsf_pi = np.matrix(rmsf_pi, dtype='float64')
    # for each feature, get max, 75% median and median (opertion on column) 11x3 = 33 features
    rmsf_max = np.round(np.max(rmsf_pi, axis=0), 3)
    rmsf_75 = np.round(np.percentile(rmsf_pi, q=75, axis=0), 3)
    rmsf_median = np.round(np.median(rmsf_pi, axis=0), 3)

    # process all distance (11x3=33 features)
    dist_pi = []
    for i in range(64):
        os.chdir(odir + '/'+ pepids[pi])
        os.chdir('pert' + str(i))
        dist_ti = []
        with open('mdDist.dat') as mdd:
            lines = mdd.readlines()
            for li in range(len(lines)):
                if li % 2 ==0:
                    continue
                line =lines[li]
                values = np.array([round(float(x), 3) for x in line.strip().split()], dtype ='float64')
                # for each feature, get mean and median, fraction and frames in which distance larger than 3.5 A (about the donor-accpetor distance )
                dist_ti_mean = np.mean(values)
                dist_ti_median = np.median(values)
                dist_ti_frac = 1.0 * len([x for x in line.strip().split() if float(x) > 3.5]) / len(line.strip().split())
                dist_ti.extend([dist_ti_mean, dist_ti_median, dist_ti_frac])
        dist_pi.append(dist_ti)
    dist_pi = np.matrix(dist_pi)
    # for each feature, get max, 75% median and median (opertion on column) 33x3 = 99 features
    dist_max = np.round(np.max(dist_pi, axis=0), 3)
    dist_75 = np.round( np.percentile(dist_pi, q=75, axis=0), 3)
    dist_median = np.round(np.median(dist_pi, axis=0), 3)

    outstr = [pepids[pi]] + map(str, rmsf_max.flatten()) + map(str, rmsf_75) + map(str, rmsf_median.flatten()) + map(str, dist_max.flatten())+ map(str, dist_75) + map(str, dist_median.flatten())
    outstr = '\t'.join(outstr)
    outf.write(outstr + '\n')
# 99 + 33 = 132 features together. Need to do PCA

outf.close()

