__author__ = 'fanzheng'

# run this in correct directory
from General import *
# import PDB,  Stability
from PDB import aaindex, t2s
from Analyze import readMultiColumn
from Stability import loadParams
from pymol import cmd

cmd.util.cbaw() # color everything to white
cmd.set('label_font_id', 7)
cmd.set('label_size', 20)
cmd.set('label_position',(2,2,2))

def showRes(chain, resi, color, mode):
    sele = 'c. ' + chain + ' and i. ' + str(resi)
    cmd.show(mode, sele)
    cmd.color(color, sele)
    cmd.util.cnc(sele)
cmd.extend( "showRes", showRes )

def showSeed(mode='sticks'):
    seedname = removePath(os.getcwd()).split('_')[-1]
    showRes(seedname[0], seedname[1:], 'yellow', mode)
cmd.extend( "showSeed", showSeed )


def showCons(mode='sticks', showcond = 0):
    conlistf = removePath(os.getcwd())+'.conlist'
    (cons, conds) = readMultiColumn(conlistf, coln=[2,3])
    for i in range(len(conds)):
        if float(conds[i]) >= 0.02:
            showRes(cons[i][0], cons[i][1:], 'cyan', mode)
            if showcond: # add a label of contact degree
                cmd.label(selection ='c. ' + cons[i][0] + ' and i. ' + cons[i][1:] + ' and n. CA', expression=format(float(conds[i]), '.3f'))
cmd.extend( "showCons", showCons )


def showFrags():
    pdbs = [x for x in os.listdir('.') if x.endswith('.pdb')]
    cmd.hide('everything')
    for p in pdbs:
        cmd.load(p)
        cmd.show('cartoon', getBase(p))
        cmd.show('lines', getBase(p))
cmd.extend( "showFrags", showFrags )


def showParamsValue(mut, matf):
    conlistf = removePath(os.getcwd())+'.conlist'
    (cons, conds, conres) = readMultiColumn(conlistf, [2, 3, -1])
    optParams = loadParams(matf)
    totalScore = 0

    # center needs to be read from PDB, since conlist may not have contacts
    seedres = []
    seedname = removePath(os.getcwd()).split('_')[-1]
    seedsele = 'c. ' + seedname[0] + ' and i. ' + seedname[1:] + ' and n. CA'
    atoms = cmd.get_model(seedsele)
    for at in atoms.atom:
        seedres.append(at.resn)
    wt = t2s(seedres[0])

    selfscore = optParams[aaindex[mut]] - optParams[aaindex[wt]]
    totalScore += selfscore
    cmd.label(selection=seedsele, expression=float(selfscore))
    cons2, conres2 = [], []
    for i in range(len(cons)):
        if conds[i] > 0.02:
            cons2.append(cons[i])
            conres2.append(conres[i])
    for i in range(len(cons2)):
        paramsi = optParams[20 + 400*i: 20+400*(i+1)].reshape(20,20)
        conind = aaindex[t2s(conres2[i])]
        consele = 'c. ' + cons2[i][0] + ' and i. ' + cons2[i][1:] + ' and n. CA'
        pair_diff_i = paramsi[conind, aaindex[mut]] - paramsi[conind, aaindex[wt]]
        totalScore += pair_diff_i
        cmd.label(selection=consele, expression= pair_diff_i)

    return seedsele, totalScore
cmd.extend( "showParamsValue", showParamsValue )


def showPrediction(mut, matf):
    seedsele, totalScore = showParamsValue(mut, matf)
    cmd.label()
    cmd.label(selection=seedsele, expression=totalScore)
cmd.extend( "showPrediction", showPrediction )


def showParamsUsage(mut, matf):
    paruse = loadParams(matf, var='paramsUsage')
    conlistf = removePath(os.getcwd())+'.conlist'
    cons, conds, conres = readMultiColumn(conlistf, [2, 3, -1])

    # center needs to be read from PDB, since conlist may not have contacts
    seedres = []
    seedname = removePath(os.getcwd()).split('_')[-1]
    seedsele = 'c. ' + seedname[0] + ' and i. ' + seedname[1:] + ' and n. CA'
    atoms = cmd.get_model(seedsele)
    for at in atoms.atom:
        seedres.append(at.resn)
    wt = t2s(seedres[0])

    print(paruse[1, 1:20])
    selfuse = [paruse[aaindex[wt]], paruse[aaindex[mut]]]
    selflabel = str(selfuse[0]) + ' ' + str(selfuse[1])
    cmd.label(selection=seedsele, expression=selflabel)

    cons2, conres2 = [], []
    for i in range(len(cons)):
        if conds[i] > 0.02:
            cons2.append(cons[i])
            conres2.append(conres[i])
    for i in range(len(cons2)):
        paramsi = paruse[20 + 400*i: 20+400*(i+1)].reshape(20,20)
        conind = aaindex[t2s(conres2[i])]
        consele = 'c. ' + cons2[i][0] + ' and i. ' + cons2[i][1:] + ' and n. CA'
        pairuse = [paramsi[conind, aaindex[wt]], paramsi[conind, aaindex[mut]]]
        pairlabel = str(pairuse[0]) + ' ' + str(pairuse[1])
        cmd.label(selection=consele, expression=pairlabel)
cmd.extend( "showParamsUsage", showParamsUsage )

