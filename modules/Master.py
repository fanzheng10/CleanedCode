from General import *

def confind(mbin = PATH_confind, *args, **kwargs):
    cmd  = [mbin]
    for arg in args:
        cmd.append(arg)

    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])
    sub.call(cmd)


def createPDS(mbin = PATH_master, dry = False, **kwargs):
    cmd = [mbin + '/createPDS']

    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])

    if dry == True:
        return ' '.join(cmd)
    else:
        sub.call(cmd)
        

def masterSearch (mbin = PATH_master, dry = False, rmsdcut = 2.0, bbrmsd = True, *args, **kwargs):
    '''most arguments is the same with master program; if dry is True, have a dry run, only return the command'''
        
    cmd = [mbin + '/master']
        
    cmd.extend(['--rmsdCut', str(rmsdcut)])
    if bbrmsd:
        cmd.append('--bbRMSD')
    
    for arg in args:
        cmd.append(arg)

    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])

    if dry == True:
        return ' '.join(cmd)
    else:
        sub.call(cmd)


def matchInFile(mbin = PATH_master, dry = False, otype = 'match', bbrmsd = True, **kwargs):
    cmd = [mbin + '/master']
    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])
    if dry == True:
        return ' '.join(cmd)
    else:
        sub.call(cmd)

    
def extractPDB(mbin = PATH_master, dry = False, **kwargs):
    cmd = [mbin + '/extractPDB']
    for key, value in kwargs.iteritems():
        cmd.extend(['--'+key, str(value)])
    if dry == True:
        return ' '.join(cmd)
    else:
        sub.call(cmd)