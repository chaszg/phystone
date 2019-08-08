#!/usr/bin/env python
import sqlite3 as lite
from dict_definitions import *
from subprocess import Popen,PIPE
from jasp import *
import os 

#sql_file = '/ihome/jkeith/kas389/p2_alchemy/alchemy-data.db'
sql_file = os.environ['SQL_PATH']

def add_refdata(energy,**args):
    '''
    Writes a reference calculation entry to the refdata table
    '''
    cov = '{}x{}'.format(*args['cov'])
    data = [args['host'],args['fac'],args['ads'],args['site'],cov,
            energy,args['root']]
    con = lite.connect(sql_file)
    cur = con.cursor()
    try:
        cur.execute("INSERT into refdata VALUES(?,?,?,?,?,?,?)",((data)))
    except lite.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)
    con.commit()
    con.close()

def energies_are_fine(slabdir, adsdir, ads, alloy):
    p = Popen('grep "energy  without entropy" {0}/{1}/OUTCAR | tail -n 1'.format(slabdir,alloy), stdin=PIPE, stdout=PIPE, stderr=PIPE,shell=True)
    sline = p.communicate()[0]
    p = Popen('grep "energy  without entropy" {0}/{1}_{2}/OUTCAR | tail -n 1'.format(adsdir,alloy,ads), stdin=PIPE, stdout=PIPE, stderr=PIPE,shell=True)
    aline = p.communicate()[0]
    slab_en, ads_en = 0, 0
    if sline:
        slab_en = sline.split()[-1]
    if aline:
        ads_en = aline.split()[-1]
    if not (slab_en == 0 or ads_en == 0):
        if len(slab_en) <= 14 and len(ads_en) <= 14:
            return True
    return False

def data_in_db(data, nullstr, calctype):
    '''
    Returns true if an entry for the alloy exists in the database
    '''
    con = lite.connect(sql_file)
    cur = con.cursor()
    if calctype == 'slab':
        if nullstr == 'is null':
            bdata = data[:2] + data[4:6] + data[-1:]
        else:
            bdata = data[:2] + data[4:]
        try:
            cur.execute('SELECT * from alloydata WHERE '
                        'host = ? AND facet = ? AND ads is null AND '
                        'bindsite is null AND cov = ? AND alloy = ? AND '
                        'root {} AND alloytype = ?'.format(nullstr),
                        (bdata))
            rows = cur.fetchall()
        except lite.Error, e:
            print "Error %s:" % e.args[0]
            sys.exit(1)

    elif calctype == 'ads':
        if nullstr == 'is null':
            bdata = data[:-2] + data[-1:]
        else:
            bdata = data
        try:
            cur.execute('SELECT * from alloydata WHERE '
                        'host = ? AND facet = ? AND ads = ? AND '
                        'bindsite = ? AND cov = ? AND alloy = ? AND '
                        'root {} AND alloytype = ?'.format(nullstr),
                        (bdata))
            rows = cur.fetchall()
        except lite.Error, e:
            print "Error %s:" % e.args[0]
            sys.exit(1)

    con.commit()
    con.close()
    if len(rows) == 0:
        return False
    else:
        return True

def write_db(data):
    '''
    Writes an alloy entry to the alloydata table
    '''
    con = lite.connect(sql_file)
    cur = con.cursor()
    try:
        cur.execute('INSERT into alloydata '
                    'VALUES(?,?,?,?,?,?,?,?,?,?,?,?)',((data)))
    except lite.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)
    con.commit()
    con.close()

def add_alloydata(**args):
    cov = args['cov']
    ads = args['ads']
    site = args['site']
    root = args['root']
    if root == None:
        nullstr = 'is null'
    else:
        nullstr = '= ?'
    alloytype = args['alloytype']
    slabdir = 'slab/{1}x{2}/'.format(args['host'],cov[0],cov[1])
    adsdir = '{2}_BE/{3}/{0}x{1}/'.format(cov[0],cov[1],ads,args['site'])
    str_index = len(ads) + 1

    for alldir in os.listdir(adsdir):
        alloy = alldir[:-str_index]
        if (not alloy.startswith(args['host'])):
            allsites = [string.split('_') for string in alloy.split('-')]
            isite = allsites[0][0]
            asite = allsites[1][0]
            hostnum = ele_dict[args['host']]
            solnum = ele_dict[allsites[0][1]]
            deltaz = abs(hostnum - solnum)
            cov = '{}x{}'.format(*args['cov'])
            efine = energies_are_fine(slabdir,adsdir,ads,alloy)

            queryargs = [args['host'],args['fac'],None,None,cov,alloy,
                         root,alloytype]
            data_exists = data_in_db(queryargs, nullstr, calctype='slab')
            if (not data_exists and efine):
                with jasp(slabdir+'/'+alloy) as calc:
                    try:
                        atoms = calc.get_atoms()
                        energy = atoms.get_potential_energy()
                        data = queryargs[:-2] + [deltaz,isite,asite,energy] + queryargs[-2:] 
                        write_db(data)
                    except:
                        print alloy, 'Something went wrong: Slab'

            queryargs = [args['host'],args['fac'],ads,site,cov,alloy,
                         root,alloytype]
            data_exists = data_in_db(queryargs, nullstr, calctype='ads')
            if efine and data_exists==False:
                with jasp(adsdir+'/'+alloy+'_'+ads) as calc:
                    try:
                        atoms = calc.get_atoms()
                        energy = atoms.get_potential_energy()
                        data = queryargs[:-2] + [deltaz,isite,asite,energy] + queryargs[-2:] 
                        write_db(data)
                    except:
                        print alloy, 'Something went wrong: Ads'
