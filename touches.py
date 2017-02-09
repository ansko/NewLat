#!/usr/bin/env python3
# encoding utf-8

#from goto import with_goto

import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

import lat.read_atoms

foldernames = ['/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Comp/5chains/1 - Quick/2 - relaxation 300K/']


def make_file_name(j):
    fname = foldernames[0] + 'co.' + str(j * 50000 + 0) + '.data'
    return fname

def phase_distr(multiplier=10):
    """builds components density profile"""
    masses_poly = [
    0,
    1.00797,
    12.0112,
    1.00797,
    12.0112,
    15.9994,
    14.0067,
    14.0067,
    12.0112,
    ]

    masses_comp = [
    0,
    26.981540, # ao
    24.305000, # mgo
    28.085500, # st
    15.999400, # ob
    15.999400, # oh
    15.999400, # obts
    15.999400, # ohs
    1.007970,  # ho
    14.006700, # n4
    12.011150, # c2
    12.011150, # c3
    1.007970,  # h
    1.007970,  # hn
    12.011150, # c'
    15.999400, # o'
    14.006700, # n2
    14.006700, # n
    ]

    masses = masses_comp

    distr = [[ 0, # all 
               0, # mmt
               0, # mod
               0, # N from mod
               0] # poly
             for i in range(10000)]
    new_distr = []

    for foldername in foldernames:
        for j in range(1, 1 + 30):
            fname = foldername + '/co.' + str((j*50 + 0)*1000) + '.data'
            (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
            for i in range(len(atoms)):
                distr[int(5000 + multiplier * atoms[i][6])][0] += masses[atoms[i][2]] * 1.66
                distr[int(5000 + multiplier * atoms[i][6])][1] += masses[atoms[i][2]] * 1.66

    for i in range(len(distr)):
        #if distr[i][0] > 0:
            new_distr.append([(i - 5000) / multiplier,
                              multiplier * distr[i][0] / 50 / (bounds[1] - bounds[0]) / (bounds[3] - bounds[2]),
                              multiplier * distr[i][1] / 50 / (bounds[1] - bounds[0]) / (bounds[3] - bounds[2]),
                              multiplier * distr[i][2] / 50 / (bounds[1] - bounds[0]) / (bounds[3] - bounds[2]),
                              multiplier * distr[i][3] / 50 / (bounds[1] - bounds[0]) / (bounds[3] - bounds[2]),
                              multiplier * distr[i][4] / 50 / (bounds[1] - bounds[0]) / (bounds[3] - bounds[2])])

    for i in range(len(new_distr)):
        for j in range(len(new_distr[i])):
            print(new_distr[i][j], end=' ')
        print('')

def touches(lowz=-45, highz=-37):
    numberoftouches = [0 for i in range(20)]
    for filenum in range(1, 51):
        fname = make_file_name(filenum)
        (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
        for cellnum in range(9):
            for molnum in range(5): #12 for mod 5/10 for poly
                #----------------smallsystems
                # mixed, segregated
                #start = cellnum * 3480 + 1560 + 192 * molnum + 1
                #end = cellnum * 3480 + 1560 + 192 * molnum + 192
                #----------------5chains
                # for mod
                #int start = cellnum * 3470 + 720 + 70 * molnum + 1
                #int end = cellnum * 3470 + 720 + 70 * molnum + 70
                #for poly
                start = cellnum * 3470 + 1560 + 382 * molnum + 1
                end = cellnum * 3470 + 1560 + 382 * molnum + 382

                #----------------10chains
                # for mod
                #int start = cellnum * 5380 + 720 + 70 * molnum + 1
                #int end = cellnum * 5380 + 720 + 70 * molnum + 70
                #for poly
                #start = cellnum * 5380 + 1560 + 382 * molnum + 1
                #end = cellnum * 5380 + 1560 + 382 * molnum + 382

                closeatoms = []
                faratoms = []
                for atomnum in range(start, end):
                    atom = atoms[atomnum]
                    if min(abs(atom[6] - lowz), abs(atom[6] - highz)) < 5:
                        closeatoms.append(atom[0])
                    if min(abs(atom[6] - lowz), abs(atom[6] - highz)) > 10:
                        faratoms.append(atom[0])

                for m in range(len(closeatoms)):
                    for i in range(len(closeatoms)):
                        if len(closeatoms) < 2:
                            break
                        flag = 0
                        j = len(closeatoms) - 1 - m
                        if j < 0:
                            break
                        tmpmin = closeatoms[j - 1]
                        tmpmax = closeatoms[j]
                        for k in range(len(faratoms)):
                            if tmpmin < faratoms[k] < tmpmax:
                                flag = 1
                        if flag == 0:
                            closeatoms.pop(j)

                            
                numberoftouches[len(closeatoms)] += 1 
    for i in range(len(numberoftouches)):
        numberoftouches[i] /= 50
    pprint(numberoftouches)
    return None
           
#phase_distr()
touches()
