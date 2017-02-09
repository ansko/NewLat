#!/usr/bin/env python3
# coding utf-8

import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint
import lat.read_atoms
import lat.define_phase
import lat.soft_hard

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

def phase_distr(filenum = 49, multiplier = 10):
    """builds components density profile"""
    foldernames = ['/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Comp/10chains/2.2 - More relaxation 500 (wiggle)/1799293 - wiggle2/']

    distr = [[ 0, # all 
               0, # mmt
               0, # mod
               0, # N from mod
               0] # poly
             for i in range(10000)]
    new_distr = []

    for foldername in foldernames:
        for j in range(1, 1 + filenum):
            fname = foldername + '/co.' + str((j*50 + 0)*1000) + '.data'
            (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
            for i in range(len(atoms)):
                distr[int(5000 + multiplier * atoms[i][6])][0] += masses[atoms[i][2]] * 1.66
                #distr[int(250 + multiplier * atoms[i][6])][lat.define_phase.define_phase(len(atoms), i + 1)] += 1
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

phase_distr()
