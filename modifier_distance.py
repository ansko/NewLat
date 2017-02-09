#!/usr/bin/env python3
# coding utf-8

import copy
import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint
import lat.read_atoms
import lat.define_phase

step = 1
stepsnumber = 50
filenum = 50
multiplier = 1

def modifier_distance():
    foldernames = ['/media/anton/Seagate Expansion Drive/Article-MMT/' +
                   'Cluster calculations for article/BiggerSystems/' +
                   'Comp/5chains/0 - relaxation/1 - 1742996', 

                   '/media/anton/Seagate Expansion Drive/Article-MMT/' +
                   'Cluster calculations for article/BiggerSystems/' +
                   'Comp/5chains/0 - relaxation/2 - 1748451',

                   '/media/anton/Seagate Expansion Drive/Article-MMT/' +
                   'Cluster calculations for article/BiggerSystems/' +
                   'Comp/5chains/0 - relaxation/3 - 1751298',

                   '/media/anton/Seagate Expansion Drive/Article-MMT' +
                   '/Cluster calculations for article/BiggerSystems/' +
                   'Comp/5chains/0 - relaxation/4 - 1752718/']

    foldernames = ['/media/anton/Seagate Expansion Drive/Article-MMT/' +
                   'Cluster calculations for article/BiggerSystems/' +
                   'Comp/10chains/0 - relaxation/1 - 1742989', 

                   '/media/anton/Seagate Expansion Drive/Article-MMT/' +
                   'Cluster calculations for article/BiggerSystems/' +
                   'Comp/10chains/0 - relaxation/2 - 1748450',

                   '/media/anton/Seagate Expansion Drive/Article-MMT/' +
                   'Cluster calculations for article/BiggerSystems/' +
                   'Comp/10chains/0 - relaxation/3 - 1751266',

                   '/media/anton/Seagate Expansion Drive/Article-MMT' +
                   '/Cluster calculations for article/BiggerSystems/' +
                   'Comp/10chains/0 - relaxation/4 - 1752706/']

    for k in range(len(foldernames)):
        foldername = foldernames[k]
        if k == 0:
            filenum = 49
            delta = 1
        elif k == 1:
            filenum = 51
            delta = 0
        else:
            filenum = 50
            delta = 0
        for j in range(1, 1 + filenum):
            distr = [0 for i in range(1000)]
            fname = foldername + '/co.' + str((j*50 + delta)*1000) + '.data'
            (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
            for i in range(len(atoms)):
                if lat.define_phase.define_phase(len(atoms), i) == 2:
                    distr[50 + int(multiplier * atoms[i][6])] += 1

            tmp = []
            for i in range(len(distr)):
                if distr[i] != 0:
                    tmp.append(distr[i])
            distr = tmp

            distr = distr[20 * multiplier:-1:1] + distr[0:20 * multiplier:1]

            #pprint(distr)
            #return None

            distribution=copy.deepcopy(distr)
            estimation = 0

            distr_rev=copy.deepcopy(distr)
            distr_rev.reverse()

            for i in range(int(len(distr)/2) + 1):
                distr[i] += distr_rev[i]

            #pprint(distr)
            #return None

            for i in range(int(len(distr)/2) + 1):
                estimation += distr[i] * i

            print(estimation/840/9)


modifier_distance()
