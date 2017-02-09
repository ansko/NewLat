#!/usr/bin/env python3

import math
import copy
import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

import lat.read_atoms

folder = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Polymer/2 - 300K relaxation/1748460/'

# формируется имя файла для чтения
def make_file_name(j):
    fname = folder + 'co.' + str(j * 50000 + 0) + '.data'
    return fname

def atomic_mobility(multiplier=10):
    displacements = [[0, 0] for i in range(5000)]

    (atoms_s, bounds_s, bonds_s, angles_s) = lat.read_atoms.read_atoms(make_file_name(1))
    (atoms_e, bounds_e, bonds_e, angles_e) = lat.read_atoms.read_atoms(make_file_name(2))

    #pprint(atoms_s)
    #return None

    for j in range(3, 51):
        (atoms_s, bounds_s, bonds_s, angles_s) = (atoms_e, bounds_e, bonds_e, angles_e)
        (atoms_e, bounds_e, bonds_e, angles_e) = lat.read_atoms.read_atoms(make_file_name(j))

        for i in range(1, len(atoms_s)):
            new = atoms_e[i]
            old = atoms_s[i]

            displacements[int(500 + multiplier * new[6])][0] += 1

            # x
            delta = new[4] - old[4]
            delta2 = bounds_s[1] - bounds_s[0] - abs(delta)
            if abs(delta) < abs(delta2):
                displacements[int(500 + multiplier * new[6])][1] += delta**2
            else:
                displacements[int(500 + multiplier * new[6])][1] += delta2**2
            # y
            delta = new[5] - old[5]
            delta2 = bounds_s[3] - bounds_s[2] - abs(delta)
            if abs(delta) < abs(delta2):
                displacements[int(500 + multiplier * new[6])][1] += delta**2
            else:
                displacements[int(500 + multiplier * new[6])][1] += delta2**2
            # z
            delta = new[6] - old[6]
            delta2 = bounds_s[5] - bounds_s[4] - abs(delta)
            if abs(delta) < abs(delta2):
                displacements[int(500 + multiplier * new[6])][1] += delta
            else:
                displacements[int(500 + multiplier * new[6])][1] += delta2**2

    for i in range(len(displacements)):
        if(displacements[i][0] != 0):
            print((i - 500) / multiplier, displacements[i][0],  displacements[i][1] / displacements[i][0] / 49)

atomic_mobility()
