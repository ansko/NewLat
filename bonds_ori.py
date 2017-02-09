#!/usr/bin/env python3

import math
import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

import lat.read_atoms
import lat.define_phase

folder = ('/media/anton/Seagate Expansion Drive/Article-MMT' +
          '/Cluster calculations for article/BiggerSystems' +
          '/Comp/5chains/0 - relaxation/4 - 1752718/')

folder = ('/media/anton/Seagate Expansion Drive/Article-MMT' +
          '/Cluster calculations for article/BiggerSystems' +
          '/Comp/10chains/0 - relaxation/4 - 1752706/')

multiplier = 5

def bond_orientation():
    orientation = [[0, 0] for i in range(400)]
    for j in range(1, 11):
        fname = folder + 'co.' + str(j * 50000) + '.data'
        (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
        for i in range(len(bonds)):
            x1 = atoms[bonds[i][2] - 1][4]
            x2 = atoms[bonds[i][3] - 1][4]
            y1 = atoms[bonds[i][2] - 1][5]
            y2 = atoms[bonds[i][3] - 1][5]
            z1 = atoms[bonds[i][2] - 1][6]
            z2 = atoms[bonds[i][3] - 1][6]
            xlen = min(abs(x2 - x1), abs(bounds[1] - bounds[0] - abs(x2 - x1)))
            ylen = min(abs(y2 - y1), abs(bounds[3] - bounds[2] - abs(y2 - y1)))
            zlen = min(abs(z2 - z1), abs(bounds[5] - bounds[4] - abs(z2 - z1)))
            bond_len = math.sqrt(xlen**2 + ylen**2 + zlen**2)
            if bond_len > 5:
                print(bond_len)
                print("Error!")

            if abs(z2 - z1) < abs(bounds[5] - bounds[4] - abs(z2 - z1)):
                center_z = (z2+ z1) / 2
            elif bounds[5] - max(z1, z2) > min(z1, z2) - bounds[4]:
                center_z = (max(z1, z2) + abs(bounds[5] - bounds[4] - abs(z2 - z1))) / 2
            else:
                center_z = (min(z1, z2) - abs(bounds[5] - bounds[4] - abs(z2 - z1))) / 2

            if bonds[i][1] == 13:
                orientation[int(center_z * multiplier)][0] += 1
                orientation[int(center_z * multiplier)][1] += (3 * (zlen / bond_len)**2 - 1) / 2

    for i in range(len(orientation)):
        if orientation[i][0] != 0:
            print(i / multiplier, orientation[i][0], orientation[i][1] / orientation[i][0])

def angle_orientation():
    orientation = [[0, 0] for i in range(400)]
    for j in range(1, 31):
        fname = folder + 'co.' + str(j * 50000) + '.data'
        (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
        pprint(angles)
        for i in range(len(bonds)):
            x1 = atoms[bonds[i][2] - 1][4]
            x2 = atoms[bonds[i][3] - 1][4]
            y1 = atoms[bonds[i][2] - 1][5]
            y2 = atoms[bonds[i][3] - 1][5]
            z1 = atoms[bonds[i][2] - 1][6]
            z2 = atoms[bonds[i][3] - 1][6]
            xlen = min(abs(x2 - x1), abs(bounds[1] - bounds[0] - abs(x2 - x1)))
            ylen = min(abs(y2 - y1), abs(bounds[3] - bounds[2] - abs(y2 - y1)))
            zlen = min(abs(z2 - z1), abs(bounds[5] - bounds[4] - abs(z2 - z1)))
            bond_len = math.sqrt(xlen**2 + ylen**2 + zlen**2)
            if bond_len > 5:
                print(bond_len)
                print("Error!")

            if abs(z2 - z1) < abs(bounds[5] - bounds[4] - abs(z2 - z1)):
                center_z = (z2+ z1) / 2
            elif bounds[5] - max(z1, z2) > min(z1, z2) - bounds[4]:
                center_z = (max(z1, z2) + abs(bounds[5] - bounds[4] - abs(z2 - z1))) / 2
            else:
                center_z = (min(z1, z2) - abs(bounds[5] - bounds[4] - abs(z2 - z1))) / 2

            if bonds[i][1] == 5:
                orientation[int(center_z * multiplier)][0] += 1
                orientation[int(center_z * multiplier)][1] += (3 * (zlen / bond_len)**2 - 1) / 2

    for i in range(len(orientation)):
        if orientation[i][0] != 0:
            print(i, orientation[i][0], orientation[i][1])

bond_orientation()
#angle_orientation()
