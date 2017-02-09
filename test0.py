#!/usr/bin/env python3
# encoding utf-8

import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

import lat
import lat.soft_hard
import lat.define_phase
import lat.read_atoms

folder = ('/media/anton/Seagate Expansion Drive/Article-MMT' +
          '/Cluster calculations for article/BiggerSystems' +
          '/Comp/5chains/0 - relaxation/4 - 1752718/')

def bond_orientation():
    orientation = [[0, 0] for i in range(400)]
    for j in range(1, 2):
        fname = folder + 'co.' + str(j * 50000) + '.data'
        (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
        for i in range(len(atoms)):
            soft_hard = lat.soft_hard.soft_hard(len(atoms), atoms[i])
            define_phase = lat.define_phase.define_phase(len(atoms), i + 1)

            if (soft_hard == 1 and define_phase != 1) or (soft_hard == 2 and define_phase == 1):
                print(soft_hard, define_phase, atoms[i])

bond_orientation()
