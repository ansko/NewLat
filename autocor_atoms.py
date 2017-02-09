#!/usr/bin/env python3
# encoding utf-8

import math
import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

import lat.get_bonds

folder = ('/media/anton/Seagate Expansion Drive/Article-MMT' +
          '/Cluster calculations for article/BiggerSystems' +
          '/Comp/5chains/2 - relaxation 300K/')

def make_file_name(j):
    fname = folder + 'co.' + str(j * 50000 + 0) + '.data'
    return fname

def autocor():
    fname_start = make_file_name(1)
    bonds_start = lat.get_bonds.get_bonds(fname_start)
    for j in range(1, 51):
        parameter = 0
        fname_current = make_file_name(j)
        bonds_current = lat.get_bonds.get_bonds(fname_current)
        for bond in range(len(bonds_current)):
            costheta = ((bonds_current[bond][0] * bonds_start[bond][0] +
                         bonds_current[bond][1] * bonds_start[bond][1] +
                         bonds_current[bond][2] * bonds_start[bond][2]) /
                        (math.sqrt(bonds_current[bond][0]**2 +
                                   bonds_current[bond][1]**2 +
                                   bonds_current[bond][2]**2) *
                        (math.sqrt(bonds_start[bond][0]**2 +
                                   bonds_start[bond][1]**2 +
                                   bonds_start[bond][2]**2))))
            parameter += (3 * costheta**2 - 1) / 2
        parameter /= len(bonds_start)
        print(j, parameter)

autocor()
