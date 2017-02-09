#!/usr/bin/env python3

import math
import copy
import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

#folder = '/home/anton/Desktop/2016-11-16/pa_mob_data2/' # 160 molecules
#folder = '/home/anton/Desktop/2016-10-12/2016-11-16/very_big_poly/' # 1280 molecules
#folder = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Polymer/2 - 300K relaxation/1748460/' # 343 molecules
#folder = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Polymer/2 - 300K relaxation/1748460/'


#folder = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Comp/5chains/0 - relaxation/4 - 1752718/'
folder = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Comp/10chains/0 - relaxation/4 - 1752706/'

#folder = 'comp_mob_data2/'
masses = [ 
    0, # there is no 0 type
    1.00797,
    12.0112,
    1.00797,
    12.0112,
    15.9994,
    14.0067,
    14.0067,
    12.0112
]
masses2 = [
0,
26.9815,
24.305,
28.0855,
15.9994,
15.9994,
15.9994,
15.9994,
1.00797,
14.0067,
12.0112,
12.0112,
1.00797,
1.00797,
12.0112,
15.9994,
14.0067,
14.0067,
]
masses = masses

# прочитать датафайл и вернуть атомы и границы в списке
def read_datafile(filename):
    f = open(filename, 'r')
    flag = 0
    atoms = []
    bonds = []
    bounds = []
    atom = 0
    bond = 0
    for line in f:
        if line.endswith('xlo xhi\n'):
            line_s = line.split()
            bounds.append(float(line_s[0]))
            bounds.append(float(line_s[1]))
        if line.endswith('ylo yhi\n'):
            line_s = line.split()
            bounds.append(float(line_s[0]))
            bounds.append(float(line_s[1]))
        if line.endswith('zlo zhi\n'):
            line_s = line.split()
            bounds.append(float(line_s[0]))
            bounds.append(float(line_s[1]))
        if line.startswith('Atoms # full'):
            flag = 1
            continue
        if line.startswith('Velocities'):
            flag = 2
        if line.startswith('Bonds'):
            flag = 3
            continue
        if line.startswith('Angles'):
            break
        if flag == 1:
            line = line.split()
            atoms.append([])
            for word in line:
                try:
                    word = int(word)
                except:
                    word = float(word)
                atoms[atom].append(word)
            atom += 1
        if flag == 3:
            line = line.split()
            bonds.append([])
            for word in line:
                try:
                    word = int(word)
                except:
                    word = float(word)
                bonds[bond].append(word)
            bond += 1
    atoms.pop(0)
    atoms.pop(len(atoms) - 1)
    atoms.sort()
    bonds.pop(0)
    bonds.pop(len(bonds) - 1)
    bonds.sort()

    return (atoms, bonds, bounds)

# рассчитать смещение с учётом возможного перехода границ задаются старые
# и новые координаты чего-то одного (атома или цм); границы ячейки
def calculate_disp(old, new, bounds=None):
    displacement=[0, 0, 0]
    bounds_def = [
        3.6574177352488846e-01, 9.3767139801735425e+01,
       -1.7173796716130170e+00, 7.8533787176849529e+01,
        4.2580775568494040e+00, 4.6615800720898385e+01
    ]
    if bounds is None:
        bounds = bounds_def
    # x
    delta = new[0] - old[0]
    delta2 = bounds[1] - bounds[0] - abs(delta)
    if abs(delta) < abs(delta2):
        displacement[0] += delta
    else:
        displacement[0] = -math.copysign(delta2, delta)
    # y
    delta = new[1] - old[1]
    delta2 = bounds[3] - bounds[2] - abs(delta)
    if abs(delta) < abs(delta2):
        displacement[1] += delta
    else:
        displacement[1] = -math.copysign(delta2, delta)
    # z
    delta = new[2] - old[2]
    delta2 = bounds[5] - bounds[4] - abs(delta)
    if abs(delta) < abs(delta2):
        displacement[2] += delta
    else:
        displacement[2] = -math.copysign(delta2, delta)

    return (displacement[0], displacement[1], displacement[2])

# формируется имя файла для чтения
def make_file_name(j):
    fname = folder + 'co.' + str(j * 50000 + 0) + '.data'
    return fname

#----------------------------------------------------------------------------
def main2(mol_number=343, mol_len=192):
#1280, 160 - poly
    coms = [[0, 0, 0] for i in range(mol_number)]
    dr2 = [[0, 0, 0] for i in range(mol_number)]
    dr2ave = [0, 0, 0]

    fname = make_file_name(1)
    (atoms0, bonds0, bounds0) = read_datafile(fname)
    com = [atoms0[i][4], atoms0[i][5], atoms0[i][6]]
    coms.append(com)

    for filenum in range(1, 51):
        fname = make_file_name(filenum)
        (atoms, bonds, bounds) = read_datafile(fname)
        for i in range(len(atoms)):
            print(atoms[i])
            return None
            com = [atoms[i][4], atoms[i][5], atoms[i][6]]
            coms.append(com)

            #x
            delta = coms[i][0] - com[0]
            delta2 = bounds[1] - bounds[0] - abs(delta)
            if abs(delta) < abs(delta2):
                dr2[i][0] += delta
            else:
                dr2[i][0] -= math.copysign(delta2, delta)
            #y
            delta = coms[i][1] - com[1]
            delta2 = bounds[3] - bounds[2] - abs(delta)
            if abs(delta) < abs(delta2):
                dr2[i][1] += delta
            else:
                dr2[i][1] -= math.copysign(delta2, delta)
            #z
            delta = coms[i][2] - com[2]
            delta2 = bounds[1] - bounds[0] - abs(delta)
            if abs(delta) < abs(delta2):
                dr2[i][2] += delta
            else:
                dr2[i][2] -= math.copysign(delta2, delta)
            ##
            coms[i] = com
            dr2ave[0] += dr2[i][0]
            dr2ave[1] += dr2[i][1]
            dr2ave[2] += dr2[i][2]
        print((dr2ave[0]**2/mol_number + 
               dr2ave[1]**2/mol_number +
               dr2ave[2]**2/mol_number))
        dr2ave = [0, 0, 0]
    return None

multiplier = 5

def main3():
    disps = [[0, 0] for i in range(200 * multiplier)]
    coms = [[]]
    old = []
    new = []
    dr2 = [0, 0, 0]

    fname = make_file_name(1)
    (atoms, bonds, bounds) = read_datafile(fname)
    for i in range(len(atoms)):
        com = [atoms[i][4], atoms[i][5], atoms[i][6]]
        old.append(com)
        #coms[0].append(com)
    #return None
    for filenum in range(1, 51):
        new = []
        disps = [[0, 0] for i in range(200 * multiplier)]
        fname = make_file_name(filenum)
        (atoms, bonds, bounds) = read_datafile(fname)
        #coms.append([])
        for i in range(len(atoms)):

            com = [atoms[i][4], atoms[i][5], atoms[i][6]]
            #coms[filenum].append(com)
            new.append(com)
            #print(new)
            #return None

            #x
            #delta = coms[filenum][i][0] - coms[filenum-1][i][0]
            delta = new[i][0] - old[i][0]
            delta2 = bounds[1] - bounds[0] - abs(delta)
            if abs(delta) < abs(delta2):
                #dr2[i][0] += delta
                dr2[0] += delta**2
            else:
                #dr2[i][0] -= math.copysign(delta2, delta)
                dr2[0] += delta2**2
            #y
            #delta = coms[i][1] - com[1]
            delta = new[i][1] - old[i][1]
            delta2 = bounds[3] - bounds[2] - abs(delta)
            if abs(delta) < abs(delta2):
                #dr2[i][1] += delta
                dr2[0] += delta**2
            else:
                #dr2[i][1] -= math.copysign(delta2, delta)
                dr2[0] += delta2**2
            #z
            #delta = coms[i][2] - com[2]
            delta = new[i][2] - old[i][2]
            delta2 = bounds[1] - bounds[0] - abs(delta)
            if abs(delta) < abs(delta2):
                #dr2[i][2] += delta
                dr2[0] += delta**2
            else:
                #dr2[i][2] -= math.copysign(delta2, delta)
                dr2[0] += delta2**2
            ##
            #coms[i] = com
            #dr2ave[0] += dr2[i][0]
            #dr2ave[1] += dr2[i][1]
            #dr2ave[2] += dr2[i][2]

            disps[int(new[i][2] * multiplier)][0] += 1
            disps[int(new[i][2] * multiplier)][1] += dr2[0] + dr2[1] + dr2[2]
        old = new
        print(filenum, "--------------------------------")
        #print((dr2ave[0]**2/mol_number + 
        #       dr2ave[1]**2/mol_number +
        #       dr2ave[2]**2/mol_number))
        #dr2ave = [0, 0, 0]
        for i in range(len(disps)):
            if disps[i][0] != 0:
                #print(disps[i][0], disps[i][1] / disps[i][0])
                print(disps[i][0] / filenum, disps[i][1] / disps[i][0] / filenum)
    return None

main3()
