#!/usr/bin/env python3
# encoding utf-8

import math

import lat.read_atoms

# 5x20
#startatom = 21
#middleatom = 10
#endatom = 22
#systemsize = 3470
#folder = ('')

# 10x20
#startatom = 21
#middleatom = 10
#endatom = 22
systemsize = 5380

folder = ('/home/anton/Cluster/polymer/2 - 300K relaxation/1748460/')


def make_file_name(j):
    fname = folder + 'co.' + str(j * 50000 + 0) + '.data'
    return fname

def autocor_ends():
    starts = []
    for filenum in range(1, 2):
        fname = make_file_name(filenum)
        parameter = 0
        for j in range(9):
            for i in range(10): #12 for mod 5/10 for poly
                #----------------5chains
                # for mod
                #start = j * systemsize + 720 + 70 * i + startatom
                #end = j * systemsize + 720 + 70 * i + endatom
                #for poly
                #start = j * systemsize + 1560 + 382 * i + 19
                #end = j * systemsize + 1560 + 382 * i + 364

                #----------------10chains
                # for mod
                #start = j * systemsize + 720 + 70 * i + startatom
                #end = j * systemsize + 720 + 70 * i + endatom
                #for poly
                start = j * systemsize + 1560 + 382 * i + 19
                end = j * systemsize + 1560 + 382 * i + 364

                (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
                atomstart = atoms[start - 1]
                atomend = atoms[end - 1]
                #x
                dx = atomend[4] - atomstart[4]
                dx2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dx) < abs(dx2)):
                    pass
                else:
                    dx = -math.copysign(dx2, dx)
                #y
                dy = atomend[5] - atomstart[5]
                dy2 = bounds[3] - bounds[2] - abs(dy)
                if (abs(dy) < abs(dy2)):
                    pass
                else:
                    dy = -math.copysign(dy2, dy)
                #z
                dz = atomend[6] - atomstart[6]
                dz2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dz) < abs(dz2)):
                    pass
                else:
                    dz = -math.copysign(dz2, dz)
                starts.append([dx, dy, dz])
    print(len(starts))

    
    for filenum in range(1, 50):
        fname = make_file_name(filenum)
        k = 0
        parameter = 0
        for j in range(9):
            for i in range(10):
                start = j * systemsize + 1560 + 382 * i + 19
                end = j * systemsize + 1560 + 382 * i + 364
                (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
                atomstart = atoms[start - 1]
                atomend = atoms[end - 1]
                #x
                dx = atomend[4] - atomstart[4]
                dx2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dx) < abs(dx2)):
                    pass
                else:
                    dx = -math.copysign(dx2, dx)
                #y
                dy = atomend[5] - atomstart[5]
                dy2 = bounds[3] - bounds[2] - abs(dy)
                if (abs(dy) < abs(dy2)):
                    pass
                else:
                    dy = -math.copysign(dy2, dy)
                #z
                dz = atomend[6] - atomstart[6]
                dz2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dz) < abs(dz2)):
                    pass
                else:
                    dz = -math.copysign(dz2, dz)

                dx0 = starts[k][0]
                dy0 = starts[k][1]
                dz0 = starts[k][2]
                k += 1
                costheta = ((dx0 * dx + dy0 * dy + dz0 * dz) /
                            math.sqrt(dx0 ** 2 + dy0 ** 2 + dz0 ** 2) /
                            math.sqrt(dx ** 2 + dy ** 2 + dz ** 2))
                parameter += (3 * costheta ** 2 - 1) / 2
        parameter /= 90
        print(parameter)

def autocor_ends_poly():
    starts = []
    for filenum in range(1, 2):
        fname = make_file_name(filenum)
        parameter = 0
        #for j in range(9):
        for i in range(343):
                start = 192 * i + 19
                end = 192 * i + 174
                (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
                atomstart = atoms[start - 1]
                atomend = atoms[end - 1]
                #x
                dx = atomend[4] - atomstart[4]
                dx2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dx) < abs(dx2)):
                    pass
                else:
                    dx = -math.copysign(dx2, dx)
                #y
                dy = atomend[5] - atomstart[5]
                dy2 = bounds[3] - bounds[2] - abs(dy)
                if (abs(dy) < abs(dy2)):
                    pass
                else:
                    dy = -math.copysign(dy2, dy)
                #z
                dz = atomend[6] - atomstart[6]
                dz2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dz) < abs(dz2)):
                    pass
                else:
                    dz = -math.copysign(dz2, dz)
                starts.append([dx, dy, dz])
    print(len(starts))

    
    for filenum in range(1, 50):
        fname = make_file_name(filenum)
        k = 0
        parameter = 0
        #for j in range(9):
        for i in range(343):
                start = 192 * i + 19
                end = 192 * i + 174
                (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
                atomstart = atoms[start - 1]
                atomend = atoms[end - 1]
                #x
                dx = atomend[4] - atomstart[4]
                dx2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dx) < abs(dx2)):
                    pass
                else:
                    dx = -math.copysign(dx2, dx)
                #y
                dy = atomend[5] - atomstart[5]
                dy2 = bounds[3] - bounds[2] - abs(dy)
                if (abs(dy) < abs(dy2)):
                    pass
                else:
                    dy = -math.copysign(dy2, dy)
                #z
                dz = atomend[6] - atomstart[6]
                dz2 = bounds[1] - bounds[0] - abs(dx)
                if (abs(dz) < abs(dz2)):
                    pass
                else:
                    dz = -math.copysign(dz2, dz)

                dx0 = starts[k][0]
                dy0 = starts[k][1]
                dz0 = starts[k][2]
                k += 1
                costheta = ((dx0 * dx + dy0 * dy + dz0 * dz) /
                            math.sqrt(dx0 ** 2 + dy0 ** 2 + dz0 ** 2) /
                            math.sqrt(dx ** 2 + dy ** 2 + dz ** 2))
                parameter += (3 * costheta ** 2 - 1) / 2
        parameter /= 343
        print(parameter)

autocor_ends_poly()
