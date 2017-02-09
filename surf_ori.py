#!/usr/bin/env python3
# encoding utf-8

import math

import lat.read_atoms

# 5x20
startatom = 21
middleatom = 10
endatom = 22
systemsize = 3470
folder = ('/media/anton/Seagate Expansion Drive/Article-MMT' +
          '/Cluster calculations for article/BiggerSystems' +
          '/Comp/5chains/0 - relaxation/4 - 1752718/')

# 10x20
startatom = 21
middleatom = 10
endatom = 22
systemsize = 5380
folder = ('/media/anton/Seagate Expansion Drive/Article-MMT' +
          '/Cluster calculations for article/BiggerSystems' +
          '/Comp/10chains/0 - relaxation/4 - 1752706/')


def make_file_name(j):
    fname = folder + 'co.' + str(j * 50000 + 0) + '.data'
    return fname

def surf_ori():
    px = 0
    py = 0
    pz = 0
    for filenum in range(1, 51):
        fname = make_file_name(filenum)
        parameterx = 0
        parametery = 0
        parameterz = 0
        for j in range(9):
            for i in range(10): #12 for mod 5/10 for poly
                #----------------5chains
                # for mod
                #start = j * systemsize + 720 + 70 * i + startatom
                #middle = j * systemsize + 720 + 70 * i + middleatom
                #end = j * systemsize + 720 + 70 * i + endatom
                #for poly
                #start = j * systemsize + 1560 + 382 * i + 19
                #end = j * systemsize + 1560 + 382 * i + 364

                #----------------10chains
                start = j * systemsize + 720 + 70 * i + startatom
                end = j * systemsize + 720 + 70 * i + endatom
                #for poly
                start = j * systemsize + 1560 + 382 * i + 19
                end = j * systemsize + 1560 + 382 * i + 364

                (atoms, bounds, bonds, angles) = lat.read_atoms.read_atoms(fname)
                atomstart = atoms[start - 1]
                #atommiddle = atoms[middle - 1]
                atomend = atoms[end - 1]

                #print(atomstart, atomend)
                #continue

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
                dz2 = bounds[5] - bounds[4] - abs(dx)
                if (abs(dz) < abs(dz2)):
                    pass
                else:
                    dz = -math.copysign(dz2, dz)

                costheta = dx / math.sqrt(dx**2 + dy**2 + dz**2)
                parameterx += (3 * costheta**2 - 1) / 2
                costheta = dy / math.sqrt(dx**2 + dy**2 + dz**2)
                parametery += (3 * costheta**2 - 1) / 2
                costheta = dz / math.sqrt(dx**2 + dy**2 + dz**2)
                parameterz += (3 * costheta**2 - 1) / 2
        parameterx /= 45
        parametery /= 45
        parameterz /= 45
        px += parameterx
        py += parametery
        pz += parameterz
        print(parameterx, parametery, parameterz)
    px /= 50
    py /= 50
    pz /= 50
    print(px, py, pz)
surf_ori()
