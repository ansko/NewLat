#!/usr/bin/env python3
# encoding utf-8

import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

#     /m11 m12 m13\
# M = |m21 m22 m23|
#     \m31 m32 m33/
def rotate(vector):
    m11 = 0.223
    m12 = 0.324
    m13 = 0.919
    m21 = 0.904
    m22 = -0.42
    m23 = -0.071
    m31 = -0.152
    m32 = 0.847
    m33 = -0.387

    for i in range(len(vector)):
        x = vector[i][0] * m11 + vector[i][1] * m21 + vector[i][2] * m31
        y = vector[i][0] * m12 + vector[i][1] * m22 + vector[i][2] * m32
        z = vector[i][0] * m13 + vector[i][1] * m23 + vector[i][2] * m33
        vector[i][0] = x
        vector[i][1] = y
        vector[i][2] = z

    return vector
        

def clean():
    coords = []
    lines = []
    fname = 'PA6x100.cif'
    f = open(fname, 'r')
    for line in f:
        lines.append(line)
        line_s = line.split()
        coords.append([float(line_s[2]), float(line_s[3]), float(line_s[4])])
    return (coords, lines)

def get_angles():
    fname = 'PA6x100.cif'
    f = open(fname, 'r')
    for line in f:
        


get_angles()
return None

(coords, lines) = clean()
ccords = rotate(coords)

new_lines = []
for i in range(len(lines)):
    line_s = lines[i].split()
    new_lines.append([line_s[0], line_s[1], coords[i][0], coords[i][1], coords[i][2], line_s[4], line_s[5], line_s[6]])

for i in range(len(new_lines)):
    for j in range(len(new_lines[i])):
        print(new_lines[i][j], end = "\t")
    print("")
