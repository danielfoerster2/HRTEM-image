#!/usr/bin/python3

import os

import numpy as np


file_xyz = os.getenv('file_xyz')
file_cel = os.getenv('file_cel')


with open(file_xyz, 'r') as xyz_file:
    lines = xyz_file.readlines()
    n_atoms = int(lines[0].strip())
    
    comment_line = lines[1].strip()
    i1 = comment_line.find('"')
    comment_line = comment_line[i1+1:]
    i2 = comment_line.find('"')
    comment_line = comment_line[:i2]
    comment_line = comment_line.split()
    a = np.linalg.norm(np.array([float(comment_line[0]), float(comment_line[1]), float(comment_line[2])]))*0.1
    b = np.linalg.norm(np.array([float(comment_line[3]), float(comment_line[4]), float(comment_line[5])]))*0.1
    c = np.linalg.norm(np.array([float(comment_line[6]), float(comment_line[7]), float(comment_line[8])]))*0.1
    
    
    
    pos = np.zeros((n_atoms, 3))
    symbols = []
    atom_lines = lines[2:2+n_atoms]

    for i_atom in range(len(atom_lines)):
        line = atom_lines[i_atom].split()
        pos[i_atom] = np.array([float(line[1]), float(line[2]), float(line[3])])*0.1
        symbols.append(line[0])

atom_types = set(symbols)

with open(file_cel, 'w') as cel_file:
    res = ""
    for elem in atom_types:
        res += elem
    print(f"{res} {n_atoms}", file=cel_file)
    print(0, a, b, c, 90.0, 90.0, 90.0, file=cel_file)
    for i_atom in range(n_atoms):
        x = pos[i_atom][0]/a
        y = pos[i_atom][1]/b
        z = pos[i_atom][2]/c
        occupancy = 1.0
        if symbols[i_atom] == 'Ag':
            debye_waller_factor = 0.008
        elif symbols[i_atom] == 'Co':
            debye_waller_factor = 0.004
        else:
            exit(f"Error: unknown element {symbols[i_atom]}")
        print(symbols[i_atom], x, y, z, occupancy, debye_waller_factor, 0, 0, 0, file=cel_file)