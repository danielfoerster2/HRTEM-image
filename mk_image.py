#!/usr/bin/python3

import os
import subprocess
import sys
from matplotlib.pyplot import plt
import numpy as np

box = 100   # EN ANGSTROM
n_variants = 10
dbf_ag = 0.008
dbf_co = 0.004
nx = 64
ny = 64
nz = int(os.getenv('nz'))
electron_energy = float(os.getenv('electron_energy'))
count_min = 500
count_max = 1000   

aberrations = np.array([[0.0, 0.0],    # image shift
              [-4.0, -8.0],  # defocus
              [0.0, 0.0],    # astigmatism
              [0.0, 0.0],    # coma
              [0.0, 0.0],    # three-lobe aberration
              [0.0, 0.0],    # spherical aberration
              [0.0, 0.0]])   # star aberration

with open(sys.argv[1], 'r') as xyz_file:
    lines = xyz_file.readlines()
    n_atoms = int(lines[0].strip())
    pos_ini = np.zeros((n_atoms, 3))
    symbols = []
    atom_lines = lines[2:2+n_atoms]

    for i_atom in range(len(atom_lines)):
        line = atom_lines[i_atom].split()
        pos_ini[i_atom] = np.array([float(line[1]), float(line[2]), float(line[3])])
        symbols.append(line[0])

    gyration_radius = np.sqrt(np.sum(np.linalg.norm(pos_ini, axis=1)**2) / n_atoms)

    radius_limit = .8*gyration_radius

    pos_atm1 = pos_ini[np.array(symbols) == 'Ag']
    pos_atm2 = pos_ini[np.array(symbols) == 'Co']

    mask_out_atm1 = np.linalg.norm(pos_atm1, axis=1) > radius_limit
    mask_out_atm2 = np.linalg.norm(pos_atm2, axis=1) > radius_limit
    
    nat1 = pos_atm1.shape[0]
    nat2 = pos_atm2.shape[0]

    nat1_out = np.sum(mask_out_atm1)
    nat2_out = np.sum(mask_out_atm2)
    nat1_in = np.sum(~mask_out_atm1)
    nat2_in = np.sum(~mask_out_atm2)

        
    if nat1 != 0:
        r_cm1 = np.sum(pos_atm1, axis=0) / nat1
    else:
        r_cm1 = np.array([np.nan, np.nan, np.nan])
    if nat2 != 0:
        r_cm2 = np.sum(pos_atm2, axis=0) / nat2
    else:
        r_cm2 = np.array([np.nan, np.nan, np.nan])

    d_com = np.linalg.norm(r_cm1 - r_cm2)





for i_variant in range(n_variants):
    phi = np.random.uniform(0, 2*np.pi)
    theta = np.random.uniform(0, np.pi)
    rot_matrix = np.array([[np.cos(phi)*np.cos(theta), -np.sin(phi), np.cos(phi)*np.sin(theta)],
                        [np.sin(phi)*np.cos(theta), np.cos(phi), np.sin(phi)*np.sin(theta)],
                        [-np.sin(theta), 0, np.cos(theta)]])
    pos = np.dot(pos_ini, rot_matrix) * np.random.uniform(0.95, 1.05)

    with open("coord.cel", 'w') as cel_file:
        print("", file=cel_file)
        print(0, box, box, box, 90.0, 90.0, 90.0, file=cel_file)
        for i_atom in range(n_atoms):
            x = pos[i_atom][0]/box
            y = pos[i_atom][1]/box
            z = pos[i_atom][2]/box
            occupancy = 1.0
            if symbols[i_atom] == 'Ag':
                debye_waller_factor = dbf_ag # XXX
            elif symbols[i_atom] == 'Co':
                debye_waller_factor = dbf_co # XXX
            else:
                exit(f"Error: unknown element {symbols[i_atom]}")
            print(symbols[i_atom], x, y, z, occupancy, debye_waller_factor, 0, 0, 0, file=cel_file)

    aberrations_values = np.zeros_like(aberrations)

    with open('wavimg.prm', "w") as f:
        print(f"'msa_{nz}.wav'", file=f)      # line 1
        print(f"{nx} {ny}", file=f)
        print(f"{0.1*box/nx} {0.1*box/ny}", file=f)
        print(f"{electron_energy}", file=f)
        print(f"{0}", file=f)               # line 5
        print(f"'image.dat'", file=f)
        print(f"{nx} {ny}", file=f)
        print(f"{0} {1.0} {1.0} {0.0}", file=f)
        print(f"{0}", file=f)
        print(f"{0.0}", file=f)            # line 10
        print(f"{0.0} {0.0}", file=f)
        print(f"{0.0}", file=f)
        print(f"{1}", file=f)
        print(f"{1} {3.8}", file=f)
        print(f"{1} {0.4}", file=f)        # line 15
        print(f"{0} {1.0} '{'mtf.dat'}'", file=f)
        print(f"{0} {1.0} {1.0} {0.0}", file=f)
        print(f"{len(aberrations)}", file=f)
        for i, abr in enumerate(aberrations):
            value = np.random.uniform(aberrations[i, 0], aberrations[i, 1])
            valx = value * np.cos(np.random.uniform(0, 2*np.pi))
            valy = value * np.sin(np.random.uniform(0, 2*np.pi))
            aberrations_values[i] = [valx, valy]
            print(f"{i} {valx} {valy}", file=f) 
        print(f"{250.0} {0.03}", file=f)   # line 19+Nabr
        print(f"{0.0} {0.0}", file=f)
        print(f"{0}", file=f)

    subprocess.run(["celslc", "-cel", "coord.cel", "-slc", "slice", "-nx", str(nx), "-ny", str(ny), "-nz", str(nz), "-ht", str(electron_energy), "-dwf", "-abs"])
    subprocess.run(["msa", "-prm", "msa.prm", "-out", f"msa.wav", "/ctem"])
    subprocess.run(["wavimg", "-prm", "wavimg.prm", "-out", "image.dat"])

    data = np.fromfile("image.dat", dtype=np.float32)
    counts = np.random.uniform(count_min, count_max)
    image = np.random.poisson(counts*data.reshape((nx, ny)))
    plt.imshow(image, cmap='gray')
    plt.axis('off')
    out_image = f"{os.path.basename(sys.argv[1]).split('.')[0]}_{i_variant}.png"
    plt.savefig(f"hrtem_images/{out_image}", bbox_inches='tight', pad_inches=0)
    plt.close()