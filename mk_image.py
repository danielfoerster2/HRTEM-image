#!/usr/bin/python3

import glob
import os
import subprocess
import sys

import matplotlib.pyplot    as plt
import numpy                as np


""" Parameters initialization """

process_id = os.getenv("process_id")
box = 100.0   # EN ANGSTROM
n_variants = 10 # Number of variants for each structure (random rotation, interatomic distance variation, aberrations)
dbf_ag = 0.019 # 0.008 XXX
dbf_co = 0.015 # 0.004 XXX
nx = 128
ny = 128
nz = int(os.getenv('nz'))
electron_energy = float(os.getenv('electron_energy'))
count_min = 100
count_max = 500
aberrations = np.array([[0.0, 0.0], # image shift
              [-4.0, -8.0],         # defocus
              [0.0, 0.0],           # astigmatism
              [0.0, 0.0],           # coma
              [0.0, 0.0],           # three-lobe aberration (Three-fold astigmatism ?)
              [0.001, 0.001],       # spherical aberration
              [0.0, 0.0]])          # star aberration

md_data = np.genfromtxt(f"tmp_{process_id}/tmp.dat", dtype=str)
first_column = np.array([row[0] for row in md_data]) # Extract the first column (id_sim)



""" Read the .xyz file and extract atomic positions and symbols """

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


for i_variant in range(n_variants):

    """ Random rotation and interatomic distance variation """

    phi = np.random.uniform(0, 2*np.pi)
    theta = np.random.uniform(0, np.pi)
    rot_matrix = np.array([[np.cos(phi)*np.cos(theta), -np.sin(phi), np.cos(phi)*np.sin(theta)],
                           [np.sin(phi)*np.cos(theta), np.cos(phi), np.sin(phi)*np.sin(theta)],
                           [-np.sin(theta), 0, np.cos(theta)]])
    pos = np.dot(pos_ini, rot_matrix) * np.random.uniform(0.95, 1.05)



    """ Computation of structural parameters """

    gyration_radius = np.sqrt(np.sum(np.linalg.norm(pos, axis=1)**2) / n_atoms)
    radius_limit = .8*gyration_radius

    pos_atm1 = pos[np.array(symbols) == 'Ag']
    pos_atm2 = pos[np.array(symbols) == 'Co']
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
    r_cm = np.sum(pos, axis=0) / n_atoms



    """ Creation of the .cel file for celslc program """

    atom_types = set(symbols)
    with open(f"tmp_{process_id}/coord.cel", 'w') as cel_file:
        res = ""
        for elem in atom_types:
            res += elem
        print(f"{res} {n_atoms}", file=cel_file)
        print(0, box, box, box, 90.0, 90.0, 90.0, file=cel_file)
        for i_atom in range(n_atoms):
            x = pos[i_atom][0]/box + 0.5
            y = pos[i_atom][1]/box + 0.5
            z = pos[i_atom][2]/box + 0.5
            occupancy = 1.0
            if symbols[i_atom] == 'Ag':
                debye_waller_factor = dbf_ag
            elif symbols[i_atom] == 'Co':
                debye_waller_factor = dbf_co
            else:
                exit(f"Error: unknown element {symbols[i_atom]}")
            print(symbols[i_atom], x, y, z, occupancy, debye_waller_factor, 0, 0, 0, file=cel_file)



    """ Creation of the .prm file for wavimg program """

    aberrations_values = np.zeros_like(aberrations)
    with open(f'tmp_{process_id}/wavimg.prm', "w") as f:
        print(f"'tmp_{process_id}/msa_sl0{nz}.wav'", file=f)      # line 1
        print(f"{nx} {ny}", file=f)
        print(f"{0.1*box/nx} {0.1*box/ny}", file=f)
        print(f"{electron_energy}", file=f)
        print(f"{0}", file=f)   # line 5
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
    


    """ DrProbe Simulation """

    subprocess.run(["celslc", "-cel", f"tmp_{process_id}/coord.cel", "-slc", f"tmp_{process_id}/slice", "-nx", str(nx), "-ny", str(ny), "-nz", str(nz), "-ht", str(electron_energy), "-dwf", "-abs"], stdout=subprocess.DEVNULL)
    subprocess.run(["rm", "-f", f"tmp_{process_id}/coord.cel"])
    subprocess.run(["msa", "-prm", f"tmp_{process_id}/msa.prm", "-out", f"tmp_{process_id}/msa.wav", "/ctem"], stdout=subprocess.DEVNULL)
    for fichier in glob.glob(f"tmp_{process_id}/*.sli"):
        os.remove(fichier)
    subprocess.run(["wavimg", "-prm", f"tmp_{process_id}/wavimg.prm", "-out", f"tmp_{process_id}/image.dat"], stdout=subprocess.DEVNULL)
    os.remove(f"tmp_{process_id}/msa_sl0{nz}.wav")



    """ Add Poisson noise to the image """

    data = np.fromfile(f"tmp_{process_id}/image.dat", dtype=np.float32)
    subprocess.run(["rm", "-f", f"tmp_{process_id}/image.dat"])
    counts = np.random.uniform(count_min, count_max)
    data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
    image = np.random.poisson(counts*data_normalized.reshape((nx, ny)))



    """ Save the HRTEM image """

    plt.imshow(image, cmap='gray')
    plt.axis('off')
    id_sim = f"{os.path.basename(sys.argv[1]).split('.')[0]}_{i_variant}"
    plt.savefig(f"hrtem_images/{id_sim}.png", bbox_inches='tight', pad_inches=0)
    plt.close()



    """ Save the simulation parameters """

    mask = first_column == os.path.basename(sys.argv[1]).split('.')[0]
    param = md_data[mask][0]
    f = open(f"tmp_{process_id}/data.dat", "a")
    print(id_sim, param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11], param[12], param[13], param[14], param[15], param[16], param[17], gyration_radius, nat1, nat2, nat1_out, nat2_out, nat1_in, nat2_in, r_cm1[0], r_cm1[1], r_cm1[2], r_cm2[0], r_cm2[1], r_cm2[2], r_cm[0], r_cm[1], r_cm[2], d_com, counts, phi, theta, aberrations_values[0][0], aberrations_values[0][1], aberrations_values[1][0], aberrations_values[1][1], aberrations_values[2][0], aberrations_values[2][1], aberrations_values[3][0], aberrations_values[3][1], aberrations_values[4][0], aberrations_values[4][1], aberrations_values[5][0], aberrations_values[5][1], aberrations_values[6][0], aberrations_values[6][1], sep="\t", file=f)
    f.close()
