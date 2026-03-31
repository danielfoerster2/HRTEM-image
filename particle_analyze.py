from ovito.io import import_file
from ovito.modifiers import *
from ovito.data import *

import glob
import os
import numpy as np

path_xyz = os.getenv("path_xyz")

md_data = np.genfromtxt(os.path.join(os.path.dirname(path_xyz), "d.dat"), dtype=None)
first_column = np.array([row[0] for row in md_data])

f = open("data.dat", "w")

for xyz_file in glob.glob("./xyz/*.xyz"):

    i_sim = os.path.basename(xyz_file).split(".")[0]
    mask = first_column == int(i_sim.split("_")[0])    
    n_steps = md_data[mask][0][3]
    initial_temperature = md_data[mask][0][4]

    node = import_file(xyz_file)

    node.modifiers.append(ConstructSurfaceModifier())
    node.modifiers.append(CommonNeighborAnalysisModifier())
    node.modifiers.append(BondAngleAnalysisModifier())
    node.modifiers.append(CentroSymmetryModifier())

    data = node.compute()
    pos = data.particle_properties.position.array
    n_atoms = pos.shape[0]
    epot = data.particle_properties['epot'].array
    box = data.cell.matrix.diagonal()

    surface_area = node.output.attributes['ConstructSurfaceMesh.surface_area']
    solid_volume = node.output.attributes['ConstructSurfaceMesh.solid_volume']
        
    cna_others = node.output.attributes['CommonNeighborAnalysis.counts.OTHER']
    cna_fcc = node.output.attributes['CommonNeighborAnalysis.counts.FCC']
    cna_hcp = node.output.attributes['CommonNeighborAnalysis.counts.HCP']
    cna_bcc = node.output.attributes['CommonNeighborAnalysis.counts.BCC']
    cna_ico = node.output.attributes['CommonNeighborAnalysis.counts.ICO']

    bond_angle_others = node.output.attributes['BondAngleAnalysis.counts.OTHER']
    bond_angle_fcc = node.output.attributes['BondAngleAnalysis.counts.FCC']
    bond_angle_hcp = node.output.attributes['BondAngleAnalysis.counts.HCP']
    bond_angle_bcc = node.output.attributes['BondAngleAnalysis.counts.BCC']
    bond_angle_ico = node.output.attributes['BondAngleAnalysis.counts.ICO']

    csp = np.mean(node.output.particle_properties['Centrosymmetry'].array)

    pos_atm1 = pos[data.particle_properties.particle_type.array == 1]
    pos_atm2 = pos[data.particle_properties.particle_type.array == 2]

    nat1 = pos_atm1.shape[0]
    nat2 = pos_atm2.shape[0]
        
    if nat1 != 0:
        r_cm1 = np.sum(pos_atm1, axis=0) / nat1
    else:
        r_cm1 = np.array([np.nan, np.nan, np.nan])
    if nat2 != 0:
        r_cm2 = np.sum(pos_atm2, axis=0) / nat2
    else:
        r_cm2 = np.array([np.nan, np.nan, np.nan])

    d_com = np.linalg.norm(r_cm1 - r_cm2)

    gyration_radius = np.sqrt(np.sum(np.linalg.norm(pos, axis=1)**2) / n_atoms)

    radius_limit = .8*gyration_radius
    mask_out_atm1 = np.linalg.norm(pos_atm1, axis=1) > radius_limit
    mask_out_atm2 = np.linalg.norm(pos_atm2, axis=1) > radius_limit
    
    nat1_out = np.sum(mask_out_atm1)
    nat2_out = np.sum(mask_out_atm2)
    nat1_in = np.sum(~mask_out_atm1)
    nat2_in = np.sum(~mask_out_atm2)

    print(i_sim, n_atoms, nat1, nat2, n_steps, initial_temperature, np.sum(epot), surface_area, solid_volume, cna_others, cna_fcc, cna_hcp, cna_bcc, cna_ico, bond_angle_others, bond_angle_fcc, bond_angle_hcp, bond_angle_bcc, bond_angle_ico, d_com, gyration_radius, nat1_out, nat2_out, nat1_in, nat2_in, csp, sep="\t", file=f)

f.close()