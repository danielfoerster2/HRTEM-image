from ovito.io import import_file, export_file
from ovito.modifiers import ClusterAnalysisModifier
from ovito.data import *

import os
import glob
import numpy as np

process_id = os.getenv("process_id")
path_xyz = os.getenv("path_xyz")

files = glob.glob(os.path.join(path_xyz, "*" + process_id + ".xyz"))

for xyz_file in files:
    id_sim = int(os.path.basename(xyz_file).split(".")[0])
    node = import_file(path_xyz + "/" + str(id_sim) + ".xyz")
    cluster_mod = ClusterAnalysisModifier(cutoff=4.085 + 4.085*.1, sort_by_size=True)
    node.modifiers.append(cluster_mod)

    data = node.compute()
    pos = data.particle_properties.position.array
    epot = data.particle_properties['epot'].array
    cell = data.cell
    box = cell.matrix.diagonal()
    box_center = np.array([box[0]/2, box[1]/2, box[2]/2])
    clusters = data.particle_properties.cluster.array
    types = data.particle_properties.particle_type.array
    typ_list = data.particle_properties.particle_type.type_list
    symbols = np.zeros(pos.shape[0], dtype="U2")
    for i in range(pos.shape[0]):
        symbols[i] = typ_list[types[i]-1].name
    n_clusters = np.unique(data.particle_properties['Cluster'].array).size

    nanoparticle_id = 1
    for cluster_id in range(1, n_clusters+1):
        mask = clusters == cluster_id
        size_cluster = np.sum(mask)
        if size_cluster < 10:
            continue
        pos_cluster = pos[mask]
        symbols_cluster = symbols[mask]
        epot_cluster = epot[mask]
        pos_cluster = pos_cluster - box * np.round((pos_cluster - pos_cluster[0]) / box)
        pos_cluster = pos_cluster - np.sum(pos_cluster, axis=0)/size_cluster
        pos_cluster = pos_cluster - box * np.round(pos_cluster / box)
        pos_cluster = pos_cluster - np.sum(pos_cluster, axis=0)/size_cluster

        with open("./xyz/%i_%i.xyz" %(nanoparticle_id, id_sim), "w") as file:
            print(size_cluster, file=file)
            print('Lattice="%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f" Properties=species:S:1:pos:R:3:epot:R:1' % (box[0], box[1], box[2]), file=file)
            for i_atom in range(size_cluster):
                print(symbols_cluster[i_atom], pos_cluster[i_atom][0], pos_cluster[i_atom][1], pos_cluster[i_atom][2], epot_cluster[i_atom], file=file)
        nanoparticle_id += 1
