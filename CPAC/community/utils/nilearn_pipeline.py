import numpy as np
import pandas as pd
from nilearn import datasets
from nilearn.input_data import NiftiMapsMasker
from nilearn import plotting
from matplotlib import pyplot as plt
from nilearn import image

import networkx as nx
import sys
sys.path.append("..")
import louvain

"""
Use MDL probabilistic atlas to have regions of functional structures
https://team.inria.fr/parietal/research/spatial_patterns/spatial-patterns-in-resting-state/
References
----------
[1] Learning and comparing functional connectomes across subjects.
     Gael Varoquaux, R.C. Craddock NeuroImage, 2013
"""


number_subjects=15


def run_mini_pipeline():
    atlas = datasets.fetch_atlas_msdl()
    atlas_img = atlas['maps']
    labels = pd.read_csv(atlas['labels'])['name']

    masker = NiftiMapsMasker(maps_img=atlas_img, standardize=True,
                               memory='/tmp/nilearn', verbose=0)

    data = datasets.fetch_adhd(number_subjects)

    figures_folder = '../figures/'
    count=0
    for func_file, confound_file in zip(data.func, data.confounds):
        
        # fit the data to the atlas mask, regress out confounds
        time_series = masker.fit_transform(func_file, confounds=confound_file)

        correlation = np.corrcoef(time_series.T)

        #plotting starts here
        plt.figure(figsize=(10, 10))
        plt.imshow(correlation, interpolation="nearest")
        x_ticks = plt.xticks(range(len(labels)), labels, rotation=90)
        y_ticks = plt.yticks(range(len(labels)), labels)
        corr_file = figures_folder+'subject_number_' + str(count) + '_correlation.pdf'
        plt.savefig(corr_file)

        atlas_region_coords = [plotting.find_xyz_cut_coords(img) for img in image.iter_img(atlas_img)]
        threshold = 0.6
        plotting.plot_connectome(correlation, atlas_region_coords, edge_threshold=threshold)
        connectome_file = figures_folder+'subject_number_' + str(count) + '_connectome.pdf'
        plt.savefig(connectome_file)


        #graph setup

        #binarize correlation matrix
        correlation[correlation<threshold] = 0
        correlation[correlation != 0] = 1

        graph = nx.from_numpy_matrix(correlation)

        partition=louvain.best_partition(graph)

        values = [partition.get(node) for node in graph.nodes()]

        plt.figure()
        nx.draw_spring(graph, cmap = plt.get_cmap('jet'), node_color = values, node_size=30, with_labels=True)
        graph_file = figures_folder+'subject_number_' + str(count) + '_community.pdf'
        plt.savefig(graph_file)

        count += 1

        plt.close('all')


if __name__ == '__main__':
    print "running mini pipeline with {number} subjects".format(number=number_subjects)
    run_mini_pipeline()

