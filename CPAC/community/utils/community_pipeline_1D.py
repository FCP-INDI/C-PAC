import numpy as np
import pandas as pd
import networkx as nx
import os
import sys
sys.path.append("louvain")
import louvain

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import logging

global_frame_list = []

def run_mini_pipeline():
	datadir = 'data'
	for root, dirs, filenames in os.walk(datadir):
		for file in filenames:
			print 'computing subject ' + str(filenames.index(file))
			logging.info('computing subject ' + str(filenames.index(file)))
			file_handle = os.path.join(root, file)
			df = pd.read_csv(file_handle, sep='\t')
			np_ts_matrix=df.as_matrix()
			correlation_matrix = np.corrcoef(np_ts_matrix.T)
			threshold = 0.7
			correlation_matrix[correlation_matrix<threshold] = 0
			correlation_matrix[correlation_matrix != 0] = 1
			graph = nx.from_numpy_matrix(correlation_matrix)
			print graph.size()
			print graph.number_of_nodes()
			print graph.number_of_edges()
			# nx.write_edgelist(graph, 'logs/test_edgelist' + str(filenames.index(file)))
			if filenames.index(file) != float('inf'):
				save_name = 'subject ' + str(filenames.index(file))
				analyse_community(graph, save_name)
		print 'done'
		logging.info('done')


def analyse_community(graph, save_name):
	partition = louvain.best_partition(graph)
	dendro = louvain.generate_dendrogram(graph)
	print 'Number of communities: ' + str(len(set(partition.values())))
	logging.info('Number of communities: ' + str(len(set(partition.values()))))
	com_node_map = {}
	for community in set(partition.values()):
		com_node_map[community] = len([nodes for nodes in partition.keys() if partition[nodes] == community])
	for k, v in com_node_map.items():
		print 'Community ' + str(k+1) + ' has ' + str(v) + ' nodes' 
		logging.info('Community ' + str(k+1) + ' has ' + str(v) + ' nodes')
	nx.set_node_attributes(graph, 'community', partition)
	no_pe_com_map_list = nodes_per_community(dendro)
	community_size_distribution(no_pe_com_map_list, save_name)
	drawNetwork(graph, save_name)



def community_size_distribution(comm_map, save_name):
	from collections import Counter
	import pandas as pd
	import seaborn as sns
	for level in range(len(comm_map)):
		# sizes=comm_map[level].values()
		# distribution=pd.DataFrame(Counter(sizes).items(), columns=['community size','number of communities'])
		# plt.clf()
		# distribution.plot(kind='scatter', x='community size', y='number of communities', title='Level: ' + str(level))
		# plt.savefig('figures/' + save_name + '_community_size_distribution_level_' + str(level))
		df = pd.DataFrame([(level,val) for level, dct in enumerate(comm_map) for val in dct.values()], columns=['level', 'size'])
		size_count = df.groupby(['level'])['size'].apply(lambda x: x.value_counts())
		size_count = size_count.reset_index()
		size_count.columns = ['level', 'community size', 'number of communities']
		global_frame_list.append(size_count)
		cmap = plt.get_cmap('jet')
		plt.clf()
		size_count.plot(kind='scatter', x='community size', y='number of communities', s=100, c='level', cmap=cmap)
		plt.savefig('figures/' + save_name + '_community_size_distribution_level_' + str(level))
	plt.clf()
	plt.close()
	# subject_frame = pd.concat(frames, ignore_index=True)
	



def nodes_per_community(dendrogram):
    result = list()
    com_node_map = dict()
    for level in range(len(dendrogram)):
        current_partition = dendrogram[level]
        for community in set(current_partition.values()):
            com_node_map[community] = len([nodes for nodes in current_partition.keys() if current_partition[nodes] == community])
        result.append(com_node_map)
        com_node_map = {}
    return result


def drawNetwork(G, save_name):
	# position map
	pos = nx.spring_layout(G)
	# community ids
	communities = [v for k,v in nx.get_node_attributes(G, 'community').items()]
	numCommunities = max(communities) + 1
	# color map from http://colorbrewer2.org/
	cmapLight = colors.ListedColormap(['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], 'indexed', numCommunities)
	cmapDark = colors.ListedColormap(['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], 'indexed', numCommunities)

	# edges
	nx.draw_networkx_edges(G, pos)

	# nodes
	nodeCollection = nx.draw_networkx_nodes(G,
		pos = pos,
		node_color = communities,
		cmap = cmapLight
	)
	# set node border color to the darker shade
	darkColors = [cmapDark(v) for v in communities]
	nodeCollection.set_edgecolor(darkColors)

	# Print node labels separately instead
	for n in G.nodes_iter():
		plt.annotate(n,
			xy = pos[n],
			textcoords = 'offset points',
			horizontalalignment = 'center',
			verticalalignment = 'center',
			xytext = [0, 2],
			color = cmapDark(communities[n])
		)

	plt.axis('off')
	plt.savefig('figures/' + save_name + '_community_network')
	# plt.show()
	plt.clf()
	plt.close()

if __name__ == '__main__':
	logging.basicConfig(filename='logs/community.log', level=logging.INFO, format="%(asctime)-15s %(levelname)-8s %(message)s", filemode='w')
	run_mini_pipeline()
	global_frame = pd.concat(global_frame_list, ignore_index=True)
	mean_result=global_frame.groupby(['level', 'community size']).mean().reset_index()
	mean_result.columns = ['level', 'community size', 'average number of communities']
	mean_result.to_csv('logs/' + datadir + 'mean_result.csv')
	print mean_result