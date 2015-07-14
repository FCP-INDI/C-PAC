import networkx as nx
import random
from itertools import groupby

CONST_MAX_PASS = -1
EPSILON_GAIN = 0.0000001

def girvan_graphs(zout) :
	"""
	Create a graph of 128 vertices, 4 communities.
	Used as ground truth to compare to.
	Girvan newman, 2002. PNAS June, vol 99 n 12

	Community is node modulo 4
	"""

	pout = float(zout)/96.
	pin = (16.-pout*96.)/31.
	graph = nx.Graph()
	graph.add_nodes_from(range(128))
	for x in graph.nodes() :
		for y in graph.nodes() :
			if x < y :
				val = random.random()
				if x % 4 == y % 4 :
					#nodes belong to the same community
					if val < pin :
						graph.add_edge(x, y)

				else :
					if val < pout :
						graph.add_edge(x, y)
	return graph


class Community(object):
	"""Class representing communities
	   A community is here a disjoint set of nodes
	"""

	def __init__(self, graph):
		self.communities = dict(zip(graph, range(graph.number_of_nodes())))
		self.m = float(graph.size())
		self.degrees_per_node = graph.degree(graph.nodes())
		for k, v in self.degrees_per_node.items():
			self.degrees_per_node[k] = float(self.degrees_per_node[k])
		self.degrees_per_community = self.degrees_per_node.copy()
		# number of links within the community
		self.internals = {i:i-i for i in range(graph.number_of_nodes())}
		self.loops = {}

	def reinitialize(self, graph):
		self.communities = {}
		self.total_weight = 0
		self.m = 0.0
		self.degrees_per_node = {}
		self.degrees_per_community = {}
		self.internals = {}

		self.communities = dict(zip(graph, range(graph.number_of_nodes())))
		self.total_weight = float(graph.size(weight = 'weight'))
		self.m = self.total_weight
		self.degrees_per_node = graph.degree(graph.nodes(), weight='weight')
		for k, v in self.degrees_per_node.items():
			self.degrees_per_node[k] = float(self.degrees_per_node[k])
		self.degrees_per_community = self.degrees_per_node.copy()
		# number of links within the community
		self.internals = {i:i-i for i in range(graph.number_of_nodes())}
		
		count=0
		for node in graph.nodes():
			self.loops[node] = float(graph.get_edge_data(node, node, {"weight":0}).get("weight", 1))
			self.internals[count] = self.loops[node]
			count += 1
		# [(communities[node1], communities[node2], attr.get('weight', 1) ) for node1, node2, attr in graph.edges(data=True)]
		# graph.get_edge_data(node, node) for node, node in graph.edges(data=True)

		return self

	def get_all_nodes_of_community(self, community):
		"""
		Parameter community = int representing number of community
		"""
		nodelist = [node for node, com in self.communities.items() if community == com]
		return nodelist

	# def renumber_communites(self, communities):
	# 	vals = set(self.communities.values())
	# 	mapping = dict(zip(vals,range(len(vals))))

	# 	for key in communities.keys():
	# 		self.communities[key] = mapping[self.communities[key]]


def renumber_communites(communities):
	ret = communities.copy()
	vals = set(communities.values())
	mapping = dict(zip(vals,range(len(vals))))

	for key in communities.keys():
		ret[key] = mapping[communities[key]]
	
	return ret

def remove(node, community_to_move, com_ref, links):
	com_ref.degrees_per_community[community_to_move] = com_ref.degrees_per_community.get(community_to_move, 0.0) - com_ref.degrees_per_node.get(node, 0.0)
	com_ref.internals[community_to_move] = float(com_ref.internals.get(community_to_move, 0.0) - links - com_ref.loops.get(node, 0.0))	
	com_ref.communities[node] = -1

def insert(node, community_to_move, com_ref, links):
	com_ref.communities[node] = community_to_move
	com_ref.degrees_per_community[community_to_move] = (com_ref.degrees_per_community.get(community_to_move, 0.0) + com_ref.degrees_per_node.get(node, 0.0))
	com_ref.internals[community_to_move] = float(com_ref.internals.get(community_to_move, 0.0) + links + com_ref.loops.get(node, 0.0))

def neighbourhood_link_strength(node, graph, com_ref):
	community_links = {}
	for neighbour in graph.neighbors(node):
		community_of_neighbour = com_ref.communities[neighbour]
		community_links[community_of_neighbour] = community_links.get(community_of_neighbour, 0) + 1
	return community_links

def compute_dendrogram(graph, part_init=None):
	current_graph = graph.copy()
	com_ref = Community(current_graph)
	mod = modularity(current_graph, com_ref)
	community_list = list()
	first_pass(current_graph, com_ref)
	new_mod = modularity(current_graph, com_ref)
	partition = renumber_communites(com_ref.communities)
	community_list.append(partition)
	mod = new_mod
	current_graph = second_pass(partition, current_graph)
	com_ref = com_ref.reinitialize(current_graph)

	while True:
		first_pass(current_graph, com_ref)
		new_mod = modularity(current_graph, com_ref)
		if new_mod - mod < EPSILON_GAIN :
			break
		partition = renumber_communites(com_ref.communities)
		community_list.append(partition)
		mod = new_mod
		current_graph = second_pass(partition, current_graph)
		com_ref = com_ref.reinitialize(current_graph)
	return community_list[:]


def first_pass(graph, com_ref):
	enhancement_possible = True
	current_modularity = modularity(graph, com_ref)
	new_modularity = current_modularity
	iterations = 0

	while enhancement_possible  and iterations != CONST_MAX_PASS :
		current_modularity = new_modularity
		enhancement_possible = False
		iterations += 1

		for node in graph.nodes():
			node_community = com_ref.communities[node]
			degc_totw = com_ref.degrees_per_node.get(node, 0.0) / (com_ref.m*2.0)
			neigh_communities =  neighbourhood_link_strength(node, graph, com_ref)
			remove(node, node_community, com_ref, neigh_communities.get(node_community, 0.0))
			best_com = node_community
			best_increase = 0
			for com, dnc in neigh_communities.items():
				incr =  dnc  - com_ref.degrees_per_community.get(com, 0.) * degc_totw
				if incr > best_increase:
					best_increase = incr
					best_com = com
			insert(node, best_com, com_ref, neigh_communities.get(best_com, 0.))
			if best_com != node_community:
				enhancement_possible = True
		new_modularity = modularity(graph, com_ref)
		if new_modularity - current_modularity < EPSILON_GAIN :
			break


def second_pass(communities, graph):
	""" Nodes belonging to the same community as detected in pass 1 are merged into a single node.
		The new graph is build up from these so called "super nodes"
	"""
	aggregated_graph = nx.Graph()

	# The new graph consists of as many "supernodes" as there are communities
	aggregated_graph.add_nodes_from(set(communities.values()))
	# make edges between communites, bundle more edges between nodes in weight attribute
	edge_list=[(communities[node1], communities[node2], attr.get('weight', 1) ) for node1, node2, attr in graph.edges(data=True)]
	sorted_edge_list = sorted(edge_list)
	sum_z = lambda tuples: sum(t[2] for t in tuples)
	weighted_edge_list = [(k[0], k[1], sum_z(g)) for k, g in groupby(sorted_edge_list, lambda t: (t[0], t[1]))]
	aggregated_graph.add_weighted_edges_from(weighted_edge_list)

	return aggregated_graph


	

def modularity(graph, communityRef):
	"""
	Compute the modularity of the graph using:

	.. math::

	M = \sum_{c=1}^{n_C} \left [ \frac{L_c}{L} - \left (  \frac{k_c}{2L} \right ) \right ]

	where L_c is the total number of links within the community C and k_c is the total degree of the noes in this community

	References:
	-----------
	.. [1] M. E. J. Newman and M. Girvan (2008).
		   Finding and evaluating community structure in networks.	
		   Phys. Rev. E 69, 026113
	"""
	q = 0.0
	m = float(communityRef.m)
	for community in set(communityRef.communities.values()):
		L_c = communityRef.internals.get(community, 0.0)
		K_c = communityRef.degrees_per_community.get(community, 0.0)
		q += (L_c / m) - ((K_c / (2.0*m))** 2.)
	return q


def delta_q(graph, community, node):
	"""Compute the gain of modularity delta_Q if node would be inserted in community

	 .. math::

		delta_Q = \left [  \frac{\sum_{in}+2k_{i_in}}{2m} -  \left (\frac{\sum_{tot+k_i}}{2m}  \right )^{2}   \right ] -  \left [  \frac{\sum_{in}}{2m} -  \left (\frac{\sum_{tot}}{2m}  \right )^{2} - \left ( \frac{k_i}{2m}  \right )^{2} \right ]

	References
	----------
	..[1] Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). 
		  Fast unfolding of communities in large networks. 
		  Journal of Statistical Mechanics: Theory and Experiment, 10(1), 008. http://doi.org/10.1088/1742-5468/2008/10/P10008
	"""	
	nodes = Community.get_all_nodes_of_community(community)

	#number of links inside the community
	sigma_in  = sum(graph.degree((nodes)).values())
	#number of inks from node i to other nodes in community
	k_i_in    = len(set(graph.neighbours(node)).intersection(nodes))
	#number links incident to nodes in community
	sigma_tot = len(set(sum([graph.neighbors(node) for node in nodes], [])) )
	#number links incident to node i
	k_i       = len(graph.neighbours(node))
	m         = graph.size()

	term1 = (((sigma_in+k_i_in)/2.0*m)-((sigma_tot+k_i)/(2.0*m))**2.0)
	term2 = sigma_in/(2.0*m)-(sigma_tot/(2.0*m)**2.0)-(k_i/(2.0*m)**2.0)
	delta_q = term1 - term2
	return delta_q


if __name__ == '__main__':
	graph = girvan_graphs(4)
	community = Community(graph)
	delta_q(graph, 2, 4)
 