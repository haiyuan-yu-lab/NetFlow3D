description='''
- Functions used by Net3D.py
- Dependencies: None
'''

import scipy.stats
import numpy as np
import os
from collections import defaultdict,Counter
import pandas as pd
import networkx as nx
import random
import copy
from os import path
import logging
import shutil
from networkx.algorithms.community.quality import modularity
import heapq


def calc_enr(c1, n1, c2, n2):
	"""
	Parameters
	----------
	c1: integer
		The count of instance A of a specific quality (e.g. mutations in a protein).
	n1: integer
		The total count of instance A (e.g. mutations in all proteins).
	c2: integer
		The count of instance B of a specific quality (e.g. length of a protein).
	n2: integer
		The total count of instance B (e.g. total length of all proteins).

	Returns
	-------
	enr: float
		Enrichment value, np.nan if invalid
	SE_enr: float
		Standard error of the log odds ratio, np.nan if invalid
	pval: float
		Two-tailed p-value for the null hypothesis that the log odds ratio is 0.
	"""
	if c1 == 0 or c2 == 0 or n1 == 0 or n2 == 0 or c1 == n1 and c2 == n2:
		return [np.nan, np.nan, np.nan]
	p1 = c1 / n1
	p2 = c2 / n2
	enr = p1 / p2
	SE_log_enr = np.sqrt(1.0 / c1 - 1.0 / n1 + 1.0 / c2 - 1.0 / n2) # Calculated by the delta method.
	SE_enr = SE_log_enr * enr
	z = np.log(enr) / SE_log_enr
	pval = scipy.stats.norm.sf(abs(z)) * 2
	return [enr, SE_enr, pval]


def get_uniprot(row, gene2uniprot, ensp2uniprot, enst2uniprot, ensg2uniprot):
	if type(row["ENSP"]) == str and row["ENSP"].split(".")[0] in ensp2uniprot:
		return ensp2uniprot[row["ENSP"].split(".")[0]]
	elif ("Transcript_ID" in row) and type(row["Transcript_ID"]) == str and row["Transcript_ID"].split(".")[0] in enst2uniprot:
		return enst2uniprot[row["Transcript_ID"].split(".")[0]][0]
	elif ("Gene" in row) and type(row["Gene"]) == str and row["Gene"].split(".")[0] in ensg2uniprot:
		return ensg2uniprot[row["Gene"].split(".")[0]][0]
	elif type(row["Hugo_Symbol"]) == str and row["Hugo_Symbol"] in gene2uniprot:
		return gene2uniprot[row["Hugo_Symbol"]]
	else:
		return np.nan


def df_to_dict(df_file, header = None, from_col = None, to_col = None):
	def df_to_dict(row, from_col, to_col):
		df2dict[row[from_col]] = row[to_col]
		return True
	df = pd.read_csv(df_file, sep = "\t", header = header)
	df2dict = {}
	df["status"] = df.apply(lambda row: df_to_dict(row, from_col, to_col), axis = 1)
	return df2dict


def split_consecutive_pos(pos):
	try:
		if pos == "-":
			return []
		elif pos.find('-') == -1:
			return [int(pos)]
		else:
			init, end = pos.split('-')
			if init == '?':
				return [int(end)]
			if end == '?':
				return [int(init)]
			return list(range(int(init), int(end) + 1))
	except:
		logging.warning("Unknown protein position: " + pos)
		return []


def get_32_trinuc(trinuc, complement = {"A":"T", "T":"A", "C":"G", "G":"C"}):
	if trinuc[1] in {"A", "G"}:
		output = "".join([complement[element] for element in trinuc][::-1])
	else:
		output = trinuc
	return output


def get_res_annotation(df, prolen_dict, output_path):
	"""
	Generate ShortVersion_mutation_data.txt, MUg.txt, and Patient_list.txt
	"""
	def gather_mutation_info(row, prolen):
		ref = row["Codons"].split("/")[0].upper()
		pos = split_consecutive_pos(row["Protein_position"])
		if len(pos) > 0 and len(ref) == 3 * len(pos):
			for i, element in enumerate(pos):
				res2mutcount[(row["UniProt"], element)] = res2mutcount[(row["UniProt"], element)] + 1
				if element <= prolen:
					uniprot2mutcount[row["UniProt"]][element] = uniprot2mutcount[row["UniProt"]][element] + 1
				
				ref_tmp = ref[(3*i):(3*i+3)]
				res2trinuc[(row["UniProt"], element)].append(ref_tmp)
				res2mutinfo[(row["UniProt"], element)].append(":".join([row["Tumor_Sample_Barcode"], row["Hugo_Symbol"], ref_tmp]))
		return True

	# gather mutation information
	res2mutcount = defaultdict(int)
	uniprot2mutcount = defaultdict(lambda: defaultdict(int))
	res2trinuc = defaultdict(list)
	res2mutinfo = defaultdict(list)
	df["status"] = df.apply(lambda row: gather_mutation_info(row, prolen_dict[row["UniProt"]]), axis = 1)
	uniprot2gene = df[["UniProt", "Hugo_Symbol"]].drop_duplicates().groupby("UniProt").count()

	# ShortVersion_mutation_data.txt
	f = open(output_path + "ShortVersion_mutation_data.txt", "w")
	f.write("\t".join(["Uniprot", "Uniprot_Position", "Mutations", "Ref_Tri", "Mutation_info"]) + "\n")
	for key in res2mutcount:
		trinuc_dict = dict(Counter(res2trinuc[key]))
		f.write("\t".join([key[0], str(key[1]), str(res2mutcount[key]), get_32_trinuc(sorted(trinuc_dict, key = trinuc_dict.get)[-1]), \
			",".join(res2mutinfo[key])]) + "\n")
	f.close()

	# Patient_list.txt
	sample_list = df["Tumor_Sample_Barcode"].unique().tolist()
	f = open(output_path + "Patient_list.txt", "w")
	f.write("#Tumor_Sample_Barcode" + "\n" + "\n".join(sample_list) + "\n")
	f.close()

	# MUg.txt
	f = open(output_path + "MUg.txt", "w")
	f.write("\t".join(["#Uniprot", "MUg"]) + "\n")
	for key in df["UniProt"].unique():
		if key in uniprot2mutcount:
			quantile = float(pd.DataFrame(list(uniprot2mutcount[key].values())).quantile(q = 0.99, axis = 0))
			cutoff = max(quantile, 20 * uniprot2gene.loc[key]["Hugo_Symbol"])
			i = 0
			j = 0
			for res in uniprot2mutcount[key]:
				if uniprot2mutcount[key][res] > cutoff:
					i = i + 1
					j = j + uniprot2mutcount[key][res]
			mut_freq = (sum(uniprot2mutcount[key].values()) - j) / (len(sample_list) * uniprot2gene.loc[key]["Hugo_Symbol"] * (prolen_dict[key] - i))
			f.write("\t".join([key, str(mut_freq)]) + "\n")
		else:
			f.write("\t".join([key, str(0.1)]) + "\n")
	f.close()

	return [output_path + "ShortVersion_mutation_data.txt", output_path + "Patient_list.txt", output_path + "MUg.txt"]


def get_res_mutability(shortversion, uniprot2gene, gene2mutability, uniprot2mutability, trinuc2mutability, output):
	df = pd.read_csv(shortversion, sep = "\t")
	mg = {}
	for element in df["Uniprot"].unique():
		gene_mut = []
		for item in uniprot2gene[element]:
			if item in gene2mutability:
				gene_mut.append(gene2mutability[item])
		if len(gene_mut) > 0:
			mg[element] = np.mean(gene_mut)
		else:
			mg[element] = np.mean(list(gene2mutability.values()))
	df["res_mutability"] = df.apply(lambda row: trinuc2mutability[row["Ref_Tri"]] / mg[row["Uniprot"]] * uniprot2mutability[row["Uniprot"]], axis = 1)
	df.to_csv(output, sep = "\t", index = None, header = True)
	return output


def get_graph_by_uniprot_binary(uniprot, uniprot2res, output_path, resource_dir_intra, resource_dir_inter = None, dist_intra = np.inf, \
	dist_inter = np.inf):
	"""
	Given all the mutated residues in a protein, generate intra- and inter-protein MR-MR contact network. (MR: mutated residue)
	"""
	output_list = []
	if path.exists(resource_dir_intra + uniprot + ".graphml.gz"):
		# load the intra-chain residue-residue distance map
		G = nx.read_graphml(resource_dir_intra + uniprot + ".graphml.gz")
		# get residues with mutation(s)
		res_list = []
		for item in uniprot2res[uniprot]:
			if uniprot + "_" + str(item) in G:
				res_list.append(uniprot + "_" + str(item))
		if len(res_list) > 0:
			G_sub = nx.Graph(G.subgraph(res_list))
			for element in list(G_sub.edges()):
				if G_sub[element[0]][element[1]]['distance'] > dist_intra:
					G_sub.remove_edge(element[0], element[1])
				else:
					G_sub[element[0]][element[1]]['type'] = 'intra'
			nx.write_graphml(G_sub, output_path + uniprot + ".graphml")
			output_list.append(output_path + uniprot + ".graphml")

		if resource_dir_inter != None and path.exists(resource_dir_inter + uniprot + ".graphml.gz"):
			# load the inter-chain residue-residue distance map
			G_inter = nx.read_graphml(resource_dir_inter + uniprot + ".graphml.gz")
			# get residues with mutation(s)
			res_list_inter = []
			for element in G_inter.nodes():
				element = element.split("_")
				if int(element[1]) in uniprot2res[element[0]]:
					res_list_inter.append("_".join(element))
			if len(res_list_inter) > 0:
				G_inter_sub = nx.Graph(G_inter.subgraph(res_list_inter))
				for element in list(G_inter_sub.edges()):
					if G_inter_sub[element[0]][element[1]]['distance'] > dist_inter:
						G_inter_sub.remove_edge(element[0], element[1])
					else:
						G_inter_sub[element[0]][element[1]]['type'] = 'inter'

				# get all interactors of the input protein
				interactor = set([item.split("_")[0] for item in res_list_inter])
				for element in interactor:
					if path.exists(output_path + "_".join(sorted([element, uniprot])) + ".graphml") == False:
						# get all edges connecting the input uniprot and this interactor
						edge_list = []
						for item in G_inter_sub.edges():
							if (item[0].split("_")[0], item[1].split("_")[0]) in {(uniprot, element), (element, uniprot)}:
								edge_list.append(item)
						if len(edge_list) > 0:
							G_tmp = G_inter_sub.edge_subgraph(edge_list)
							G_final = nx.Graph()
							G_final.add_edges_from(G_tmp.edges(data = True))
							G_final.add_nodes_from(G_tmp.nodes(data = True))

							# for each binary interaction, combine intra- and inter-chain edges
							if len(res_list) > 0:
								G_final.add_edges_from(G_sub.edges(data = True))
								G_final.add_nodes_from(G_sub.nodes(data = True))

							# for a heterodimer, the intra-protein edges in the interactor shoud be added
							if element != uniprot and path.exists(resource_dir_intra + element + ".graphml.gz"):
								G_interactor = nx.read_graphml(resource_dir_intra + element + ".graphml.gz")
								# get nodes with mutation(s)
								res_list_interactor = []
								for item in uniprot2res[element]:
									if element + "_" + str(item) in G_interactor:
										res_list_interactor.append(element + "_" + str(item))
								if len(res_list_interactor) > 0:
									G_interactor_sub = nx.Graph(G_interactor.subgraph(res_list_interactor))
									for item in list(G_interactor_sub.edges()):
										if G_interactor_sub[item[0]][item[1]]["distance"] > dist_intra:
											G_interactor_sub.remove_edge(item[0], item[1])
										else:
											G_interactor_sub[item[0]][item[1]]["type"] = "intra"
									G_final.add_edges_from(G_interactor_sub.edges(data = True))
									G_final.add_nodes_from(G_interactor_sub.nodes(data = True))
							nx.write_graphml(G_final, output_path + "_".join(sorted([element, uniprot])) + ".graphml")
							output_list.append(output_path + "_".join(sorted([element, uniprot])) + ".graphml")
	return output_list


def format_graph(G):
	"""
	Sort the nodes and edges in a graph.
	"""
	H = nx.Graph()
	H.add_nodes_from(sorted(G.nodes(data = True)))
	H.add_edges_from(sorted([(sorted([u,v])[0], sorted([u,v])[1], info) for u,v,info in G.edges(data = True)]))
	return H


class MappedQueue:
	"""The MappedQueue class implements an efficient minimum heap. The
	smallest element can be popped in O(1) time, new elements can be pushed
	in O(log n) time, and any element can be removed or updated in O(log n)
	time. The queue cannot contain duplicate elements and an attempt to push an
	element already in the queue will have no effect.

	MappedQueue complements the heapq package from the python standard
	library. While MappedQueue is designed for maximum compatibility with
	heapq, it has slightly different functionality.

	Examples
	--------

	A `MappedQueue` can be created empty or optionally given an array of
	initial elements. Calling `push()` will add an element and calling `pop()`
	will remove and return the smallest element.

	>>> q = MappedQueue([916, 50, 4609, 493, 237])
	>>> q.push(1310)
	True
	>>> x = [q.pop() for i in range(len(q.h))]
	>>> x
	[50, 237, 493, 916, 1310, 4609]

	Elements can also be updated or removed from anywhere in the queue.

	>>> q = MappedQueue([916, 50, 4609, 493, 237])
	>>> q.remove(493)
	>>> q.update(237, 1117)
	>>> x = [q.pop() for i in range(len(q.h))]
	>>> x
	[50, 916, 1117, 4609]

	References
	----------
	.. [1] Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2001).
	   Introduction to algorithms second edition.
	.. [2] Knuth, D. E. (1997). The art of computer programming (Vol. 3).
	   Pearson Education.
	"""

	def __init__(self, data=[]):
		"""Priority queue class with updatable priorities.
		"""
		self.h = list(data)
		self.d = dict()
		self._heapify()

	def __len__(self):
		return len(self.h)

	def _heapify(self):
		"""Restore heap invariant and recalculate map."""
		heapq.heapify(self.h)
		self.d = {elt: pos for pos, elt in enumerate(self.h)}
		if len(self.h) != len(self.d):
			raise AssertionError("Heap contains duplicate elements")

	def push(self, elt):
		"""Add an element to the queue."""
		# If element is already in queue, do nothing
		if elt in self.d:
			return False
		# Add element to heap and dict
		pos = len(self.h)
		self.h.append(elt)
		self.d[elt] = pos
		# Restore invariant by sifting down
		self._siftdown(pos)
		return True

	def pop(self):
		"""Remove and return the smallest element in the queue."""
		# Remove smallest element
		elt = self.h[0]
		del self.d[elt]
		# If elt is last item, remove and return
		if len(self.h) == 1:
			self.h.pop()
			return elt
		# Replace root with last element
		last = self.h.pop()
		self.h[0] = last
		self.d[last] = 0
		# Restore invariant by sifting up, then down
		pos = self._siftup(0)
		self._siftdown(pos)
		# Return smallest element
		return elt

	def update(self, elt, new):
		"""Replace an element in the queue with a new one."""
		# Replace
		pos = self.d[elt]
		self.h[pos] = new
		del self.d[elt]
		self.d[new] = pos
		# Restore invariant by sifting up, then down
		pos = self._siftup(pos)
		self._siftdown(pos)

	def remove(self, elt):
		"""Remove an element from the queue."""
		# Find and remove element
		try:
			pos = self.d[elt]
			del self.d[elt]
		except KeyError:
			# Not in queue
			raise
		# If elt is last item, remove and return
		if pos == len(self.h) - 1:
			self.h.pop()
			return
		# Replace elt with last element
		last = self.h.pop()
		self.h[pos] = last
		self.d[last] = pos
		# Restore invariant by sifting up, then down
		pos = self._siftup(pos)
		self._siftdown(pos)

	def _siftup(self, pos):
		"""Move element at pos down to a leaf by repeatedly moving the smaller
		child up."""
		h, d = self.h, self.d
		elt = h[pos]
		# Continue until element is in a leaf
		end_pos = len(h)
		left_pos = (pos << 1) + 1
		while left_pos < end_pos:
			# Left child is guaranteed to exist by loop predicate
			left = h[left_pos]
			try:
				right_pos = left_pos + 1
				right = h[right_pos]
				# Out-of-place, swap with left unless right is smaller
				if right < left:
					h[pos], h[right_pos] = right, elt
					pos, right_pos = right_pos, pos
					d[elt], d[right] = pos, right_pos
				else:
					h[pos], h[left_pos] = left, elt
					pos, left_pos = left_pos, pos
					d[elt], d[left] = pos, left_pos
			except IndexError:
				# Left leaf is the end of the heap, swap
				h[pos], h[left_pos] = left, elt
				pos, left_pos = left_pos, pos
				d[elt], d[left] = pos, left_pos
			# Update left_pos
			left_pos = (pos << 1) + 1
		return pos

	def _siftdown(self, pos):
		"""Restore invariant by repeatedly replacing out-of-place element with
		its parent."""
		h, d = self.h, self.d
		elt = h[pos]
		# Continue until element is at root
		while pos > 0:
			parent_pos = (pos - 1) >> 1
			parent = h[parent_pos]
			if parent > elt:
				# Swap out-of-place element with parent
				h[parent_pos], h[pos] = elt, parent
				parent_pos, pos = pos, parent_pos
				d[elt] = pos
				d[parent] = parent_pos
			else:
				# Invariant is satisfied
				break
		return pos


def greedy_modularity_communities(G, weight=None):
	"""Find communities in graph using Clauset-Newman-Moore greedy modularity
	maximization. This method currently supports the Graph class.

	Greedy modularity maximization begins with each node in its own community
	and joins the pair of communities that most increases modularity until no
	such pair exists.

	Parameters
	----------
	G : NetworkX graph

	Returns
	-------
	Yields sets of nodes, one for each community.

	Examples
	--------
	>>> from networkx.algorithms.community import greedy_modularity_communities
	>>> G = nx.karate_club_graph()
	>>> c = list(greedy_modularity_communities(G))
	>>> sorted(c[0])
	[8, 14, 15, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]

	References
	----------
	.. [1] M. E. J Newman 'Networks: An Introduction', page 224
	   Oxford University Press 2011.
	.. [2] Clauset, A., Newman, M. E., & Moore, C.
	   "Finding community structure in very large networks."
	   Physical Review E 70(6), 2004.
	"""

	# Count nodes and edges
	N = len(G.nodes())
	m = sum([d.get("weight", 1) for u, v, d in G.edges(data=True)])
	q0 = 1.0 / (2.0 * m)

	# Map node labels to contiguous integers
	label_for_node = {i: v for i, v in enumerate(G.nodes())}
	node_for_label = {label_for_node[i]: i for i in range(N)}

	# Calculate degrees
	k_for_label = G.degree(G.nodes(), weight=weight)
	k = [k_for_label[label_for_node[i]] for i in range(N)]

	# Initialize community and merge lists
	communities = {i: frozenset([i]) for i in range(N)}
	merges = []

	# Initial modularity
	partition = [[label_for_node[x] for x in c] for c in communities.values()]
	q_cnm = modularity(G, partition)

	# Initialize data structures
	# CNM Eq 8-9 (Eq 8 was missing a factor of 2 (from A_ij + A_ji)
	# a[i]: fraction of edges within community i
	# dq_dict[i][j]: dQ for merging community i, j
	# dq_heap[i][n] : (-dq, i, j) for communitiy i nth largest dQ
	# H[n]: (-dq, i, j) for community with nth largest max_j(dQ_ij)
	a = [k[i] * q0 for i in range(N)]
	dq_dict = {
		i: {
			j: 2 * q0 - 2 * k[i] * k[j] * q0 * q0
			for j in [node_for_label[u] for u in G.neighbors(label_for_node[i])]
			if j != i
		}
		for i in range(N)
	}
	dq_heap = [
		MappedQueue([(-dq, i, j) for j, dq in dq_dict[i].items()]) for i in range(N)
	]
	H = MappedQueue([dq_heap[i].h[0] for i in range(N) if len(dq_heap[i]) > 0])

	# Merge communities until we can't improve modularity
	while len(H) > 1:
		# Find best merge
		# Remove from heap of row maxes
		# Ties will be broken by choosing the pair with lowest min community id
		try:
			dq, i, j = H.pop()
		except IndexError:
			break
		dq = -dq
		# Remove best merge from row i heap
		dq_heap[i].pop()
		# Push new row max onto H
		if len(dq_heap[i]) > 0:
			H.push(dq_heap[i].h[0])
		# If this element was also at the root of row j, we need to remove the
		# duplicate entry from H
		if dq_heap[j].h[0] == (-dq, j, i):
			H.remove((-dq, j, i))
			# Remove best merge from row j heap
			dq_heap[j].remove((-dq, j, i))
			# Push new row max onto H
			if len(dq_heap[j]) > 0:
				H.push(dq_heap[j].h[0])
		else:
			# Duplicate wasn't in H, just remove from row j heap
			dq_heap[j].remove((-dq, j, i))
		# Stop when change is non-positive
		if dq <= 0:
			break

		# Perform merge
		communities[j] = frozenset(communities[i] | communities[j])
		del communities[i]
		merges.append((i, j, dq))
		# New modularity
		q_cnm += dq
		# Get list of communities connected to merged communities
		i_set = set(dq_dict[i].keys())
		j_set = set(dq_dict[j].keys())
		all_set = (i_set | j_set) - {i, j}
		both_set = i_set & j_set
		# Merge i into j and update dQ
		for k in all_set:
			# Calculate new dq value
			if k in both_set:
				dq_jk = dq_dict[j][k] + dq_dict[i][k]
			elif k in j_set:
				dq_jk = dq_dict[j][k] - 2.0 * a[i] * a[k]
			else:
				# k in i_set
				dq_jk = dq_dict[i][k] - 2.0 * a[j] * a[k]
			# Update rows j and k
			for row, col in [(j, k), (k, j)]:
				# Save old value for finding heap index
				if k in j_set:
					d_old = (-dq_dict[row][col], row, col)
				else:
					d_old = None
				# Update dict for j,k only (i is removed below)
				dq_dict[row][col] = dq_jk
				# Save old max of per-row heap
				if len(dq_heap[row]) > 0:
					d_oldmax = dq_heap[row].h[0]
				else:
					d_oldmax = None
				# Add/update heaps
				d = (-dq_jk, row, col)
				if d_old is None:
					# We're creating a new nonzero element, add to heap
					dq_heap[row].push(d)
				else:
					# Update existing element in per-row heap
					dq_heap[row].update(d_old, d)
				# Update heap of row maxes if necessary
				if d_oldmax is None:
					# No entries previously in this row, push new max
					H.push(d)
				else:
					# We've updated an entry in this row, has the max changed?
					if dq_heap[row].h[0] != d_oldmax:
						H.update(d_oldmax, dq_heap[row].h[0])

		# Remove row/col i from matrix
		i_neighbors = dq_dict[i].keys()
		for k in i_neighbors:
			# Remove from dict
			dq_old = dq_dict[k][i]
			del dq_dict[k][i]
			# Remove from heaps if we haven't already
			if k != j:
				# Remove both row and column
				for row, col in [(k, i), (i, k)]:
					# Check if replaced dq is row max
					d_old = (-dq_old, row, col)
					if dq_heap[row].h[0] == d_old:
						# Update per-row heap and heap of row maxes
						dq_heap[row].remove(d_old)
						H.remove(d_old)
						# Update row max
						if len(dq_heap[row]) > 0:
							H.push(dq_heap[row].h[0])
					else:
						# Only update per-row heap
						dq_heap[row].remove(d_old)

		del dq_dict[i]
		# Mark row i as deleted, but keep placeholder
		dq_heap[i] = MappedQueue()
		# Merge i into j and update a
		a[j] += a[i]
		a[i] = 0

	communities = [
		frozenset([label_for_node[i] for i in c]) for c in communities.values()
	]
	return sorted(communities, key=len, reverse=True)


def CNM_clustering(G, m, dist_intra = None, dist_inter = None):
	output = []
	for element in nx.connected_components(G):
		if len(element) > 1:
			G_sub = nx.Graph(G.subgraph(element))
			for u,v,info in G_sub.edges(data = True):
				if info["type"] == "intra":
					distance = dist_intra
				else:
					distance = dist_inter
				if info["distance"] > 0:
					G_sub[u][v]["weight"] = (distance / info["distance"]) ** m
				else:
					G_sub[u][v]["weight"] = np.inf
			# for thing in greedy_modularity_communities(format_graph(G_sub), weight = "weight"):
			for thing in greedy_modularity_communities(format_graph(G_sub)):
				G_cluster = G_sub.subgraph(set(thing))
				flag = 0
				for u,v,info in G_cluster.edges(data = True):
					if info["type"] == "inter":
						flag = 1
						break
				if flag == 1:
					output.append((set(thing), "inter"))
				else:
					output.append((set(thing), "intra"))
		else:
			output.append((element, "intra"))
	return output


def get_cluster_from_graph(graph_file, output, dist_intra, dist_inter, m):
	G = nx.read_graphml(graph_file)
	f = open(output, "w")
	if len(G.edges()) > 0:
		clusters = CNM_clustering(G, m, dist_intra = dist_intra, dist_inter = dist_inter)
		for x,y in clusters:
			f.write("\t".join([",".join(sorted(x)), y]) + "\n")
	else:
		for u in G.nodes():
			f.write("\t".join([u, "intra"]) + "\n")
	f.close()
	return output


def get_cluster_from_interface(uniprot2res, ires_file, binary_interactome, output):
	def get_mutated_cluster(row, uniprot2res):
		output1 = []
		for element in row["P1_IRES"].split(","):
			uniprot, res = element.split("_")
			if int(res) in uniprot2res[uniprot]:
				output1.append(element)
		output2 = []
		for element in row["P2_IRES"].split(","):
			uniprot, res = element.split("_")
			if int(res) in uniprot2res[uniprot]:
				output2.append(element)
		if len(output1) > 0 and len(output2) > 0:
			return ",".join(sorted(set(output1 + output2)))
		else:
			return np.nan

	# get interface residues
	df_pioneer = pd.read_csv(ires_file, sep = "\t", dtype = str).dropna(subset = ["P1_IRES", "P2_IRES"])
	df_pioneer = df_pioneer[df_pioneer["Interaction"].isin(binary_interactome)]

	# get mutated clusters
	if df_pioneer.shape[0] == 0:
		return None
	else:
		df_pioneer["Mutated_residues"] = df_pioneer.apply(lambda x: get_mutated_cluster(x, uniprot2res), axis = 1)
		df_pioneer = df_pioneer.dropna(subset = ["Mutated_residues"])
		if df_pioneer.shape[0] > 0:
			df_pioneer["Structure_source"] = "PIONEER"
			df_pioneer[["Mutated_residues", "Structure_source", "Interaction"]].to_csv(output, sep = "\t", header = None, index = None)
			return output
		else:
			return None


def get_cluster_pval(mutres_file, cluster_file, sample_size, output):
	def get_uniprot(res):
		uniprot = sorted(set([item.split("_")[0] for item in res.split(",")]))
		return ",".join(uniprot)
	def get_cluster_value(res, dictionary):
		values = []
		for element in res.split(","):
			values.append(str(dictionary[element]))
		return ",".join(values)
	def get_mutated_samples(res_mutatedsample):
		total_sample = []
		for element in res_mutatedsample.split(","):
			total_sample.extend(element.split("|"))
		return len(set(total_sample))
	def get_cluster_mutability(res_mutability):
		final_pval = 1
		for element in res_mutability.split(","):
			final_pval = final_pval * (1 - float(element))
		return 1 - final_pval

	# load mutated residues
	df_mut = pd.read_csv(mutres_file, sep = "\t", usecols = ["Uniprot", "Uniprot_Position", "Mutations", "res_mutability", "Ref_Tri", "Mutation_info"], \
		dtype = {"Uniprot_Position": int, "Mutations": str})
	
	# set a lower boundary for the mutability of residues in rarely mutated genes
	p_lowerlimit = df_mut["res_mutability"].quantile(q = 0.2)
	df_mut.loc[df_mut["res_mutability"] < p_lowerlimit, "res_mutability"] = p_lowerlimit

	# annotate background mutability and observed mutation count for the clusters
	df_mut["res"] = df_mut.apply(lambda row: row["Uniprot"] + "_" + str(row["Uniprot_Position"]), axis = 1)
	res2mutability = {}
	res2mutcount = {}
	res2trinuc = {}
	res2sample = {}
	res_mutability = df_mut["res_mutability"].tolist()
	res_mutcount = df_mut["Mutations"].tolist()
	res_trinuc = df_mut["Ref_Tri"].tolist()
	res_mutinfo = df_mut["Mutation_info"].tolist()
	for i, element in enumerate(df_mut["res"].tolist()):
		res2mutability[element] = res_mutability[i]
		res2mutcount[element] = res_mutcount[i]
		res2trinuc[element] = res_trinuc[i]
		res2sample[element] = "|".join([item.split(":")[0] for item in res_mutinfo[i].split(",")])
	df_cluster = pd.read_csv(cluster_file, sep = "\t", header = None).rename(columns = {0: "Residues", 1: "Structure_source", 2: "Uniprots"})
	if "Uniprots" not in df_cluster.columns:
		df_cluster["Uniprots"] = df_cluster["Residues"].apply(get_uniprot)
	df_cluster["res_mutability"] = df_cluster["Residues"].apply(lambda x: get_cluster_value(x, res2mutability))
	df_cluster["Mutation_count"] = df_cluster["Residues"].apply(lambda x: get_cluster_value(x, res2mutcount))
	df_cluster["res_trinuc"] = df_cluster["Residues"].apply(lambda x: get_cluster_value(x, res2trinuc))
	df_cluster["Mutated_sample"] = df_cluster["Residues"].apply(lambda x: get_cluster_value(x, res2sample))
	df_cluster["sample_mutcount"] = df_cluster["Mutated_sample"].apply(get_mutated_samples)
	df_cluster["cluster_mutability"] = df_cluster["res_mutability"].apply(get_cluster_mutability)
	df_cluster["Raw_pvalue"] = df_cluster.apply(lambda row: scipy.stats.binom_test(row["sample_mutcount"], sample_size, row["cluster_mutability"], alternative = "greater"), axis = 1)
	df_cluster.to_csv(output, sep = "\t", index = None, header = True)
	return output


def get_sig_cluster_intra(pval_file_PDB, pval_file_AlphaFold2, pdb_intra_resource, sig_cutoff = 0.05):
	def whether_retain(row):
		G = nx.read_graphml(pdb_intra_resource + row["Uniprots"] + ".graphml.gz")
		if set(row["Residues"].split(",")).intersection(set(G.nodes())) == set():
			return True
		else:
			return False

	if pval_file_PDB != None:
		df_pdb = pd.read_csv(pval_file_PDB, sep = "\t", dtype = {"Mutation_count": str})
		df_pdb["Adjusted_pvalue"] = df_pdb["Raw_pvalue"].apply(lambda x: x * df_pdb.shape[0])
		df_pdb = df_pdb[df_pdb["Adjusted_pvalue"] < sig_cutoff]
		pdb_uniprots = df_pdb["Uniprots"].tolist()
	else:
		df_pdb = pd.DataFrame()
		pdb_uniprots = []
	logging.info(str(df_pdb.shape[0]) + " significant intra-protein clusters are identfied using PDB structures.")
	
	if pval_file_AlphaFold2 != None:
		df_af2 = pd.read_csv(pval_file_AlphaFold2, sep = "\t", dtype = {"Mutation_count": str})
		df_af2["Adjusted_pvalue"] = df_af2["Raw_pvalue"].apply(lambda x: x * df_af2.shape[0])
		df_af2 = df_af2[df_af2["Adjusted_pvalue"] < sig_cutoff]
		df1 = df_af2[~df_af2["Uniprots"].isin(pdb_uniprots)]
		df2 = df_af2[df_af2["Uniprots"].isin(pdb_uniprots)]
		df2 = df2[df2.apply(whether_retain, axis = 1)]
		df_af2 = pd.concat([df1, df2])
	else:
		df_af2 = pd.DataFrame()
	logging.info(str(df_af2.shape[0]) + " significant intra-protein clusters are identfied using AlphaFold2 structures.")

	if pd.concat([df_pdb, df_af2]).shape[0] > 0:
		df_final = pd.concat([df_pdb, df_af2])[["Structure_source", "Uniprots", "Residues", "Mutated_sample", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
		df_final["Type"] = "NonTrunc_IntraProtein_Cluster"
	else:
		df_final = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutated_sample": [], "Mutation_count": [], \
			"Raw_pvalue": [], "Adjusted_pvalue": []})
	return df_final


def get_sig_cluster_inter(pval_file_PDB, pval_file_PIONEER, sig_cutoff = 0.05):
	if pval_file_PDB != None:
		df_pdb = pd.read_csv(pval_file_PDB, sep = "\t", dtype = {"Mutation_count": str})
		df_pdb["Adjusted_pvalue"] = df_pdb["Raw_pvalue"].apply(lambda x: x * df_pdb.shape[0])
		df_pdb = df_pdb[df_pdb["Adjusted_pvalue"] < sig_cutoff]
	else:
		df_pdb = pd.DataFrame()
	logging.info(str(df_pdb.shape[0]) + " significant inter-protein clusters are identfied using PDB structures.")
	
	if pval_file_PIONEER != None:
		df_pio = pd.read_csv(pval_file_PIONEER, sep = "\t", dtype = {"Mutation_count": str})
		df_pio["Adjusted_pvalue"] = df_pio["Raw_pvalue"].apply(lambda x: x * df_pio.shape[0]) 
		df_pio = df_pio[df_pio["Adjusted_pvalue"] < sig_cutoff]
	else:
		df_pio = pd.DataFrame()
	logging.info(str(df_pio.shape[0]) + " significant inter-protein clusters are identfied using PIONEER-predicted interfaces.")
	
	if pd.concat([df_pdb, df_pio]).shape[0] > 0:
		df_final = pd.concat([df_pdb, df_pio])[["Structure_source", "Uniprots", "Residues", "Mutated_sample", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
		df_final["Type"] = "NonTrunc_InterProtein_Cluster"
	else:
		df_final = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutated_sample": [], "Mutation_count": [], \
			"Raw_pvalue": [], "Adjusted_pvalue": []})
	return df_final


def binary_interaction(uniprots):
	if len(uniprots.split(",")) == 1:
		return uniprots + "," + uniprots
	else:
		return uniprots


def get_heat_diffusion_matrix(G, beta):
	'''
	Get the diffusion matrix of a network. 
	G: undirected graph (no isolated nodes)
	beta: restart probability / the fraction of retained heat 
	return:
	  --F: heat diffusion matrix
	  --all_nodes: a list of nodes corresponding to the index of heat diffusion matrix
	'''
	all_nodes = list(G.nodes())
	adjacency_matrix = nx.adjacency_matrix(G).todense()
	W = adjacency_matrix / np.sum(adjacency_matrix, axis = 0)
	F = beta * np.linalg.inv(np.identity(len(all_nodes)) - (1 - beta) * W)
	return [F, all_nodes]


def identify_hot_modules(G, graph_output, beta, delta):
	"""
	Parameters
	----------
	G: undirected graph. Each node has an attribute **heat_score**. heat_score = 0 means that there is no heat from that node.
	beta: restart probability
	delta: minimum amount of heat diffused from one node to another so that a directed edge will be added
	
	Return
	------
	a list of sets, the nodes in each set belong to the same hot module
	"""
	def get_hot_subnetworks(F, nodes, Dh, delta):
		"""
		F: heat diffusion matrix
		nodes: a list of nodes corresponding to the index of heat diffusion matrix
		Dh: diagonal matrix with the heat scores of **nodes**
		return: a list of sets, the nodes in each set belong to the same hot subnetwork (min size of a hot subnetwork is 2)
		"""
		idx2uniprot = {}
		for i in range(len(nodes)):
			idx2uniprot[i] = nodes[i]
		E = np.matmul(F, Dh)
		E[E <= delta] = 0
		G = nx.from_numpy_matrix(E.transpose(), create_using = nx.DiGraph)
		G = nx.relabel_nodes(G, idx2uniprot)
		return G

	all_modules = []
	G_final = nx.DiGraph()
	for element in nx.connected_components(G):
		heat_scores = [G.nodes[item]["heat_score"] for item in element]
		if max(heat_scores) > 0:
			if len(element) > 1:
				G_sub = G.subgraph(element)
				F, nodes = get_heat_diffusion_matrix(G_sub, beta = beta)
				diag_elements = [G_sub.nodes[item]["heat_score"] for item in nodes]
				G_one = get_hot_subnetworks(F, nodes, np.diag(diag_elements), delta = delta)
				G_final.add_nodes_from(G_one.nodes(data = True))
				G_final.add_edges_from(G_one.edges(data = True))
				for x in nx.strongly_connected_components(G_one):
					if len(x) > 1 or G.nodes[list(x)[0]]["heat_score"] > 0:
						all_modules.append(x)
			else:
				all_modules.append(element)
	nx.write_graphml(G_final, graph_output)
	return all_modules


def get_initial_distribution(G_original, pairs, ones, lof, heat_source_score, edge_weight, output):
	"""
	Purpose
	----------
	generate an undirected graph whose nodes have attribute **heat_score** and whose edges have attribute **weight**
	
	
	Parameters
	----------
	G: PPI network. Node IDs are UniProt IDs.
	pairs: list of protein pairs with hot clusters at the interaction interface, e.g. [("P12I92", "Q13924"), ("P01116", "O29389"), ...]
	ones: list of proteins with intra-protein hot clusters, e.g. ["Q03202", "P02933", ...]
	lof: list of proteins enriched of LoF mutations, e.g. ["O23802", "P29031", ...]
	"""
	# initialization
	G = copy.deepcopy(G_original)
	for u,v in G.edges():
		G[u][v]["weight"] = 1.0
	for u in G.nodes():
		G.nodes[u]["heat_score"] = 0.0

	# set enhanced edge conductivity and assign heat sources
	for pro1, pro2 in pairs:
		G[pro1][pro2]["weight"] = edge_weight
		G.nodes[pro1]["heat_score"] = heat_source_score
		G.nodes[pro2]["heat_score"] = heat_source_score
	for pro in ones:
		if pro in G:
			G.nodes[pro]["heat_score"] = heat_source_score
	for pro in lof:
		if pro in G:
			G.nodes[pro]["heat_score"] = G.nodes[pro]["heat_score"] + heat_source_score
	G = format_graph(G)
	nx.write_graphml(G, output)
	return G


def get_suitable_delta(G, beta, size_cutoff, output, upper_bound = 1, lower_bound = 0):
	"""
	Purpose
	--------- 
	- To find the appropriate delta for a network with labeled heat sources
	- delta: minimum amount of heat diffused from one node to another so that a directed edge will be added
	- "Appropriate" means that using this delta value there won't be giant modules. The maximum module size would be **size_cutoff**
	
	Parameters
	----------
	G: undirected graph. Each node has an attribute **heat_score**. heat_score = 0 means that there is no heat placed on that node.
	beta: restart probability / the fraction of retained heat
	size_cutoff: the maximum size of a hot module
	

	Returns
	-------
	final_delta: the appropriate delta value
	final_size: the list of hot module sizes using **final_delta** as cutoff
	"""
	# get edge values 
	heat_edges = []
	E_list = []
	for element in nx.connected_components(G):
		heat_scores = [G.nodes[item]["heat_score"] for item in element]
		if len(element) > 1 and max(heat_scores) > 0:
			G_sub = G.subgraph(element)
			F, nodes = get_heat_diffusion_matrix(G_sub, beta = beta)
			diag_elements = [G_sub.nodes[item]["heat_score"] for item in nodes]
			E = np.matmul(F, np.diag(diag_elements))
			E_list.append(E)
			heat_edges.extend(E[E.nonzero()].tolist()[0])
	heat_edges = sorted(heat_edges, reverse = True)

	# get the maximum difference between upper_bound and lower_bound
	diff = max(2 / len(heat_edges), 0.00001)

	# find appropriate delta
	module_size = []
	for E_copy in copy.deepcopy(E_list):
		G = nx.from_numpy_matrix(E_copy.transpose(), create_using = nx.DiGraph)
		for element in list(nx.strongly_connected_components(G)):
			module_size.append(len(element))
	if max(module_size) <= size_cutoff:
		final_delta = -1
		final_midpoint = 1
	else:
		while upper_bound - lower_bound > diff:
			midpoint = (upper_bound + lower_bound) / 2
			delta = heat_edges[int(midpoint * len(heat_edges))]
			module_size = []
			for E_copy in copy.deepcopy(E_list):
				E_copy[E_copy <= delta] = 0
				G = nx.from_numpy_matrix(E_copy.transpose(), create_using = nx.DiGraph)
				for element in list(nx.strongly_connected_components(G)):
					module_size.append(len(element))
			if max(module_size) <= size_cutoff:
				lower_bound = midpoint
			else:
				upper_bound = midpoint
		final_midpoint = lower_bound
		final_delta = heat_edges[int(final_midpoint * len(heat_edges)) + 1]

	# write output
	f = open(output, "a+")
	f.write("\t".join([str(final_delta), str(final_midpoint), ",".join([str(thing) for thing in module_size])]) + "\n")
	f.close()
	
	return


def get_one_delta(G_original, delta_output, restart_prob, max_subnetwork_size, seed):
	G = copy.deepcopy(G_original)
	nx.double_edge_swap(G, nswap = len(G.edges()), max_tries = 200000, seed = seed)
	for u,v in G.edges():
		G[u][v]["weight"] = 1.0
	get_suitable_delta(G, beta = restart_prob, size_cutoff = max_subnetwork_size, output = delta_output)
	return


def remove_whole_dir(dir_name):
	shutil.rmtree(dir_name)
	return


