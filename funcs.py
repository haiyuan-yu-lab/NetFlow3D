description='''
- Functions used by NetFlow3D.py
- Dependencies: None
'''

import scipy.stats
import numpy as np
import os
from collections import defaultdict
import pandas as pd
import networkx as nx
import random
import copy
from os import path
import logging
import shutil
from scipy.stats import poisson


def whether_canonical(uniprots, canonical_isoforms):
	"""
	Check if a list of comma-separated UniProt IDs are the canonical isoforms of their encoding gene
	"""
	# Initialize an empty list to store the results
	output = []
	
	# Split the input string of UniProt IDs by comma
	for x in uniprots.split(","):
		# Check if the current UniProt ID is in the list of canonical isoforms
		if x in canonical_isoforms:
			# If it is, append "Yes" to the output list
			output.append("Yes")
		else:
			# If it is not, append "No" to the output list
			output.append("No")
	
	# Join the output list into a comma-separated string and return it
	return ",".join(output)


def extract_uniprot_ids(fasta_file):
	"""
	Extracts UniProt IDs from a given FASTA file.
	"""
	# Initialize an empty list to store UniProt IDs
	uniprot_ids = []
	
	# Open the given FASTA file in read mode
	with open(fasta_file, 'r') as file:
		# Iterate through each line in the file
		for line in file:
			# Check if the line starts with '>', indicating a header line
			if line.startswith('>'):
				# Split the header line by '|' and extract the UniProt ID (second element)
				parts = line.split('|')
				uniprot_id = parts[1]
				# Append the extracted UniProt ID to the list
				uniprot_ids.append(uniprot_id)
	
	# Return the list of extracted UniProt IDs
	return uniprot_ids


def split_consecutive_pos(pos):
	"""
	Splits a string representing a range of positions into a list of integers.
	"""
	try:
		# Check if the position is a single hyphen, indicating no positions
		if pos == "-":
			return []
		# Check if the position is a single number, without any hyphen
		elif pos.find('-') == -1:
			return [int(pos)]
		else:
			# Split the position range into start (init) and end (end)
			init, end = pos.split('-')
			# Handle case where the start position is unknown
			if init == '?':
				return [int(end)]
			# Handle case where the end position is unknown
			if end == '?':
				return [int(init)]
			# Return the list of positions in the specified range
			return list(range(int(init), int(end) + 1))
	except:
		# Log a warning if an error occurs and return an empty list
		logging.warning("Unknown protein position: " + pos)
		return []


def get_res_annotation(df, output):
	"""
	Gathers mutation information for each residue from a DataFrame and writes it to an output file.
	"""
	def gather_mutation_info(row):
		"""
		Processes each row to gather mutation information and update dictionaries.
		"""
		for element in split_consecutive_pos(row["Protein_position"]):
			# Update mutation count for the (UniProt, position) key
			res2mutcount[(row["UniProt"], element)] += 1
			# Append mutation information to the corresponding key
			res2mutinfo[(row["UniProt"], element)].append("_".join([row["Chromosome"], row["Start_Position"], row["Tumor_Sample_Barcode"]]))
		return True

	# Initialize dictionaries to store mutation counts and information
	res2mutcount = defaultdict(int)
	res2mutinfo = defaultdict(list)
	
	# Apply the gather_mutation_info function to each row in the DataFrame
	df["status"] = df.apply(gather_mutation_info, axis=1)

	# Open the output file in write mode
	with open(output, "w") as f:
		# Write the header line to the output file
		f.write("\t".join(["UniProt", "UniProt_Position", "Mutations", "Mutation_info"]) + "\n")
		# Write mutation data for each key (i.e. residue) in the dictionaries
		for key in res2mutcount:
			f.write("\t".join([key[0], str(key[1]), str(res2mutcount[key]), ",".join(sorted(set(res2mutinfo[key])))]) + "\n")
	return


def get_graph_by_uniprot_binary(uniprot, uniprot2res, output_path, resource_dir_intra, resource_dir_inter=None, dist_intra=np.inf, dist_inter=np.inf):
	"""
	Given all the mutated residues (MR) in a protein, generate intra- and inter-chain MR-MR contact network.

	Args:
		uniprot (str): UniProt ID of the protein.
		uniprot2res (dict): Dictionary mapping UniProt IDs to lists of mutated residues.
		output_path (str): Path to save the output graph files.
		resource_dir_intra (str): Directory containing intra-chain residue-residue distance maps.
		resource_dir_inter (str, optional): Directory containing inter-chain residue-residue distance maps. Defaults to None.
		dist_intra (float, optional): Maximum distance for intra-chain interactions. Defaults to infinity.
		dist_inter (float, optional): Maximum distance for inter-chain interactions. Defaults to infinity.

	Returns:
		list: List of file paths to the generated graph files.
	"""
	output_list = []

	# Check if intra-chain residue-residue distance map exists
	if path.exists(resource_dir_intra + uniprot + ".graphml.gz"):
		# Load the intra-chain residue-residue distance map
		G = nx.read_graphml(resource_dir_intra + uniprot + ".graphml.gz")
		
		# Get residues with mutations that exist in the intra-chain graph
		res_list = []
		for item in uniprot2res[uniprot]:
			if uniprot + "_" + str(item) in G:
				res_list.append(uniprot + "_" + str(item))
		
		# Create a subgraph induced by the residues with mutations
		if len(res_list) > 0:
			G_sub = nx.Graph(G.subgraph(res_list))
			for element in list(G_sub.edges()):
				# Remove edges with distances greater than dist_intra
				if G_sub[element[0]][element[1]]['distance'] > dist_intra:
					G_sub.remove_edge(element[0], element[1])
				else:
					G_sub[element[0]][element[1]]['type'] = 'intra'
			# Save the intra-chain graph
			nx.write_graphml(G_sub, output_path + uniprot + ".graphml")
			output_list.append(output_path + uniprot + ".graphml")

		# If inter-chain directory is provided and file exists
		if resource_dir_inter is not None and path.exists(resource_dir_inter + uniprot + ".graphml.gz"):
			# Load the inter-chain residue-residue distance map
			G_inter = nx.read_graphml(resource_dir_inter + uniprot + ".graphml.gz")
			
			# Get residues with mutations that exist in the inter-chain graph
			res_list_inter = []
			for element in G_inter.nodes():
				element_split = element.split("_")
				if int(element_split[1]) in uniprot2res[element_split[0]]:
					res_list_inter.append("_".join(element_split))
			
			# Create a subgraph induced by the residues with mutations
			if len(res_list_inter) > 0:
				G_inter_sub = nx.Graph(G_inter.subgraph(res_list_inter))
				for element in list(G_inter_sub.edges()):
					# Remove edges with distances greater than dist_inter
					if G_inter_sub[element[0]][element[1]]['distance'] > dist_inter:
						G_inter_sub.remove_edge(element[0], element[1])
					else:
						G_inter_sub[element[0]][element[1]]['type'] = 'inter'

				# Get all interactors of the input protein
				interactors = set([item.split("_")[0] for item in res_list_inter])
				for interactor in interactors:
					interactor_file = output_path + "_".join(sorted([interactor, uniprot])) + ".graphml"
					if not path.exists(interactor_file):
						# Get all edges connecting the input uniprot and this interactor
						edge_list = [item for item in G_inter_sub.edges() if {item[0].split("_")[0], item[1].split("_")[0]} == {uniprot, interactor}]
						
						if edge_list:
							# Create a subgraph induced by those inter-chain edges
							G_tmp = G_inter_sub.edge_subgraph(edge_list)
							G_final = nx.Graph()
							G_final.add_edges_from(G_tmp.edges(data=True))
							G_final.add_nodes_from(G_tmp.nodes(data=True))

							# Add intra-chain edges in the input uniprot
							if res_list:
								G_final.add_edges_from(G_sub.edges(data=True))
								G_final.add_nodes_from(G_sub.nodes(data=True))

							# Add intra-chain edges in the interactor if it's a heterodimer
							if interactor != uniprot and path.exists(resource_dir_intra + interactor + ".graphml.gz"):
								G_interactor = nx.read_graphml(resource_dir_intra + interactor + ".graphml.gz")
								res_list_interactor = [interactor + "_" + str(item) for item in uniprot2res[interactor] if interactor + "_" + str(item) in G_interactor]
								
								if res_list_interactor:
									G_interactor_sub = nx.Graph(G_interactor.subgraph(res_list_interactor))
									for item in list(G_interactor_sub.edges()):
										if G_interactor_sub[item[0]][item[1]]["distance"] > dist_intra:
											G_interactor_sub.remove_edge(item[0], item[1])
										else:
											G_interactor_sub[item[0]][item[1]]["type"] = "intra"
									G_final.add_edges_from(G_interactor_sub.edges(data=True))
									G_final.add_nodes_from(G_interactor_sub.nodes(data=True))

							# Save the combined graph
							nx.write_graphml(G_final, interactor_file)
							output_list.append(interactor_file)
	
	return output_list


def get_cluster_from_graph(graph_file, output):
	"""
	Extracts clusters (i.e. connected components) from a graph and writes them to an output file.
	"""
	# Read the graph from the GraphML file
	G = nx.read_graphml(graph_file)
	
	# Open the output file in write mode
	with open(output, "w") as f:
		# Iterate through each connected component in the graph
		for element in nx.connected_components(G):
			# Create a subgraph induced by the connected component
			G_sub = G.subgraph(element)
			
			# Assume the cluster is intra-chain unless an inter-chain edge is found
			flag = "intra"
			
			# Check edges in the subgraph to determine if it contains any inter-chain connections
			for u, v, info in G_sub.edges(data=True):
				if info["type"] == "inter":
					flag = "inter"
					break
			
			# Write the cluster and its type to the output file
			f.write("\t".join([",".join(sorted(element)), flag]) + "\n")
	
	# Return the path to the output file
	return output


def get_cluster_from_interface(uniprot2res, ires_file, binary_interactome, output):
	"""
	Generates mutated clusters from interface residues and writes them to an output file.
	"""
	
	def get_mutated_cluster(row, uniprot2res):
		"""
		Identifies mutated residues on the PPI interface for a given row.

		Args:
			row (pd.Series): A row of the DataFrame.
			uniprot2res (dict): Dictionary mapping UniProt IDs to lists of mutated residues.

		Returns:
			str: Comma-separated string of mutated residues on the interface.
		"""
		output1 = []
		# Check and process interface residues on P1
		if type(row["P1_IRES"]) == str:
			for element in row["P1_IRES"].split(","):
				uniprot, res = element.split("_")
				if int(res) in uniprot2res[uniprot]:
					output1.append(element)
		
		output2 = []
		# Check and process interface residues on P2
		if type(row["P2_IRES"]) == str:
			for element in row["P2_IRES"].split(","):
				uniprot, res = element.split("_")
				if int(res) in uniprot2res[uniprot]:
					output2.append(element)
		
		# Combine and sort the mutated residues on P1-P2 interfaces
		return ",".join(sorted(set(output1 + output2)))

	# Load the interface residues file, filter out rows with NaNs in both P1_IRES and P2_IRES
	df_pioneer = pd.read_csv(ires_file, sep="\t", dtype=str).dropna(subset=["P1_IRES", "P2_IRES"], how="all")
	
	# Filter interactions based on the binary interactome
	df_pioneer = df_pioneer[df_pioneer["Interaction"].isin(binary_interactome)].rename(columns={"Interaction": "Uniprots"})
	
	# Apply the get_mutated_cluster function to each row to get mutated clusters
	df_pioneer["Residues"] = df_pioneer.apply(lambda x: get_mutated_cluster(x, uniprot2res), axis=1)
	
	# Filter out rows where no mutated residues were found
	df_pioneer = df_pioneer[df_pioneer["Residues"] != ""]
	
	# Add a column indicating the structure source
	df_pioneer["Structure_source"] = "PIONEER"
	
	# Write the resulting DataFrame to the output file
	df_pioneer[["Residues", "Structure_source", "Uniprots"]].to_csv(output, sep="\t", header=True, index=None)

	return


def get_cluster_pval(mutres_file, cluster_file, sample_size, uniprot2inframe_mutrate):
	"""
	Annotates clusters with observed mutation counts, expected mutation counts, and p-values.

	Args:
		mutres_file (str): Path to the input file containing mutated residues.
		cluster_file (str): Path to the input file containing cluster information.
		sample_size (int): The sample size to be used in mutation rate calculations.
		uniprot2inframe_mutrate (dict): Dictionary mapping UniProt IDs to in-frame mutation rates.
	"""
	
	def get_uniprot(res):
		"""
		Extracts and returns a comma-separated string of unique UniProt IDs from residue identifiers.

		Args:
			res (str): Comma-separated residue identifiers.

		Returns:
			str: Comma-separated UniProt IDs.
		"""
		uniprot = sorted(set([item.split("_")[0] for item in res.split(",")]))
		return ",".join(uniprot)

	def get_cluster_value(res, dictionary):
		"""
		Retrieves mutation information for residues in a cluster from a dictionary.

		Args:
			res (str): Comma-separated residue identifiers.
			dictionary (dict): Dictionary mapping residue identifiers to mutation information.

		Returns:
			str: Comma-separated mutation information of residues in a cluster.
		"""
		if isinstance(res, str):
			values = []
			for element in res.split(","):
				values.append(str(dictionary[element]))
			return ",".join(values)
		else:
			return np.nan

	def get_total_mut(res, dictionary):
		"""
		Calculates total number of mutations in a cluster.

		Args:
			res (str): Comma-separated residue identifiers.
			dictionary (dict): Dictionary mapping residue identifiers to mutation information.

		Returns:
			int: Total number of mutations.
		"""
		if isinstance(res, str):
			output = []
			for element in res.split(","):
				output.extend(dictionary[element])
			return len(set(output))
		else:
			return 0

	def get_expected_mutcount(res):
		"""
		Calculate the expected mutation count for a cluster.
		"""
		# Generate a dictionary mapping UniProt identifiers to mutated residues in a cluster
		uniprot2res = defaultdict(list)
		if isinstance(res, str):
			for x in res.split(","):
				uniprot, pos = x.split("_")
				uniprot2res[uniprot].append(pos)

		# Calculate expected mutation count
		expected = 0
		# Iterate over each UniProt ID in the dictionary
		for x in uniprot2res:
			# Expected count = per base pair mutation rate * number of mutated residues * 3 (for codon) * sample size
			expected += uniprot2inframe_mutrate[x] * len(uniprot2res[x]) * 3 * sample_size

		return expected

	# Load mutated residues file
	df_mut = pd.read_csv(mutres_file, sep="\t", dtype={"UniProt_Position": int, "Mutations": str})
	df_mut["res"] = df_mut.apply(lambda row: row["UniProt"] + "_" + str(row["UniProt_Position"]), axis=1)
	
	# Create dictionaries to store mutation counts and mutation info for each residue
	res2mutcount = {}
	res2mutinfo = {}
	res_mutcount = df_mut["Mutations"].tolist()
	res_mutinfo = df_mut["Mutation_info"].tolist()
	for i, element in enumerate(df_mut["res"].tolist()):
		res2mutcount[element] = res_mutcount[i]
		res2mutinfo[element] = res_mutinfo[i].split(",")

	# Load cluster file
	df_cluster = pd.read_csv(cluster_file, sep="\t")
	if "Uniprots" not in df_cluster.columns:
		df_cluster["Uniprots"] = df_cluster["Residues"].apply(get_uniprot)
	
	# Calculate observed mutation count, expected mutation count, and per-residue mutation information for each cluster
	df_cluster["Observed_mut"] = df_cluster["Residues"].apply(lambda x: get_total_mut(x, res2mutinfo))
	df_cluster["Expected_mut"] = df_cluster["Residues"].apply(get_expected_mutcount)
	df_cluster["Mutation_count"] = df_cluster["Residues"].apply(lambda x: get_cluster_value(x, res2mutcount))
	
	# Calculate raw p-values for each cluster using the Poisson distribution
	df_cluster["Raw_pvalue"] = df_cluster.apply(lambda row: poisson.sf(row["Observed_mut"] - 1, row["Expected_mut"]), axis=1)
	
	# Save the annotated cluster information to the output file
	df_cluster.to_csv(cluster_file, sep="\t", index=None, header=True)
	
	return


def get_sig_cluster_structure(pval_file_PDB_intra, pval_file_PDB_inter, pval_file_AlphaFold2):
	"""
	Process and combine 3D clusters identified using atomic-resolution 3D structures.
	"""
	
	def whether_retain(res, uniprot2res):
		"""
		Checks if any residue in a cluster identified using AlphaFold2 structure is already present in clusters identified using PDB structure. If so, exclude this cluster.
		"""
		for x in res.split(","):
			uniprot, pos = x.split("_")
			if pos in uniprot2res[uniprot]:
				return False
		return True

	# Initialize a dictionary to store residues in the clusters identified using PDB structure.
	uniprot2res = defaultdict(list)

	# Process PDB intra-chain clusters if the file is provided
	if pval_file_PDB_intra is not None:
		df_pdb_intra = pd.read_csv(pval_file_PDB_intra, sep="\t", dtype={"Mutation_count": str})
		df_pdb_intra["Type"] = "InFrame_IntraProtein_Cluster"
		for x in df_pdb_intra["Residues"]:
			for y in x.split(","):
				uniprot, pos = y.split("_")
				uniprot2res[uniprot].append(pos)
	else:
		df_pdb_intra = pd.DataFrame()

	# Process AlphaFold2 intra-chain clusters if the file is provided
	if pval_file_AlphaFold2 is not None:
		df_af2 = pd.read_csv(pval_file_AlphaFold2, sep="\t", dtype={"Mutation_count": str})
		# Exclude AlphaFold2 intra-chain clusters that overlap with PDB intra-chain clusters
		df1 = df_af2[~df_af2["Uniprots"].isin(uniprot2res)]
		df2 = df_af2[df_af2["Uniprots"].isin(uniprot2res)]
		df2 = df2[df2["Residues"].apply(lambda x: whether_retain(x, uniprot2res))]
		df_af2 = pd.concat([df1, df2])
		df_af2["Type"] = "InFrame_IntraProtein_Cluster"
	else:
		df_af2 = pd.DataFrame()

	# Process PDB inter-chain clusters if the file is provided
	if pval_file_PDB_inter is not None:
		df_pdb_inter = pd.read_csv(pval_file_PDB_inter, sep="\t", dtype={"Mutation_count": str})
		df_pdb_inter["Type"] = "InFrame_InterProtein_Cluster"
	else:
		df_pdb_inter = pd.DataFrame()

	# Combine all the dataframes
	df_final = pd.concat([df_pdb_intra, df_af2, df_pdb_inter])

	# Adjust p-values using Bonferroni correction
	if df_final.shape[0] > 0:
		df_final["Adjusted_pvalue"] = df_final["Raw_pvalue"].apply(lambda x: min(x * df_final.shape[0], 1))
		df_final = df_final[["Type", "Structure_source", "Uniprots", "Residues", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
	else:
		df_final = pd.DataFrame({x: [] for x in ["Type", "Structure_source", "Uniprots", "Residues", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]})

	return df_final
	

def binary_interaction(uniprots):
	"""
	For homodimers, format the involved UniProt identifier as uniprot-uniprot.
	"""
	# Split the input string by commas and check the length of the resulting list
	if len(uniprots.split(",")) == 1:
		# If there is only one UniProt identifier, format it as a homodimer (uniprot,uniprot)
		return uniprots + "," + uniprots
	else:
		# If there are two UniProt identifiers, indicating a heterodimer, return the original string
		return uniprots


def get_heat_diffusion_matrix(G, beta):
	'''
	Get the diffusion matrix of a network.
	
	Args:
		G (networkx.Graph): An undirected graph with no isolated nodes.
		beta (float): Restart probability or the fraction of retained heat.
	
	Returns:
		tuple: A tuple containing the heat diffusion matrix and a list of nodes.
			- F (numpy.ndarray): The heat diffusion matrix.
			- all_nodes (list): A list of nodes corresponding to the index of the heat diffusion matrix.
	'''
	
	# Get the list of all nodes in the graph
	all_nodes = list(G.nodes())
	
	# Get the adjacency matrix of the graph and convert it to a dense format
	adjacency_matrix = nx.adjacency_matrix(G).todense()
	
	# Normalize the adjacency matrix by the sum of each column
	W = adjacency_matrix / np.sum(adjacency_matrix, axis=0)
	
	# Calculate the heat diffusion matrix: F = beta * (I - (1 - beta) * W)^-1
	F = beta * np.linalg.inv(np.identity(len(all_nodes)) - (1 - beta) * W)
	
	return [F, all_nodes]


def identify_hot_modules(G, graph_output, beta, delta):
	"""
	Identifies interconnected modules in a graph based on heat diffusion.

	Parameters
	----------
	G: networkx.Graph
		Undirected graph where each node has an attribute 'heat_score'.
		A heat_score of 0 means there is no heat from that node.
	beta: float
		Restart probability.
	delta: float
		Minimum amount of heat diffused from one node to another to add a directed edge.

	Returns
	-------
	list of sets
		Each set contains nodes belonging to the same interconnected module.
	"""
	
	def get_hot_subnetworks(F, nodes, Dh, delta):
		"""
		Generates a directed graph representing the heat diffusion network.

		Parameters
		----------
		F: numpy.ndarray
			Diffusion matrix of a network.
		nodes: list
			List of nodes corresponding to the index of the heat diffusion matrix.
		Dh: numpy.ndarray
			Diagonal matrix with the heat scores of nodes.
		
		Returns
		-------
		networkx.DiGraph
			Directed graph representing the heat diffusion network.
		"""
		# Create a mapping from index to UniProt ID (node)
		idx2uniprot = {i: nodes[i] for i in range(len(nodes))}
		
		# Multiply diffusion matrix with the diagonal heat score matrix to obtain the heat diffusion matrix
		E = np.matmul(F, Dh)
		
		# Zero out values less than or equal to delta to filter insignificant heat transfers
		E[E <= delta] = 0
		
		# Construct a directed graph from the filtered heat diffusion matrix
		# G = nx.from_numpy_matrix(E.transpose(), create_using=nx.DiGraph)
		G = nx.from_numpy_array(E.transpose(), create_using=nx.DiGraph)
		
		# Relabel the graph nodes with UniProt IDs
		G = nx.relabel_nodes(G, idx2uniprot)
		
		return G

	all_modules = []
	G_final = nx.DiGraph()
	
	# Iterate through each connected component in the graph
	for element in nx.connected_components(G):
		# Get the heat scores for the nodes in the component
		heat_scores = [G.nodes[item]["heat_score"] for item in element]
		
		# Process the component if it contains any heat
		if max(heat_scores) > 0:
			if len(element) > 1:
				# Create a subgraph for the connected component
				G_sub = G.subgraph(element)
				
				# Get the diffusion matrix and corresponding nodes of the connected component
				F, nodes = get_heat_diffusion_matrix(G_sub, beta=beta)
				
				# Create a diagonal matrix of heat scores
				diag_elements = [G_sub.nodes[item]["heat_score"] for item in nodes]
				Dh = np.diag(diag_elements)
				
				# Get the heat diffusion network for the connected component
				G_one = get_hot_subnetworks(F, nodes, Dh, delta=delta)
				
				# Add nodes and edges of the heat diffusion subnetwork to the final graph
				G_final.add_nodes_from(G_one.nodes(data=True))
				G_final.add_edges_from(G_one.edges(data=True))
				
				# Identify strongly connected components in the heat diffusion network
				for x in nx.strongly_connected_components(G_one):
					if len(x) > 1 or G.nodes[list(x)[0]]["heat_score"] > 0:
						all_modules.append(x)
			else:
				all_modules.append(element)
	
	# Write the final graph to an output file
	nx.write_graphml(G_final, graph_output)
	
	return all_modules


def format_graph(G):
	"""
	Sorts the nodes and edges in a graph.

	Parameters
	----------
	G : networkx.Graph
		The input graph to be formatted.
	
	Returns
	-------
	networkx.Graph
		A new graph with sorted nodes and edges.
	"""
	
	# Create a new empty graph
	H = nx.Graph()
	
	# Add nodes to the new graph, sorted by node identifiers and preserving node attributes
	H.add_nodes_from(sorted(G.nodes(data=True)))
	
	# Sort edges by node identifiers and add them to the new graph, preserving edge attributes
	H.add_edges_from(sorted([(sorted([u, v])[0], sorted([u, v])[1], info) for u, v, info in G.edges(data=True)]))
	
	return H


def get_initial_distribution(G_original, final_output_intra_pdb, final_output_intra_af2, final_output_inter_pdb, final_output_inter_pioneer, final_output_intra_lof, output, intercept):
	"""
	Initializes the distribution of heat scores and edge weights in a graph based on p-values from various data sources.

	Parameters
	----------
	G_original : networkx.Graph
		The original undirected graph with nodes and edges.
	final_output_intra_pdb : str
		Path to the file containing intra-chain cluster p-values from PDB.
	final_output_intra_af2 : str
		Path to the file containing intra-chain cluster p-values from AlphaFold2.
	final_output_inter_pdb : str
		Path to the file containing inter-chain cluster p-values from PDB.
	final_output_inter_pioneer : str
		Path to the file containing inter-chain cluster p-values from PIONEER.
	final_output_intra_lof : str
		Path to the file containing loss-of-function p-values.
	output : str
		Path to the output file where the initialized graph will be saved.
	intercept : float
		Baseline edge weight to be set for all edges.

	Returns
	-------
	networkx.Graph
		The initialized graph with updated heat scores and edge weights.
	"""
	
	def get_uniprot2pval(row):
		"""
		Populates the uniprot2pval dictionary with p-values for each UniProt ID.
		"""
		if isinstance(row["Residues"], str):
			for uniprot in {x.split("_")[0] for x in row["Residues"].split(",")}:
				uniprot2pval[uniprot].append(row["Raw_pvalue"])
		return True

	def get_interaction2pval(row):
		"""
		Populates the interaction2pval dictionary with p-values for each interaction.
		"""
		interaction2pval[tuple(sorted(row["Uniprots"].split(",")))].append(row["Raw_pvalue"])
		return True

	# Initialization: Deep copy the original graph
	G = copy.deepcopy(G_original)
	
	# Set initial edge weights to the baseline value
	for u, v in G.edges():
		G[u][v]["weight"] = intercept
	
	# Set initial node heat scores to 0.0
	for u in G.nodes():
		G.nodes[u]["heat_score"] = 0.0

	# Update edge weights based on interaction p-values
	interaction2pval = defaultdict(list)
	for file in [final_output_inter_pdb, final_output_inter_pioneer]:
		if file is not None:
			df = pd.read_csv(file, sep="\t")
			df["status"] = df.apply(get_interaction2pval, axis=1)
	for u, v in interaction2pval:
		if G.has_edge(u, v):
			# Adjust edge weight by adding the negative log of the minimum p-value, capped at 300
			G[u][v]["weight"] += min(-np.log10(min(interaction2pval[(u, v)])), 300)

	# Set node weights based on UniProt p-values
	uniprot2pval = defaultdict(list)
	for file in [final_output_intra_pdb, final_output_intra_af2, final_output_inter_pdb, final_output_inter_pioneer]:
		if file is not None:
			df = pd.read_csv(file, sep="\t")
			df["status"] = df.apply(get_uniprot2pval, axis=1)
	for uniprot in set(uniprot2pval).intersection(set(G.nodes())):
		# Adjust node heat score by adding the negative log of the minimum p-value, capped at 300
		G.nodes[uniprot]["heat_score"] = min(-np.log10(min(uniprot2pval[uniprot])), 300)
	
	if final_output_intra_lof is not None:
		df = pd.read_csv(final_output_intra_lof, sep="\t").set_index("Uniprots")
		for uniprot in set(df.index).intersection(set(G.nodes())):
			# Add the negative log of the p-value from loss-of-function data to the node heat score
			G.nodes[uniprot]["heat_score"] += min(-np.log10(df.loc[uniprot, "Raw_pvalue"]), 300)

	# Write the modified graph to the output file
	G = format_graph(G)
	nx.write_graphml(G, output)
	
	return G


def get_suitable_delta(G, beta, size_cutoff, output, upper_bound=1, lower_bound=0):
	"""
	Purpose:
	---------
	- To determine the optimal delta value for an edge weight cutoff in a heat diffusion network.
	- delta: the minimum amount of heat diffused from one node to another to add a directed edge.
	- "Optimal" implies that using this delta value will prevent the formation of overly large interconnected modules, with the maximum module size limited to **size_cutoff**.

	Parameters
	----------
	G : networkx.Graph
		Undirected graph where each node has a **heat_score** attribute. A heat_score of 0 indicates no heat placed on that node.
	beta : float
		Restart probability or the fraction of retained heat.
	size_cutoff : int
		The maximum permissible size for an interconnected module in a heat diffusion network derived from a randomly shuffled graph.
	output : str
		Path to the output file where the delta values and module sizes will be saved.
	upper_bound : float, optional
		Upper bound for the delta search (default is 1).
	lower_bound : float, optional
		Lower bound for the delta search (default is 0).

	Returns
	-------
	final_delta : float
		The optimal delta value.
	final_size : list
		The list of interconnected module sizes using **final_delta** as the cutoff.
	"""
	
	# Calculate heat diffusion for each connected component and collect edge values
	heat_edges = []
	E_list = []
	for element in nx.connected_components(G):
		heat_scores = [G.nodes[item]["heat_score"] for item in element]
		if len(element) > 1 and max(heat_scores) > 0:
			# Extract the subgraph for the current connected component
			G_sub = G.subgraph(element)
			# Compute the diffusion matrix for the subgraph
			F, nodes = get_heat_diffusion_matrix(G_sub, beta=beta)
			# Create a diagonal matrix of heat scores for the nodes
			diag_elements = [G_sub.nodes[item]["heat_score"] for item in nodes]
			E = np.matmul(F, np.diag(diag_elements))
			# Store the heat diffusion matrix and non-zero edge values
			E_list.append(E)
			# heat_edges.extend(E[E.nonzero()].tolist()[0])
			heat_edges.extend(E[E.nonzero()].tolist())
	
	# Sort the collected edge values in descending order
	heat_edges = sorted(heat_edges, reverse=True)

	# Calculate the maximum allowable difference for binary search termination
	diff = max(2 / len(heat_edges), 0.00001)

	# Check if the initial module sizes are within the size cutoff
	module_size = []
	for E_copy in copy.deepcopy(E_list):
		# G_tmp = nx.from_numpy_matrix(E_copy.transpose(), create_using=nx.DiGraph)
		G_tmp = nx.from_numpy_array(E_copy.transpose(), create_using=nx.DiGraph)
		for element in list(nx.strongly_connected_components(G_tmp)):
			module_size.append(len(element))

	# If initial module sizes meet the cutoff, set final_delta and final_midpoint accordingly
	if max(module_size) <= size_cutoff:
		final_delta = -1
		final_midpoint = 1
	else:
		# Use binary search to determine the optimal delta value
		while upper_bound - lower_bound > diff:
			midpoint = (upper_bound + lower_bound) / 2
			delta = heat_edges[int(midpoint * len(heat_edges))]
			
			# Evaluate the module sizes with the current delta
			module_size = []
			for E_copy in copy.deepcopy(E_list):
				E_copy[E_copy <= delta] = 0
				# G_tmp = nx.from_numpy_matrix(E_copy.transpose(), create_using=nx.DiGraph)
				G_tmp = nx.from_numpy_array(E_copy.transpose(), create_using=nx.DiGraph)
				for element in list(nx.strongly_connected_components(G_tmp)):
					module_size.append(len(element))
			
			# Adjust search bounds based on module size
			if max(module_size) <= size_cutoff:
				lower_bound = midpoint
			else:
				upper_bound = midpoint

		# Determine final delta value
		final_midpoint = lower_bound
		final_delta = heat_edges[int(final_midpoint * len(heat_edges)) + 1]

	# Write the results to the output file
	with open(output, "a+") as f:
		f.write("\t".join([str(final_delta), str(final_midpoint), ",".join([str(thing) for thing in module_size])]) + "\n")
	
	return


def get_one_delta(G_original, delta_output, restart_prob, max_subnetwork_size, seed):
	"""
	Determines an appropriate delta for a given initialized graph by shuffling edge weights and performing heat diffusion.

	Parameters
	----------
	G_original : networkx.Graph
		The original undirected graph with nodes and edges.
	delta_output : str
		Path to the output file where the determined delta value will be saved.
	restart_prob : float
		Restart probability for the heat diffusion process.
	max_subnetwork_size : int
		The maximum allowable size for an interconnected module in a heat diffusion network derived from a randomly shuffled graph.
	seed : int
		Seed for the random number generator to ensure reproducibility.
	"""
	
	# Create a deep copy of the original graph to preserve it
	G = copy.deepcopy(G_original)
	
	# Extract the list of edge weights from the graph
	weight_list = [info["weight"] for u, v, info in G.edges(data=True)]
	
	# Shuffle the list of edge weights with a given seed for reproducibility
	random.seed(seed)
	random.shuffle(weight_list)
	
	# Perform double edge swap to randomize the graph structure while preserving degree sequence
	nx.double_edge_swap(G, nswap=len(G.edges()), max_tries=500000, seed=seed)
	
	# Assign the shuffled weights back to the edges of the graph
	for i, (u, v) in enumerate(G.edges()):
		G[u][v]["weight"] = weight_list[i]
	
	# Determine the suitable delta value for the randomized graph
	get_suitable_delta(G, beta=restart_prob, size_cutoff=max_subnetwork_size, output=delta_output)
	
	return


def remove_whole_dir(dir_name):
	"""
	Removes an entire directory and all its contents.
	"""
	shutil.rmtree(dir_name) 
	return



