description='''
- Functions used by NetFlow3D.py
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
from urllib.request import urlopen
from scipy.stats import poisson

def extract_uniprot_ids(fasta_file):
	uniprot_ids = []
	with open(fasta_file, 'r') as file:
		for line in file:
			if line.startswith('>'):
				# Split the header line and extract the UniProt ID
				parts = line.split('|')
				uniprot_id = parts[1]
				uniprot_ids.append(uniprot_id)
	return uniprot_ids

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
	# pval = scipy.stats.norm.sf(abs(z)) * 2 # Two-sided p-value
	pval = scipy.stats.norm.sf(z)  # One-sided p-value for greater
	return [enr, SE_enr, pval]


def df_to_dict(df_file, header = None, from_col = None, to_col = None):
	def df_to_dict(row, from_col, to_col):
		df2dict[row[from_col]] = row[to_col]
		return True
	df = pd.read_csv(df_file, sep = "\t", header = header)
	df2dict = {}
	df["status"] = df.apply(lambda row: df_to_dict(row, from_col, to_col), axis = 1)
	return df2dict


def get_prolen(uniprot_list):
	uniprot2len = {"O60344":811,"P62158":149,"P01233":165,"Q13748":450,"Q9P0W5":563,"Q5VU13":414,"Q9Y2S0":122,"Q3BBV1":5207,"Q30KQ2":79,"P30042":268,"A2BFH1":164, "Q9H2Q1":129, "Q494X1":144, "Q9UGB4":50, "Q9UBA6":126, "Q5JXX5":376, "Q9BZ01":182, "P0C7V5":454, "Q8WXC7":389, "Q4KN68":188, "Q5T1B1": 145, "Q9Y4X1": 527, "Q8IVM7": 164, "Q15513": 63, "B9ZVM9":353, "P04745": 511}
	for element in uniprot_list:
		if element not in uniprot2len:
			try:
				html = urlopen("http://www.uniprot.org/uniprot/" + element + ".fasta")
				file = html.read().decode('utf-8')
				length = float(len(''.join(file.strip('\n').split('\n')[1:])))
				if length > 0:
					a = length
				else:
					logging.warning("Unable to find the sequence length of " + element + "!")
					a = np.nan
			except:
				logging.warning("Unable to find the sequence length of " + element + "!")
				a = np.nan
			uniprot2len[element] = a
	return(uniprot2len)


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


def get_cluster_mutations(res, df_mis):
	uniprot2res = defaultdict(set)
	for x in res.split(","):
		uniprot, pos = x.split("_")
		uniprot2res[uniprot] = uniprot2res[uniprot].union({int(pos)})
	df_mis = df_mis[df_mis["UniProt"].isin(uniprot2res) & df_mis.apply(lambda x: set(split_consecutive_pos(x["Protein_position"])).intersection(uniprot2res[x["UniProt"]]) != set())]
	return ",".join(sorted(set(df_mis["#Uploaded_variation"])))


def get_32_trinuc(trinuc, complement = {"A":"T", "T":"A", "C":"G", "G":"C"}):
	if trinuc[1] in {"A", "G"}:
		output = "".join([complement[element] for element in trinuc][::-1])
	else:
		output = trinuc
	return output


def get_res_annotation(df, output_path):
	"""
	Generate ShortVersion_mutation_data.txt, MUg.txt, and Patient_list.txt
	"""
	def gather_mutation_info(row):
		for element in split_consecutive_pos(row["Protein_position"]):
			res2mutcount[(row["UniProt"], element)] = res2mutcount[(row["UniProt"], element)] + 1
			res2mutinfo[(row["UniProt"], element)].append("_".join([row["Chromosome"], row["Start_Position"], row["Tumor_Sample_Barcode"]]))
		return True

	# gather mutation information
	res2mutcount = defaultdict(int)
	res2mutinfo = defaultdict(list)
	df["status"] = df.apply(gather_mutation_info, axis = 1)

	# ShortVersion_mutation_data.txt
	f = open(output_path + "ShortVersion_mutation_data.txt", "w")
	f.write("\t".join(["Uniprot", "Uniprot_Position", "Mutations", "Mutation_info"]) + "\n")
	for key in res2mutcount:
		f.write("\t".join([key[0], str(key[1]), str(res2mutcount[key]), ",".join(sorted(set(res2mutinfo[key])))]) + "\n")
	f.close()

	# # MUg.txt
	# mutrate = pd.read_csv(bg_file, sep = "\t").dropna()
	# average_rate = mutrate["MUg"].mean()
	# mutrate = pd.merge(df[["Hugo_Symbol", "UniProt"]].dropna().drop_duplicates(), mutrate, on = "Hugo_Symbol")
	# df_output = {"UniProt": [], "MUg": []}
	# for x in set(df["UniProt"]):
	# 	df_output["UniProt"].append(x)
	# 	df_one = mutrate[mutrate["UniProt"] == x]
	# 	if df_one.shape[0] == 0:
	# 		df_output["MUg"].append(average_rate)
	# 	elif df_one.shape[0] > 1:
	# 		df_output["MUg"].append(df_one["MUg"].mean())
	# 	else:
	# 		df_output["MUg"].append(df_one["MUg"].tolist()[0])
	# df_output = pd.DataFrame(df_output)
	# df_output.to_csv(output_path + "MUg.txt", sep = "\t", header = True, index = None)
	
	# f = open(output_path + "MUg.txt", "w")
	# f.write("\t".join(["UniProt", "MUg"]) + "\n")
	# for key in df["UniProt"].unique():
	# 	if key in uniprot2mutcount:
	# 		quantile = float(pd.DataFrame(list(uniprot2mutcount[key].values())).quantile(q = 0.99, axis = 0))
	# 		cutoff = max(quantile, 20 * uniprot2gene.loc[key]["Gene"])
	# 		i = 0
	# 		j = 0
	# 		for res in uniprot2mutcount[key]:
	# 			if uniprot2mutcount[key][res] > cutoff:
	# 				i = i + 1
	# 				j = j + uniprot2mutcount[key][res]
	# 		mut_freq = (sum(uniprot2mutcount[key].values()) - j) / (len(sample_list) * uniprot2gene.loc[key]["Gene"] * (prolen_dict[key] - i))
	# 		f.write("\t".join([key, str(mut_freq)]) + "\n")
	# 	else:
	# 		f.write("\t".join([key, str(0.1)]) + "\n")
	# f.close()

	return output_path + "ShortVersion_mutation_data.txt"


# def get_res_mutability(shortversion, uniprot2gene, gene2mutability, uniprot2mutability, trinuc2mutability, output):
# 	df = pd.read_csv(shortversion, sep = "\t")
# 	mg = {}
# 	for element in df["Uniprot"].unique():
# 		gene_mut = []
# 		for item in uniprot2gene[element]:
# 			if item in gene2mutability:
# 				gene_mut.append(gene2mutability[item])
# 		if len(gene_mut) > 0:
# 			mg[element] = np.mean(gene_mut)
# 		else:
# 			mg[element] = np.mean(list(gene2mutability.values()))
# 	df["res_mutability"] = df.apply(lambda row: trinuc2mutability[row["Ref_Tri"]] / mg[row["Uniprot"]] * uniprot2mutability[row["Uniprot"]], axis = 1)
# 	df.to_csv(output, sep = "\t", index = None, header = True)
# 	return output


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


def get_cluster_from_graph(graph_file, output):
	G = nx.read_graphml(graph_file)
	f = open(output, "w")
	for element in nx.connected_components(G):
		G_sub = G.subgraph(element)
		flag = "intra"
		for u,v,info in G_sub.edges(data = True):
			if info["type"] == "inter":
				flag = "inter"
				break
		f.write("\t".join([",".join(sorted(element)), flag]) + "\n")
	f.close()
	return output


def get_cluster_from_interface(uniprot2res, ires_file, binary_interactome, output):
	def get_mutated_cluster(row, uniprot2res):
		output1 = []
		if type(row["P1_IRES"]) == str:
			for element in row["P1_IRES"].split(","):
				uniprot, res = element.split("_")
				if int(res) in uniprot2res[uniprot]:
					output1.append(element)
		output2 = []
		if type(row["P2_IRES"]) == str:
			for element in row["P2_IRES"].split(","):
				uniprot, res = element.split("_")
				if int(res) in uniprot2res[uniprot]:
					output2.append(element)
		return ",".join(sorted(set(output1 + output2)))

	# get interface residues
	df_pioneer = pd.read_csv(ires_file, sep = "\t", dtype = str).dropna(subset = ["P1_IRES", "P2_IRES"], how = "all")
	df_pioneer = df_pioneer[df_pioneer["Interaction"].isin(binary_interactome)]

	# get mutated clusters
	if df_pioneer.shape[0] == 0:
		return None
	else:
		df_pioneer["Mutated_residues"] = df_pioneer.apply(lambda x: get_mutated_cluster(x, uniprot2res), axis = 1)
		df_pioneer = df_pioneer[df_pioneer["Mutated_residues"] != ""]
		df_pioneer["Structure_source"] = "PIONEER"
		df_pioneer[["Mutated_residues", "Structure_source", "Interaction"]].to_csv(output, sep = "\t", header = None, index = None)
		return output


# def get_cluster_from_interface(uniprot2res, ires_file, binary_interactome, output):
# 	def get_mutated_cluster(row, uniprot2res):
# 		output1 = []
# 		for element in row["P1_IRES"].split(","):
# 			uniprot, res = element.split("_")
# 			if int(res) in uniprot2res[uniprot]:
# 				output1.append(element)
# 		output2 = []
# 		for element in row["P2_IRES"].split(","):
# 			uniprot, res = element.split("_")
# 			if int(res) in uniprot2res[uniprot]:
# 				output2.append(element)
# 		if len(output1) > 0 and len(output2) > 0:
# 			return ",".join(sorted(set(output1 + output2)))
# 		else:
# 			return np.nan

# 	# get interface residues
# 	df_pioneer = pd.read_csv(ires_file, sep = "\t", dtype = str).dropna(subset = ["P1_IRES", "P2_IRES"])
# 	df_pioneer = df_pioneer[df_pioneer["Interaction"].isin(binary_interactome)]
# 	test_num = df_pioneer.shape[0]

# 	# get mutated clusters
# 	if df_pioneer.shape[0] == 0:
# 		return None
# 	else:
# 		df_pioneer["Mutated_residues"] = df_pioneer.apply(lambda x: get_mutated_cluster(x, uniprot2res), axis = 1)
# 		df_pioneer = df_pioneer.dropna(subset = ["Mutated_residues"])
# 		if df_pioneer.shape[0] > 0:
# 			df_pioneer["Structure_source"] = "PIONEER"
# 			df_pioneer[["Mutated_residues", "Structure_source", "Interaction"]].to_csv(output, sep = "\t", header = None, index = None)
# 			return [output, test_num]
# 		else:
# 			return None


def get_cluster_pval(mutres_file, cluster_file, sample_size, output, uniprot2mutcount, prolen_dict, factor_for_unknown_genes):
	def get_uniprot(res):
		uniprot = sorted(set([item.split("_")[0] for item in res.split(",")]))
		return ",".join(uniprot)
	def get_cluster_value(res, dictionary):
		if type(res) == str:
			values = []
			for element in res.split(","):
				values.append(str(dictionary[element]))
			return ",".join(values)
		else:
			return np.nan
	def get_total_mut(res, dictionary):
		if type(res) == str:
			output = []
			for element in res.split(","):
				output.extend(dictionary[element])
			return len(set(output))
		else:
			return 0
	def get_expected_mutcount(res):
		uniprot2res = defaultdict(list)
		if type(res) == str:
			for x in res.split(","):
				uniprot, pos = x.split("_")
				uniprot2res[uniprot].append(pos)
		expected = 0
		for x in uniprot2res:
			if x in uniprot2mutcount:
				expected = expected + uniprot2mutcount[x] * sample_size * (len(uniprot2res[x]) / prolen_dict[x])
			else:
				expected = expected + factor_for_unknown_genes * len(uniprot2res[x]) * 3 * sample_size
		return expected

	# load mutated residues
	df_mut = pd.read_csv(mutres_file, sep = "\t", dtype = {"Uniprot_Position": int, "Mutations": str})

	# annotate background mutability and observed mutation count for the clusters
	df_mut["res"] = df_mut.apply(lambda row: row["Uniprot"] + "_" + str(row["Uniprot_Position"]), axis = 1)
	res2mutcount = {}
	res2mutinfo = {}
	res_mutcount = df_mut["Mutations"].tolist()
	res_mutinfo = df_mut["Mutation_info"].tolist()
	for i, element in enumerate(df_mut["res"].tolist()):
		res2mutcount[element] = res_mutcount[i]
		res2mutinfo[element] = res_mutinfo[i].split(",")
	df_cluster = pd.read_csv(cluster_file, sep = "\t", header = None).rename(columns = {0: "Residues", 1: "Structure_source", 2: "Uniprots"})
	if "Uniprots" not in df_cluster.columns:
		df_cluster["Uniprots"] = df_cluster["Residues"].apply(get_uniprot)
	df_cluster["Mutation_count"] = df_cluster["Residues"].apply(lambda x: get_cluster_value(x, res2mutcount))
	df_cluster["Observed_mut"] = df_cluster["Residues"].apply(lambda x: get_total_mut(x, res2mutinfo))
	df_cluster["Expected_mut"] = df_cluster["Residues"].apply(get_expected_mutcount)
	df_cluster["Raw_pvalue"] = df_cluster.apply(lambda row: poisson.sf(row["Observed_mut"] - 1, row["Expected_mut"]), axis = 1)
	# df_cluster[["InFrame_enrichment", "se", "Raw_pvalue"]] = pd.DataFrame(df_cluster.apply(lambda row: calc_enr(row["Observed_mut"], total_mut_observed, row["Expected_mut"], total_mut_expected), axis = 1).tolist(), index = df_cluster.index)
	df_cluster.to_csv(output, sep = "\t", index = None, header = True)
	return output


def get_sig_cluster_structure(pval_file_PDB_intra, pval_file_PDB_inter, pval_file_AlphaFold2, all_uniprots, prolen_dict):
	def whether_retain(res, uniprot2res):
		for x in res.split(","):
			uniprot, pos = x.split("_")
			if pos in uniprot2res[uniprot]:
				return False
		return True

	uniprot2res = defaultdict(list)
	if pval_file_PDB_intra != None:
		df_pdb_intra = pd.read_csv(pval_file_PDB_intra, sep = "\t", dtype = {"Mutation_count": str})
		df_pdb_intra["Type"] = "InFrame_IntraProtein_Cluster"
		for x in df_pdb_intra["Residues"]:
			for y in x.split(","):
				uniprot, pos = y.split("_")
				uniprot2res[uniprot].append(pos)
	else:
		df_pdb_intra = pd.DataFrame()
	if pval_file_AlphaFold2 != None:
		df_af2 = pd.read_csv(pval_file_AlphaFold2, sep = "\t", dtype = {"Mutation_count": str})
		df1 = df_af2[~df_af2["Uniprots"].isin(uniprot2res)]
		df2 = df_af2[df_af2["Uniprots"].isin(uniprot2res)]
		df2 = df2[df2["Residues"].apply(lambda x: whether_retain(x, uniprot2res))]
		df_af2 = pd.concat([df1, df2])
		df_af2["Type"] = "InFrame_IntraProtein_Cluster"
		for x in df_af2["Residues"]:
			for y in x.split(","):
				uniprot, pos = y.split("_")
				uniprot2res[uniprot].append(pos)
	else:
		df_af2 = pd.DataFrame()
	if pval_file_PDB_inter != None:
		df_pdb_inter = pd.read_csv(pval_file_PDB_inter, sep = "\t", dtype = {"Mutation_count": str})
		df_pdb_inter["Type"] = "InFrame_InterProtein_Cluster"
	else:
		df_pdb_inter = pd.DataFrame()
	df_final = pd.concat([df_pdb_intra, df_af2, df_pdb_inter])
	test_num = df_final.shape[0]
	# # -------------------------------
	# # consider whether add this
	# for x in all_uniprots:
	# 	test_num = test_num + (prolen_dict[x] - len(uniprot2res[x]))
	# # -------------------------------
	if df_final.shape[0] > 0:
		df_final["Adjusted_pvalue"] = df_final["Raw_pvalue"].apply(lambda x: min(x * test_num, 1))
		df_final = df_final[["Type", "Structure_source", "Uniprots", "Residues", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
	else:
		df_final = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutation_count": [], \
			"Raw_pvalue": [], "Adjusted_pvalue": []})
	return df_final


def get_sig_cluster_interface(pval_file_PIONEER):
	df_pio = pd.read_csv(pval_file_PIONEER, sep = "\t", dtype = {"Mutation_count": str})
	df_pio["Adjusted_pvalue"] = df_pio["Raw_pvalue"].apply(lambda x: min(x * df_pio.shape[0], 1))
	df_pio["Type"] = "InFrame_InterProtein_Cluster"
	return df_pio[["Type", "Structure_source", "Uniprots", "Residues", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]


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


def format_graph(G):
	"""
	Sort the nodes and edges in a graph.
	"""
	H = nx.Graph()
	H.add_nodes_from(sorted(G.nodes(data = True)))
	H.add_edges_from(sorted([(sorted([u,v])[0], sorted([u,v])[1], info) for u,v,info in G.edges(data = True)]))
	return H


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

def get_initial_distribution_alternate(G_original, final_output_intra_pdb, final_output_intra_af2, final_output_inter_pdb, final_output_inter_pioneer, final_output_intra_lof, output, intercept, weighted_edge = True):
	def get_uniprot2pval(row):
		if type(row["Residues"]) == str:
			for uniprot in {x.split("_")[0] for x in row["Residues"].split(",")}:
				uniprot2pval[uniprot].append(row["Raw_pvalue"])
		return True
	def get_interaction2pval(row):
		interaction2pval[tuple(sorted(row["Uniprots"].split(",")))].append(row["Raw_pvalue"])
		return True

	# initialization
	G = copy.deepcopy(G_original)
	for u,v in G.edges():
		G[u][v]["weight"] = intercept
	for u in G.nodes():
		G.nodes[u]["heat_score"] = 0.0

	# set edge weight
	if weighted_edge == True:
		interaction2pval = defaultdict(list)
		for file in [final_output_inter_pdb, final_output_inter_pioneer]:
			df = pd.read_csv(file, sep = "\t")
			df["status"] = df.apply(get_interaction2pval, axis = 1)
		for u,v in interaction2pval:
			if G.has_edge(u,v):
				G[u][v]["weight"] = G[u][v]["weight"] + min(-np.log10(min(interaction2pval[(u,v)])), 300)

	# set node weight
	uniprot2pval = defaultdict(list)
	for file in [final_output_intra_pdb, final_output_intra_af2, final_output_inter_pdb, final_output_inter_pioneer]:
		df = pd.read_csv(file, sep = "\t")
		df["status"] = df.apply(get_uniprot2pval, axis = 1)
	for uniprot in set(uniprot2pval).intersection(set(G.nodes())):
		G.nodes[uniprot]["heat_score"] = min(-np.log10(min(uniprot2pval[uniprot])), 300)
	uniprot2pval = df_to_dict(final_output_intra_lof, header = 0, from_col = "Uniprots", to_col = "Raw_pvalue")
	for uniprot in set(uniprot2pval).intersection(set(G.nodes())):
		G.nodes[uniprot]["heat_score"] = G.nodes[uniprot]["heat_score"] + min(-np.log10(uniprot2pval[uniprot]), 300)

	# write output
	G = format_graph(G)
	nx.write_graphml(G, output)
	return G


def get_initial_distribution_alternate_sig(G_original, df_intra, df_inter, df_lof, intercept, output):
	def get_uniprot2pval(row):
		for uniprot in row["Uniprots"].split(","):
			uniprot2pval[uniprot].append(row["Raw_pvalue"])
		return True
	def get_interaction2pval(row):
		interaction2pval[tuple(sorted(row["Uniprots"].split(",")))].append(row["Raw_pvalue"])
		return True

	# initialization
	G = copy.deepcopy(G_original)
	for u,v in G.edges():
		G[u][v]["weight"] = intercept
	for u in G.nodes():
		G.nodes[u]["heat_score"] = 0.0

	# set edge weight
	interaction2pval = defaultdict(list)
	df_inter["status"] = df_inter.apply(get_interaction2pval, axis = 1)
	for u,v in interaction2pval:
		if G.has_edge(u,v):
			G[u][v]["weight"] = G[u][v]["weight"] + min(-np.log10(min(interaction2pval[(u,v)])), 300)

	# set node weight
	df_inframe = pd.concat([df_intra, df_inter])
	uniprot2pval = defaultdict(list)
	df_inframe["status"] = df_inframe.apply(get_uniprot2pval, axis = 1)
	for uniprot in uniprot2pval:
		if uniprot in G:
			G.nodes[uniprot]["heat_score"] = G.nodes[uniprot]["heat_score"] + min(-np.log10(min(uniprot2pval[uniprot])), 300)
	uniprot2pval = defaultdict(list)
	df_lof["status"] = df_lof.apply(get_uniprot2pval, axis = 1)
	for uniprot in uniprot2pval:
		if uniprot in G:
			G.nodes[uniprot]["heat_score"] = G.nodes[uniprot]["heat_score"] + min(-np.log10(min(uniprot2pval[uniprot])), 300)

	# write output
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
	weight_list = [info["weight"] for u,v,info in G.edges(data = True)]
	random.seed(seed)
	random.shuffle(weight_list)
	nx.double_edge_swap(G, nswap = len(G.edges()), max_tries = 500000, seed = seed)
	for i,(u,v) in enumerate(G.edges()):
		G[u][v]["weight"] = weight_list[i]
	get_suitable_delta(G, beta = restart_prob, size_cutoff = max_subnetwork_size, output = delta_output)
	return


def remove_whole_dir(dir_name):
	shutil.rmtree(dir_name)
	return


