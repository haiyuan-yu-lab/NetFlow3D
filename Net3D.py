description = '''
- Net3D: identification of functional modules through structurome-guided network diffusion (currently only the human structurome is available)
- Dependencies:
	funcs.py
- 
'''

import os
import argparse
from multiprocessing import Pool
import networkx as nx
import pandas as pd
from collections import defaultdict
import json
import logging
import numpy as np
import funcs
from os import path
import copy
import statsmodels.stats.multitest
import sys

def get_clusterID_for_a_mis(mut_uniprot, mut_pos, df_cluster):
	def get_cluster_with_mut(res, mut_uniprot, mut_pos):
		for x in res.split(","):
			uniprot, pos = x.split("_")
			if uniprot == mut_uniprot and int(pos) in funcs.split_consecutive_pos(mut_pos):
				return True
		return False

	df = df_cluster[df_cluster["Uniprots"].apply(lambda x: mut_uniprot in x.split(","))]
	df = df[df["Residues"].apply(lambda x: get_cluster_with_mut(x, mut_uniprot, mut_pos))]
	return ",".join(df["Signature_ID"].tolist())

def get_clusterID_for_a_lof(uniprot, df_cluster):
	return df_cluster[df_cluster["Uniprots"] == uniprot]["Signature_ID"].tolist()[0]

def get_graph_by_uniprot_binary_multirun(args):
	return funcs.get_graph_by_uniprot_binary(*args)

def get_cluster_from_graph_multirun(args):
	return funcs.get_cluster_from_graph(*args)

def get_clusterID_for_a_mis_multirun(args):
	return get_clusterID_for_a_mis(*args)

def get_clusterID_for_a_lof_multirun(args):
	return get_clusterID_for_a_lof(*args)

def get_one_delta_multirun(args):
	return funcs.get_one_delta(*args)

if __name__ == "__main__":
	# user input
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument('-m','--input_maf', required = True, type = str, help = 'Mutation data in MAF format (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)')
	parser.add_argument('-R','--resolution', required = True, type = str, help = 'high or low. If set to high, smaller subnetworks are favored. If set to low, larger subnetworks are favored.')
	parser.add_argument('-I','--job_name', required = True, type = str, help = 'Please specify a name for your job. The output files will be stored in a folder with this name')
	parser.add_argument('-t','--threads', type = int, default = 5, help = '[OPTIONAL] The number of threads to use. By default 5')
	parser.add_argument('-X','--expressed_genes', type = str, default = None, help = '[OPTIONAL] A text file with the genes expressed in your context. Each row stands for a gene (HUGO/ENSP/UniProt IDs are accepted). By default all genes are considered as expressed')
	parser.add_argument('-n','--binary_interactome', type = str, default = "./metadata/HomoSapiens_interactome_HINThq.txt", help = '[OPTIONAL] A text file with protein-protein interactions. Each row stands for an interaction (IDs should be separated by tab. HUGO/ENSP/UniProt IDs are accepted). Please include all existing interactions in your context. By default we use the human binary interactome from HINT (http://hint.yulab.org/download/HomoSapiens/binary/hq/)')
	parser.add_argument('-o','--output_path', type = str, default = "./output/", help = '[OPTIONAL] Path to the job folder. If not specified, the job folder will be stored in the default output folder')
	parser.add_argument('-L','--logfile_path', type = str, default = "./log/", help = '[OPTIONAL] Path to the log file. If not specified, the log file will be stored in the default log folder')
	
	# native files
	parser.add_argument('-M','--gene_mutability_file', type = str, default = "./metadata/Gene_mutability.txt", help = '[OPTIONAL] Average mutability of trinucleotides in each gene (native file)')
	parser.add_argument('-T','--trinuc_mutability_file', type = str, default = "./metadata/Tri_Mutability.txt", help = '[OPTIONAL] Exome-wide mutability of 32 trinucleotides (native file)')
	parser.add_argument('-l','--prolen_file', type = str, default = "./metadata/uniprot2prolen.json", help = '[OPTIONAL] A sequence length map of human proteins (native file)')
	parser.add_argument('-g','--gene2uniprot_file', type = str, default = "./metadata/gene2uniprot.json", help = '[OPTIONAL] A UniProt ID map of human genes (native file)')
	parser.add_argument('-P','--ensp2uniprot_file', type = str, default = "./metadata/ensp2uniprot.json", help = '[OPTIONAL] ID conversion between ENSP and UniProt (native file)')
	parser.add_argument('-G','--ensg2uniprot_file', type = str, default = "./metadata/ensg2uniprot.json", help = '[OPTIONAL] ID conversion between ENSG and UniProt (native file)')
	parser.add_argument('-c','--enst2uniprot_file', type = str, default = "./metadata/enst2uniprot.json", help = '[OPTIONAL] ID conversion between ENST and UniProt (native file)')
	parser.add_argument('-i','--ires_file', type = str, default = "./metadata/HomoSapiens_interfaces_PIONEER.txt", help = '[OPTIONAL] Protein-protein interaction interfaces predicted by PIONEER (native file)')
	parser.add_argument('-a','--PDB_intra_resource', type = str, default = "./graph/PDB_intra/", help = '[OPTIONAL] Intra-chain residue-residue distances derived from PDB structures (native resource)')
	parser.add_argument('-e','--PDB_inter_resource', type = str, default = "./graph/PDB_inter/", help = '[OPTIONAL] Inter-chain residue-residue distances derived from PDB structures (native resource)')
	parser.add_argument('-d','--AF2_intra_resource', type = str, default = "./graph/AF2_pLDDT0/", help = '[OPTIONAL] Intra-chain residue-residue distances derived from AlphaFold structures (native resource)')
	
	# default settings
	parser.add_argument('-s','--significance_level', type = float, default = 0.05, help = '[OPTIONAL] Significance level. By default 0.05')
	parser.add_argument('-A','--intra_dist_upperlimit', type = float, default = 6.0, help = '[OPTIONAL] Intra-chain distance cutoff. By default 6 angstrom')
	parser.add_argument('-E','--inter_dist_upperlimit', type = float, default = 9.0, help = '[OPTIONAL] Inter-chain distance cutoff. By default 9 angstrom')
	parser.add_argument('-H','--designated_score', type = float, default = 1.0, help = '[OPTIONAL] The designated score per node per type of selection signature. By default 1')
	parser.add_argument('-w','--enhanced_conductivity', type = float, default = 5.0, help = '[OPTIONAL] The weight of edges with selection signatures. By default 5')
	parser.add_argument('-r','--restart_prob', type = float, default = 0.5, help = '[OPTIONAL] Restart probability in the diffusion model. By default 0.5')
	parser.add_argument('-x','--max_subnetwork_size', type = int, default = 5, help = '[OPTIONAL] Subnetwork size cutoff used for determining appropriate delta. By default 5')
	parser.add_argument('-D','--delta_trial', type = int, default = 20, help = '[OPTIONAL] Round of permutations when determining the value of delta. By default 20')
	args = parser.parse_args()

	def get_uniprot2res(row):
		uniprot2res[row["Uniprot"]] = uniprot2res[row["Uniprot"]].union({row["Uniprot_Position"]})
		return True
	def get_affected_genes(uniprots, uniprot2gene):
		output = []
		for x in uniprots.split(","):
			output.append("/".join(sorted(uniprot2gene[x])))
		return ",".join(sorted(output))
	def get_mis_mutfreq(row, uniprot2gene):
		res = row["Residues"].split(",")
		mutcount = row["Mutation_count"].split(",")
		output = []
		for i in range(len(res)):
			uniprot, pos = res[i].split("_")
			output.append("/".join(sorted(uniprot2gene[uniprot])) + "_" + pos + ":" + mutcount[i])
		return ",".join(output)
	def get_module_for_a_mut(uniprot, uniprot2module):
		if uniprot in uniprot2module:
			size = len(uniprot2module[uniprot].split("_"))
			return [uniprot2module[uniprot], size]
		else:
			return ["[NA]", "[NA]"]
	def get_subnetwork(affected_genes, subnetwork_list):
		for x in subnetwork_list:
			if set(affected_genes.split(",")) - set(x.split(",")) == set():
				return x
		return np.nan	

	############################################################################################################################################
	# Check input files
	# -----------------
	if os.path.isdir(args.output_path) == False:
		sys.exit("No such directory: " + args.output_path)
	if os.path.isdir(args.logfile_path) == False:
		sys.exit("No such directory: " + args.logfile_path)
	if path.exists(args.input_maf) == False:
		sys.exit("No such file or directory: " + args.input_maf)
	if path.exists(args.binary_interactome) == False:
		sys.exit("No such file or directory: " + args.binary_interactome)
	if args.expressed_genes != None and path.exists(args.expressed_genes) == False:
		sys.exit("No such file or directory: " + args.expressed_genes)
	if args.resolution not in {"high", "low"}:
		sys.exit("Please check --resolution! Can only be set to high or low.")
	############################################################################################################################################

	

	############################################################################################################################################
	# Create the job
	# --------------
	if args.output_path[-1] == "/":
		output_path = args.output_path + args.job_name + "/"
	else:
		output_path = args.output_path + "/" + args.job_name + "/"
	try:
		os.mkdir(output_path)
	except:
		sys.exit("The job name is already in use! Please specify another name.")
	if args.logfile_path[-1] == "/":
		logging.basicConfig(filename = args.logfile_path + args.job_name + '.log', filemode = 'w', level = logging.DEBUG, format = '%(levelname)s - %(message)s')
	else:
		logging.basicConfig(filename = args.logfile_path + "/" + args.job_name + '.log', filemode = 'w', level = logging.DEBUG, format = '%(levelname)s - %(message)s')
	############################################################################################################################################



	############################################################################################################################################
	# Load native data
	# ----------------
	with open(args.prolen_file) as f:
		prolen_dict = json.load(f)
	with open(args.ensp2uniprot_file) as f:
		ensp2uniprot = json.load(f)
	with open(args.gene2uniprot_file) as f:
		gene2uniprot = json.load(f)
	with open(args.ensg2uniprot_file) as f:
		ensg2uniprot = json.load(f)
	with open(args.enst2uniprot_file) as f:
		enst2uniprot = json.load(f)
	############################################################################################################################################



	############################################################################################################################################
	# Prepare mutation data
	# ---------------------
	try:
		df = pd.read_csv(args.input_maf, sep = "\t", dtype = str).dropna(subset = ["Hugo_Symbol"])
	except:
		funcs.remove_whole_dir(output_path)
		sys.exit("Please check the input MAF file! (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
	if {"Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "Protein_position", "Codons", "ENSP"} - set(df.columns.tolist()) != set():
		logging.error("Please check the input MAF file! The following columns must be included: Hugo_Symbol, Tumor_Sample_Barcode, Codons, ENSP, Protein_position, Variant_Classification (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
		funcs.remove_whole_dir(output_path)
		sys.exit()
	if set(df["Variant_Classification"].dropna().tolist()) - {"Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "RNA", "Targeted_Region", "Intron", "IGR", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame"} != set():
		logging.error("Please check the content in the column named Variant_Classification in the input MAF file! (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
		funcs.remove_whole_dir(output_path)
		sys.exit()
	df = df.dropna(subset = ["Tumor_Sample_Barcode"])
	if df.shape[0] == 0:
		logging.error("Please specify the sample ID for each mutation!")
		funcs.remove_whole_dir(output_path)
		sys.exit()
	df["UniProt"] = df.apply(lambda x: funcs.get_uniprot(x, gene2uniprot, ensp2uniprot, enst2uniprot, ensg2uniprot), axis = 1)
	df = df.dropna(subset = ["UniProt"])
	if df.shape[0] == 0:
		logging.warning("Could not map the mutations to any protein!")
		funcs.remove_whole_dir(output_path)
		sys.exit()
	############################################################################################################################################
	
	

	############################################################################################################################################
	# Build protein-protein interaction network
	# -----------------------------------------
	try:
		df_interactome = pd.read_csv(args.binary_interactome, sep = "\t", header = None)
		pro1_list = df_interactome[0].tolist()
		pro2_list = df_interactome[1].tolist()
	except:
		funcs.remove_whole_dir(output_path)
		sys.exit("Please check the file with binary interactome! Each row must contain two HUGO/ENSP/UniProt IDs separated by tab. See example for details.")
	edges = set()
	for i in range(df_interactome.shape[0]):
		pro1 = pro1_list[i]
		pro2 = pro2_list[i]
		if pro1.split(".")[0] in ensp2uniprot:
			pro1 = ensp2uniprot[pro1.split(".")[0]]
		elif pro1 in gene2uniprot:
			pro1 = gene2uniprot[pro1]
		else:
			pro1 = pro1.split(".")[0].split("-")[0]
		if pro2.split(".")[0] in ensp2uniprot:
			pro2 = ensp2uniprot[pro2.split(".")[0]]
		elif pro2 in gene2uniprot:
			pro2 = gene2uniprot[pro2]
		else:
			pro2 = pro2.split(".")[0].split("-")[0]
		edges = edges.union({tuple(sorted([pro1, pro2]))})
	G = nx.Graph(edges)
	############################################################################################################################################

	
	
	############################################################################################################################################
	# Filter out mutations in unexpressed genes and remove unexpressed proteins from the network
	# ------------------------------------------------------------------------------------------
	if args.expressed_genes != None:
		try:
			expr_genes = set(pd.read_csv(args.expressed_genes, sep = "\t", header = None)[0].tolist())
		except:
			funcs.remove_whole_dir(output_path)
			sys.exit("Please check the format of the file with expressed genes! Each row should only contain one HUGO/ENSP/UniProt ID. See example for details.")
		expr_uniprots = set()
		for x in expr_genes:
			if x.split(".")[0] in ensp2uniprot:
				expr_uniprots = expr_uniprots.union({ensp2uniprot[x.split(".")[0]]})
			elif x.split(".")[0] in enst2uniprot:
				expr_uniprots = expr_uniprots.union(set(enst2uniprot[x.split(".")[0]]))
			elif x.split(".")[0] in ensg2uniprot:
				expr_uniprots = expr_uniprots.union(set(ensg2uniprot[x.split(".")[0]]))
			elif x in gene2uniprot:
				expr_uniprots = expr_uniprots.union({gene2uniprot[x]})
			else:
				expr_uniprots = expr_uniprots.union({x.split(".")[0].split("-")[0]})
		G = nx.Graph(G.subgraph(expr_uniprots))
		df = df[df["UniProt"].isin(expr_uniprots)]
		if df.shape[0] == 0:
			logging.warning("No mutation is in expressed genes!")
			funcs.remove_whole_dir(output_path)
			sys.exit()
	############################################################################################################################################

	

	############################################################################################################################################
	# Map genes to UniProts
	# ---------------------
	uniprot2gene = defaultdict(set)
	gene_list = df["Hugo_Symbol"].tolist()
	uniprot_list = df["UniProt"].tolist()
	for i in range(df.shape[0]):
		uniprot2gene[uniprot_list[i]] = uniprot2gene[uniprot_list[i]].union({gene_list[i]})
	uniprot2gene = dict(uniprot2gene)

	# define loss-of-function(LOF) and non-truncating mutations
	lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
	mis = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}
	############################################################################################################################################


	
	############################################################################################################################################
	# Identify selection signatures formed by LOF mutations
	# -----------------------------------------------------
	df_lof = df[df["Variant_Classification"].isin(lof)]
	df_lof = df_lof[df_lof["UniProt"].isin(prolen_dict)]
	df_lof_mut = copy.deepcopy(df_lof)
	if df_lof.shape[0] > 0:
		logging.info(str(df_lof.shape[0]) + " available loss-of-function mutation(s).")
		# protein-specific enrichment analysis
		total_len = sum([prolen_dict[x] for x in df_lof["UniProt"].unique()])
		f = open(output_path + "All_intra_LoF_pvalue.txt", "w")
		f.write("\t".join(["Uniprots", "Mutated_sample", "Mutation_count", "LoF_enrichment", "Raw_pvalue"]) + "\n")
		for key, df_one in df_lof.groupby("UniProt"):
			enr, _, pval = funcs.calc_enr(df_one.shape[0], df_lof.shape[0], prolen_dict[key], total_len)  
			sample_list = []
			mutcount_list = []
			for sample, df_onesample in df_one.groupby("Tumor_Sample_Barcode"):
				sample_list.append(sample)
				mutcount_list.append(str(df_onesample.shape[0]))
			f.write("\t".join([key, ",".join(sample_list), ",".join(mutcount_list), str(enr), str(pval)]) + "\n")
		f.close()

		# multiple test correction
		df_lof = pd.read_csv(output_path + "All_intra_LoF_pvalue.txt", sep = "\t", dtype = {"Mutation_count": str})
		_, df_lof["Adjusted_pvalue"], _, _ = statsmodels.stats.multitest.multipletests(df_lof["Raw_pvalue"].tolist(), alpha = args.significance_level, \
			method = 'fdr_bh')

		# identify the proteins significantly enriched of LOF mutations
		df_lof = df_lof[(df_lof["Adjusted_pvalue"] < args.significance_level) & (df_lof["LoF_enrichment"] > 1)]
		df_lof["Type"] = "LoF_IntraProtein_Enriched"
		df_lof["LoF_enrichment"] = df_lof["LoF_enrichment"].apply(lambda x: str(round(x, 2)))
		df_lof = df_lof[["Type", "Uniprots", "Mutated_sample", "Mutation_count", "LoF_enrichment", "Raw_pvalue", "Adjusted_pvalue"]]
		lof_uniprots = df_lof["Uniprots"].tolist()
		logging.info(str(df_lof.shape[0]) + " proteins are significantly enriched of loss-of-function mutations.")
	else:
		logging.warning("No available loss-of-function mutation!")
		lof_uniprots = []
		df_lof = pd.DataFrame({"Type": [], "Uniprots": [], "Mutated_sample": [], "Mutation_count": [], "LoF_enrichment": [], "Raw_pvalue": [], "Adjusted_pvalue": []})
	############################################################################################################################################
	


	############################################################################################################################################
	# Identify selection signatures formed by non-truncating mutations
	# ----------------------------------------------------------------
	df_mis = df[df["Variant_Classification"].isin(mis)].dropna(subset = ["Protein_position", "Codons"])
	df_mis = df_mis[df_mis["Codons"].apply(lambda x: set(x.split("/")[0].upper()) - {"A", "T", "C", "G"} == set())]
	df_mis["Protein_position"] = df_mis["Protein_position"].apply(lambda x: x.split("/")[0])
	df_mis_mut = copy.deepcopy(df_mis)
	if df_mis.shape[0] == 0:
		logging.warning("No available non-truncating mutation!")
		skip_flag = 1
	else:
		shortversion, patient_list, mug = funcs.get_res_annotation(df_mis, prolen_dict, output_path)
		df = pd.read_csv(shortversion, sep = "\t", usecols = ["Uniprot", "Uniprot_Position", "Mutation_info"], dtype = {"Uniprot_Position": int})
		if df.shape[0] == 0:
			logging.warning("No available non-truncating mutation!")
			skip_flag = 1
		else:
			logging.info(str(df.shape[0]) + " available non-truncating mutation(s)!")
			skip_flag = 0

			# estimate the mutability of each mutated residue
			gene2mutability = funcs.df_to_dict(args.gene_mutability_file, header = None, from_col = 0, to_col = 1)
			trinuc2mutability = funcs.df_to_dict(args.trinuc_mutability_file, header = 0, from_col = "tri", to_col = "mu")
			uniprot2mutability = funcs.df_to_dict(mug, header = 0, from_col = "#Uniprot", to_col = "MUg")
			output = funcs.get_res_mutability(shortversion, uniprot2gene, gene2mutability, uniprot2mutability, trinuc2mutability, shortversion)

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
			# identify spatial clusters and evaluate their significance
			
			uniprot2res = defaultdict(set)
			df["status"] = df.apply(get_uniprot2res, axis = 1)
			n = pd.read_csv(patient_list, sep = "\t").shape[0]
			df = pd.read_csv(output, sep = "\t", usecols = ['Uniprot', 'Uniprot_Position', 'Mutations', 'res_mutability'], dtype = {"Uniprot_Position": str})
			
			# PDB
			os.mkdir(output_path + "PDB_graph")
			function_input = []
			for key in uniprot2res:
				function_input.append([key, uniprot2res, output_path + "PDB_graph/", args.PDB_intra_resource, args.PDB_inter_resource, args.intra_dist_upperlimit, args.inter_dist_upperlimit])
			p = Pool(args.threads)
			output = p.map(get_graph_by_uniprot_binary_multirun, function_input)
			p.close()
			function_input = []
			for element in output:
				function_input.extend(element)
			if len(function_input) > 0:
				os.mkdir(output_path + "PDB_cluster")
				function_input = [[x, output_path + "PDB_cluster/" + os.path.split(x)[1].split(".")[0] + ".txt"] for x in set(function_input)]
				p = Pool(args.threads)
				output = p.map(get_cluster_from_graph_multirun, function_input)
				p.close()
				df_intra = []
				df_inter = []
				for x in output:
					df = pd.read_csv(x, sep = "\t", header = None)
					uniprot_id = os.path.split(x)[1].split(".")[0].split("_")
					if len(uniprot_id) == 1:
						df_intra.append(df)
					else:
						df = df[df[1] == "inter"]
						if df.shape[0] > 0:
							df_inter.append(df)
				if len(df_intra) > 0:
					df_intra = pd.concat(df_intra)
					df_intra[1] = "PDB"
					df_intra.to_csv(output_path + "Intra_protein_clusters_PDB.txt", sep = "\t", header = None, index = None)
					final_output_intra_pdb = funcs.get_cluster_pval(shortversion, output_path + "Intra_protein_clusters_PDB.txt", n, output_path + "Atom_PDB_intra_cnm_pvalue.txt")
				else:
					final_output_intra_pdb = None
				if len(df_inter) > 0:
					df_inter = pd.concat(df_inter)
					df_inter[1] = "PDB"
					df_inter.to_csv(output_path + "Inter_protein_clusters_PDB.txt", sep = "\t", header = None, index = None)
					final_output_inter_pdb = funcs.get_cluster_pval(shortversion, output_path + "Inter_protein_clusters_PDB.txt", n, output_path + "Atom_PDB_inter_cnm_pvalue.txt")
					df = pd.read_csv(final_output_inter_pdb, sep = "\t")
					df["Uniprots"] = df["Uniprots"].apply(funcs.binary_interaction)
					df = df[df["Uniprots"].apply(lambda x: G.has_edge(x.split(",")[0], x.split(",")[1]))]
					if df.shape[0] > 0:
						df.to_csv(final_output_inter_pdb, sep = "\t", header = True, index = None)
					else:
						final_output_inter_pdb = None
				else:
					final_output_inter_pdb = None
			else:
				final_output_intra_pdb = None
				final_output_inter_pdb = None

			# AlphaFold2
			function_input = []
			os.mkdir(output_path + "AlphaFold2_graph_pLDDT0")
			for key in uniprot2res:
				function_input.append([key, uniprot2res, output_path + "AlphaFold2_graph_pLDDT0/", args.AF2_intra_resource, None, args.intra_dist_upperlimit, args.inter_dist_upperlimit])
			p = Pool(args.threads)
			output = p.map(get_graph_by_uniprot_binary_multirun, function_input)
			p.close()

			# cluster generation (AlphaFold2)
			function_input = []
			for element in output:
				function_input.extend(element)
			if len(function_input) > 0:
				os.mkdir(output_path + "AlphaFold2_cluster")
				function_input = [[x, output_path + "AlphaFold2_cluster/" + os.path.split(x)[1].split(".")[0] + ".txt"] for x in set(function_input)]
				p = Pool(args.threads)
				output = p.map(get_cluster_from_graph_multirun, function_input)
				p.close()
				df_intra = []
				for element in output:
					df_intra.append(pd.read_csv(element, sep = "\t", header = None))
				df_intra = pd.concat(df_intra)
				df_intra[1] = "AlphaFold2"
				df_intra.to_csv(output_path + "Intra_protein_clusters_AlphaFold2.txt", sep = "\t", header = None, index = None)
				final_output_intra_af2 = funcs.get_cluster_pval(shortversion, output_path + "Intra_protein_clusters_AlphaFold2.txt", n, output_path + "Atom_AlphaFold2_intra_cnm_pvalue_pLDDT0.txt")
			else:
				final_output_intra_af2 = None

			# PIONEER
			INSIDER_clusters = funcs.get_cluster_from_interface(uniprot2res, args.ires_file, [",".join(sorted([u,v])) for u,v in G.edges()], output_path + "PIONEER_clusters.txt")
			if INSIDER_clusters != None:
				final_output_inter_pioneer = funcs.get_cluster_pval(shortversion, INSIDER_clusters, n, output_path + "PIONEER_inter_pvalue.txt")
			else:
				final_output_inter_pioneer = None

			df_intra = funcs.get_sig_cluster_intra(final_output_intra_pdb, final_output_intra_af2, args.PDB_intra_resource, sig_cutoff = args.significance_level)
			df_inter = funcs.get_sig_cluster_inter(final_output_inter_pdb, final_output_inter_pioneer, sig_cutoff = args.significance_level)
	if skip_flag == 1:
		df_intra = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutated_sample": [], "Mutation_count": [], "Raw_pvalue": [], "Adjusted_pvalue": []})
		df_inter = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutated_sample": [], "Mutation_count": [], "Raw_pvalue": [], "Adjusted_pvalue": []})
	df_all = pd.concat([df_intra, df_inter, df_lof])[["Type", "Structure_source","Uniprots","Residues","Mutated_sample","Mutation_count", "LoF_enrichment", "Raw_pvalue", "Adjusted_pvalue"]]
	df_all["Affected_genes"] = df_all["Uniprots"].apply(lambda x: get_affected_genes(x, uniprot2gene))
	if df_all.shape[0] == 0:
		logging.warning("No selection signature can be found!")
		df_all["Signature_ID"] = np.nan
		df_all["Mutation_frequency"] = np.nan
	else:
		df_all["LoF_enrichment"] = df_all["LoF_enrichment"].fillna("[NA]")
		df_all["Structure_source"] = df_all["Structure_source"].fillna("[NA]")
		df_all["Residues"] = df_all["Residues"].fillna("[NA]")
		df1 = df_all[df_all["Type"] == "LoF_IntraProtein_Enriched"].sort_values(by = "Raw_pvalue")
		if df1.shape[0] > 0:
			df1["Mutation_frequency"] = df1.apply(lambda x: x["Affected_genes"] + ":" + str(sum([int(y) for y in x["Mutation_count"].split(",")])), axis = 1)
			df1_final = []
			for key,df_one in df1.groupby("Affected_genes"):
				df_one = df_one.sort_values(by = "Raw_pvalue").reset_index().drop(columns = ["index"]).reset_index().rename(columns = {"index": "Signature_ID"})
				df_one["Signature_ID"] = df_one.apply(lambda x: x["Affected_genes"] + "_LoF" + str(x["Signature_ID"] + 1), axis = 1)
				df1_final.append(df_one)
			df1 = pd.concat(df1_final).sort_values(by = "Raw_pvalue")
		else:
			df1["Mutation_frequency"] = np.nan
			df1["Signature_ID"] = np.nan
		df2 = df_all[df_all["Type"] != "LoF_IntraProtein_Enriched"]
		if df2.shape[0] > 0:
			df2["Mutation_frequency"] = df2.apply(lambda x: get_mis_mutfreq(x, uniprot2gene), axis = 1)
			df2_final = []
			for key,df_one in df2.groupby("Affected_genes"):
				df_one = df_one.sort_values(by = "Raw_pvalue").reset_index().drop(columns = ["index"]).reset_index().rename(columns = {"index": "Signature_ID"})
				df_one["Signature_ID"] = df_one.apply(lambda x: "_".join(x["Affected_genes"].split(",") + ["NonTrunc" + str(x["Signature_ID"] + 1)]), axis = 1)
				df2_final.append(df_one)
			df2 = pd.concat(df2_final).sort_values(by = "Raw_pvalue")
		else:
			df2["Mutation_frequency"] = np.nan
			df2["Signature_ID"] = np.nan
		df_all = pd.concat([df2, df1])
	df_all["Raw_pvalue"] = df_all["Raw_pvalue"].apply(lambda x: '{:0.2e}'.format(x))
	df_all["Adjusted_pvalue"] = df_all["Adjusted_pvalue"].apply(lambda x: '{:0.2e}'.format(x))
	df_all[["Signature_ID", "Type", "Affected_genes", "Structure_source", "Mutation_frequency", "LoF_enrichment", "Raw_pvalue", "Adjusted_pvalue"]].to_csv(args.output_path + args.job_name + "_signatures.txt", sep = "\t", header = True, index = None)
	############################################################################################################################################
	


	############################################################################################################################################
	# Identify potentially functional mutations
	# -----------------------------------------
	
	# LOF mutations
	df_lof_mut = df_lof_mut[df_lof_mut["UniProt"].isin(lof_uniprots)]
	if df_lof_mut.shape[0] > 0:
		function_input = []
		for x in df_lof_mut["UniProt"]:
			function_input.append([x, df1])
		p = Pool(args.threads)
		df_lof_mut["Signature_ID"] = p.map(get_clusterID_for_a_lof_multirun, function_input)
		p.close()
	else:
		df_lof_mut["Signature_ID"] = np.nan

	# Non-truncating mutations
	uniprot2res = defaultdict(set)
	for x in df2["Residues"]:
		for y in x.split(","):
			uniprot, pos = y.split("_")
			uniprot2res[uniprot] = uniprot2res[uniprot].union({int(pos)})
	if len(uniprot2res) > 0:
		df_mis_mut = df_mis_mut[df_mis_mut["UniProt"].isin(uniprot2res)]
		df_mis_mut = df_mis_mut[df_mis_mut.apply(lambda x: set(funcs.split_consecutive_pos(x["Protein_position"])).intersection(uniprot2res[x["UniProt"]]) != set(), axis = 1)]
		uniprot_list = df_mis_mut["UniProt"].tolist()
		pos_list = df_mis_mut["Protein_position"].tolist()
		function_input = []
		for i in range(len(uniprot_list)):
			function_input.append([uniprot_list[i], pos_list[i], df2])
		p = Pool(args.threads)
		df_mis_mut["Signature_ID"] = p.map(get_clusterID_for_a_mis_multirun, function_input)
		p.close()
	else:
		df_mis_mut["Signature_ID"] = np.nan
	pd.concat([df_mis_mut, df_lof_mut]).sort_values(by = "Signature_ID").to_csv(args.output_path + args.job_name + "_drivers.txt", sep = "\t", header = True, index = None)
	############################################################################################################################################
		
		

	############################################################################################################################################
	# selection signature-oriented network diffusion
	# ----------------------------------------------
	if set(lof_uniprots).intersection(set(G.nodes())) == set() and set(uniprot2res).intersection(set(G.nodes())) == set():
		logging.warning("No available selection signature in the network! Skip the diffusion process!")
		pd.DataFrame({"Subnetwork_genes": [], "Subnetwork_size": []}).to_csv(args.output_path + args.job_name + "_subnetworks.txt", sep = "\t", header = True, index = None)
	else:
		# get heat sources
		ones = df_intra["Uniprots"].unique().tolist()
		pairs = [tuple(item.split(",")) for item in df_inter["Uniprots"].unique()]
		
		# get heat diffusion network
		G_heat_diffusion = funcs.get_initial_distribution(G, pairs, ones, lof_uniprots, args.designated_score, args.enhanced_conductivity, output_path + "initial_state.graphml.gz")
		
		# find appropriate delta
		if len(G_heat_diffusion.edges()) < args.max_subnetwork_size:
			delta = -np.inf
		else:
			f = open(output_path + "choose_delta.txt", "w")
			f.write("\t".join(["delta", "top_edge_cutoff", "subnetwork_sizes"]) + "\n")
			f.close()
			os.mkdir(output_path + "simulation")
			function_input = []
			for i in range(args.delta_trial):
				function_input.append([G_heat_diffusion, output_path + "choose_delta.txt", args.restart_prob, args.max_subnetwork_size, 2062 + i])
			p = Pool(args.threads)
			output = p.map(get_one_delta_multirun, function_input)
			p.close()
			if args.resolution == "low":
				delta = pd.read_csv(output_path + "choose_delta.txt", sep = "\t")["delta"].min()
			else:
				delta = pd.read_csv(output_path + "choose_delta.txt", sep = "\t")["delta"].median()

		# identify subnetworks
		all_subnetworks = funcs.identify_hot_modules(G_heat_diffusion, output_path + "final_state.graphml.gz", beta = args.restart_prob, delta = delta)
		f = open(args.output_path + args.job_name + "_subnetworks.txt", "w")
		f.write("\t".join(["Subnetwork_genes", "Subnetwork_size"]) + "\n")
		for x in sorted(all_subnetworks, key = lambda z: len(z), reverse = True):
			module_id = []
			for y in x:
				module_id.append("/".join(sorted(uniprot2gene[y])))
			f.write("\t".join([",".join(sorted(module_id)), str(len(x))]) + "\n")
		f.close()
		sizes = [len(x) for x in all_subnetworks]
		logging.info(str(len(all_subnetworks)) + " subnetworks are identified, whose sizes range from " + str(min(sizes)) + " to " + str(max(sizes)) + ".")
	df = pd.read_csv(args.output_path + args.job_name + "_subnetworks.txt", sep = "\t")
	df = df[df["Subnetwork_size"] > 1]
	df.to_csv(args.output_path + args.job_name + "_subnetworks.txt", sep = "\t", header = True, index = None)
	df_sig = pd.read_csv(args.output_path + args.job_name + "_signatures.txt", sep = "\t", dtype = str)
	df_sig["Subnetwork_ID"] = df_sig["Affected_genes"].apply(lambda x: get_subnetwork(x, df["Subnetwork_genes"].tolist())).fillna("[NA]")
	index = df_sig["Subnetwork_ID"] == "[NA]"
	df_sig = pd.concat([df_sig[~index], df_sig[index]])
	df_sig.to_csv(args.output_path + args.job_name + "_signatures.txt", sep = "\t", header = True, index = None)
	############################################################################################################################################
	


	############################################################################################################################################
	# remove intermediate files
	# -------------------------
	funcs.remove_whole_dir(output_path)
	############################################################################################################################################


