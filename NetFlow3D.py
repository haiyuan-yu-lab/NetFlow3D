description = '''
- NetFlow3D: identification of functional modules through structurome-guided network diffusion (currently only the human structurome is available)
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
from itertools import combinations
from scipy.stats import poisson

def eliminate_redundancy(df_one):
	if df_one.shape[0] > 1:
		if df_one.dropna(subset = ["Protein_position"]).shape[0] > 0:
			return df_one.dropna(subset = ["Protein_position"]).sort_values(by = "Protein_position")[0:1]
		else:
			return df_one.sort_values(by = ["Gene"])[0:1]
	else:
		return df_one

def get_mutcount(row, mutrate, flag):
	if row["Hugo_Symbol"] in mutrate.index:
		bmr = mutrate.loc[row["Hugo_Symbol"], "BMR"]
	else:
		bmr = mutrate["BMR"].median()
	if flag == "lof":
		return row["BMR"] * row["Indel_Coverage_Count"] + bmr * (row["Nonsense_Coverage_Count"] + row["SpliceSite_Coverage_Count"])
	else:
		return bmr * row["Missense_Coverage_Count"]
	
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

def format_id(ensembl_id):
	if type(ensembl_id) == str:
		return ensembl_id.split(".")[0]
	else:
		return ensembl_id

def get_mis_mutfreq(row):
	if row["Residues"] != "[NA]":
		res = row["Residues"].split(",")
		mutcount = row["Mutation_count"].split(",")
		output = []
		for i in range(len(res)):
			output.append(res[i] + ":" + mutcount[i])
		return ",".join(output)
	else:
		return "[NA]"

def get_module_for_a_mut(uniprot, uniprot2module):
	if uniprot in uniprot2module:
		size = len(uniprot2module[uniprot].split("_"))
		return [uniprot2module[uniprot], size]
	else:
		return ["[NA]", "[NA]"]

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

def whether_canonical(uniprots, canonical_isoforms):
	output = []
	for x in uniprots.split(","):
		if x in canonical_isoforms:
			output.append("Yes")
		else:
			output.append("No")
	return ",".join(output)

def get_one_delta_multirun(args):
	return funcs.get_one_delta(*args)

file_path = os.path.dirname(os.path.abspath(__file__))
lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
inframe = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}

if __name__ == "__main__":
	# user input
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument('-m','--input_maf', required = True, type = str, help = 'Mutation data in MAF format (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)')
	parser.add_argument('-R','--resolution', required = True, type = str, help = 'high or low. If set to high, smaller subnetworks are favored. If set to low, larger subnetworks are favored.')
	parser.add_argument('-I','--job_name', required = True, type = str, help = 'Please specify a name for your job. The output files will be stored in a folder with this name')
	parser.add_argument('-t','--threads', type = int, default = 5, help = '[OPTIONAL] The number of threads to use. By default 5')
	parser.add_argument('-X','--expressed_genes', type = str, default = None, help = '[OPTIONAL] A text file with the genes expressed in your context. Each row stands for a gene (HUGO/ENSP/UniProt IDs are accepted). By default all genes are considered as expressed')
	parser.add_argument('-n','--binary_interactome', type = str, default = file_path + "/metadata/HomoSapiens_binary_HINThq.txt", help = '[OPTIONAL] A text file with protein-protein interactions. Each row stands for an interaction (IDs should be separated by tab. HUGO/ENSP/UniProt IDs are accepted). Please include all existing interactions in your context. By default we use the human binary interactome from HINT (http://hint.yulab.org/download/HomoSapiens/binary/hq/)')
	parser.add_argument('-o','--output_path', type = str, default = file_path + "/output/", help = '[OPTIONAL] Path to the job folder. If not specified, the job folder will be stored in the default output folder')
	parser.add_argument('-L','--logfile_path', type = str, default = file_path + "/log/", help = '[OPTIONAL] Path to the log file. If not specified, the log file will be stored in the default log folder')
	
	# native files
	parser.add_argument('-B','--background_mutability_file', type = str, default = file_path + "/metadata/background.txt", help = '[OPTIONAL] Background mutation rate of each gene (native file)')
	parser.add_argument('-q','--mutrate_quantile', type = float, default = 0.01, help = '[OPTIONAL] Quantile of the gene-specific mutation rate. By default 0.01.')
	parser.add_argument('-M','--missense_factor', type = float, default = 0.4, help = '[OPTIONAL] Missense factor. By default 0.4.')
	parser.add_argument('-f','--lof_factor', type = float, default = 0.55, help = '[OPTIONAL] LOF factor. By default 0.55.')
	parser.add_argument('-l','--prolen_file', type = str, default = file_path + "/metadata/uniprot2prolen.json", help = '[OPTIONAL] A sequence length map of human proteins (native file)')
	parser.add_argument('-b','--proseq_file', type = str, default = file_path + "/metadata/uniprot2proseq.json", help = '[OPTIONAL] A sequence map of human proteins (native file)')
	parser.add_argument('-i','--ires_file', type = str, default = file_path + "/metadata/HomoSapiens_interfaces_PIONEER_veryhigh.txt", help = '[OPTIONAL] Protein-protein interaction interfaces predicted by PIONEER (native file)')
	parser.add_argument('-a','--PDB_intra_resource', type = str, default = file_path + "/graph/PDB_intra/", help = '[OPTIONAL] Intra-chain residue-residue distances derived from PDB structures (native resource)')
	parser.add_argument('-e','--PDB_inter_resource', type = str, default = file_path + "/graph/PDB_inter/", help = '[OPTIONAL] Inter-chain residue-residue distances derived from PDB structures (native resource)')
	parser.add_argument('-d','--AF2_intra_resource', type = str, default = file_path + "/graph/AF2_pLDDT0/", help = '[OPTIONAL] Intra-chain residue-residue distances derived from AlphaFold structures (native resource)')
	parser.add_argument('-P','--id_mapping', type = str, default = file_path + "/metadata/HUMAN_9606_idmapping.dat.gz", help = '[OPTIONAL] ID conversion file (native file)')
	parser.add_argument('-c','--canonical_isoform', type = str, default = file_path + "/metadata/UP000005640_9606.fasta", help = '[OPTIONAL] Canonical isoforms. (native file)')
	
	# default settings
	parser.add_argument('-S', '--no_structurome', action = 'store_true', help = "[OPTIONAL] Do not leverage the Human Protein Structurome.")
	parser.add_argument('-N','--no_network', action = 'store_true', help = '[OPTIONAL] Do not perform network analysis.')
	parser.add_argument('-s','--significance_level', type = float, default = 0.001, help = '[OPTIONAL] Significance level. By default 0.001')
	parser.add_argument('-A','--intra_dist_upperlimit', type = float, default = 6.0, help = '[OPTIONAL] Intra-chain distance cutoff. By default 6 angstrom')
	parser.add_argument('-E','--inter_dist_upperlimit', type = float, default = 9.0, help = '[OPTIONAL] Inter-chain distance cutoff. By default 9 angstrom')
	parser.add_argument('-r','--restart_prob', type = float, default = 0.5, help = '[OPTIONAL] Restart probability in the diffusion model. By default 0.5')
	parser.add_argument('-x','--max_subnetwork_size', type = int, default = 5, help = '[OPTIONAL] Subnetwork size cutoff used for determining appropriate delta. By default 5')
	parser.add_argument('-D','--delta_trial', type = int, default = 20, help = '[OPTIONAL] Round of permutations when determining the value of delta. By default 20')
	parser.add_argument('-H','--alternative_heat_score', type = str, default = "no_cutoff", help = "[OPTIONAL] Use alternative definition of heat scores. Can be set to no_cutoff, with_cutoff, or none. By default no_cutoff.")
	parser.add_argument('-W', '--no_edge_weight', action = 'store_true', help = "[OPTIONAL] No edge weight in the network propagation.")
	parser.add_argument('-C','--intercept', type = float, default = 1.0, help = "[OPTIONAL] Intercept when defining edge weight. By default 1.")
	parser.add_argument('-w','--enhanced_conductivity', type = float, default = 5.0, help = '[OPTIONAL] The weight of edges with selection signatures. By default 5')
	args = parser.parse_args()


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
		logging.warning("The job name is already in use!")
		sys.exit("The job name is already in use! Please specify another name.")
	if args.logfile_path[-1] == "/":
		logging.basicConfig(filename = args.logfile_path + args.job_name + '.log', filemode = 'w', level = logging.DEBUG, format = '%(levelname)s - %(message)s')
	else:
		logging.basicConfig(filename = args.logfile_path + "/" + args.job_name + '.log', filemode = 'w', level = logging.DEBUG, format = '%(levelname)s - %(message)s')
	###########################################################################################################################################



	############################################################################################################################################
	# Load native data
	# ----------------
	with open(args.prolen_file) as f:
		prolen_dict = json.load(f)
	with open(args.proseq_file) as f:
		proseq_dict = json.load(f)
	df_annot = pd.read_csv(args.id_mapping, sep = "\t", dtype = str, header = None)
	df_annot[2] = df_annot[2].apply(lambda x: x.split(".")[0])
	df_annot[0] = df_annot[0].apply(lambda x: x.split("-")[0])
	df_annot = df_annot[[0, 2]].drop_duplicates().rename(columns = {0: "UniProt"})
	# BMR
	mutrate = pd.read_csv(args.background_mutability_file, sep = "\t")
	mutrate["LOF_Coverage_Count"] = mutrate["Nonsense_Coverage_Count"] + mutrate["SpliceSite_Coverage_Count"]
	mutrate_indel = mutrate[mutrate["Indel_Coverage_Count"] > 0]
	indcov_density = mutrate_indel["Indel_Coverage_Count"].sum() / mutrate_indel["Coding_Length"].sum()
	mutrate_lof = mutrate[mutrate["LOF_Coverage_Count"] > 0]
	lofcov_density = mutrate_lof["LOF_Coverage_Count"].sum() / mutrate_lof["Coding_Length"].sum()
	mutrate_mis = mutrate[mutrate["Missense_Coverage_Count"] > 0]
	miscov_density = mutrate_mis["Missense_Coverage_Count"].sum() / mutrate_mis["Coding_Length"].sum()
	# gene-specific BMR
	mutrate = mutrate[mutrate["Background_Mutation_Count"] < mutrate["Background_Coverage_Count"]]
	mutrate["BMR"] = mutrate["Background_Mutation_Count"] / mutrate["Background_Coverage_Count"]
	quantile = mutrate["BMR"].quantile(q = args.mutrate_quantile)
	mutrate.loc[mutrate["BMR"] < quantile, "BMR"] = quantile
	mutrate = mutrate.set_index("Hugo_Symbol")
	mis_factor_for_unknown_genes = mutrate["BMR"].median() * miscov_density
	mutrate_mis["Missense_Mutation_Count"] = mutrate_mis.apply(lambda x: get_mutcount(x, mutrate, "mis"), axis = 1)
	mutrate_mis = pd.merge(mutrate_mis, df_annot.rename(columns = {2: "Hugo_Symbol"}), on = "Hugo_Symbol")
	# gene-specific BMR (indel)
	mutrate_indel["BMR"] = mutrate_indel["Background_Mutation_Count_Indel"] / mutrate_indel["Background_Coverage_Count_Indel"]
	quantile = mutrate_indel["BMR"].quantile(q = args.mutrate_quantile)
	mutrate_indel.loc[mutrate_indel["BMR"] < quantile, "BMR"] = quantile
	mutrate_indel["LOF_Mutation_Count"] = mutrate_indel.apply(lambda x: get_mutcount(x, mutrate, "lof"), axis = 1)
	lof_factor_for_unknown_genes = mutrate_indel["BMR"].median() * indcov_density + mutrate["BMR"].median() * lofcov_density
	mutrate_indel = pd.merge(mutrate_indel, df_annot.rename(columns = {2: "Hugo_Symbol"}), on = "Hugo_Symbol")
	# adjustment
	mutrate_mis["Missense_Mutation_Count"] = mutrate_mis["Missense_Mutation_Count"] * args.missense_factor
	mis_factor_for_unknown_genes = mis_factor_for_unknown_genes * args.missense_factor
	mutrate_indel["LOF_Mutation_Count"] = mutrate_indel["LOF_Mutation_Count"] * args.lof_factor
	lof_factor_for_unknown_genes = lof_factor_for_unknown_genes * args.lof_factor
	############################################################################################################################################



	############################################################################################################################################
	# Prepare mutation data
	# ---------------------
	try:
		df = pd.read_csv(args.input_maf, sep = "\t", dtype = str)
	except:
		funcs.remove_whole_dir(output_path)
		sys.exit("Please check the input MAF file! (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
	if {"Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Tumor_Sample_Barcode", "Protein_position", "ENSP", "Transcript_ID", "Gene"} - set(df.columns.tolist()) != set():
		logging.error("Please check the input MAF file! The following columns must be included: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Protein_position, Variant_Classification, ENSP, Transcript_ID, Gene. (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
		funcs.remove_whole_dir(output_path)
		sys.exit()
	if set(df["Variant_Classification"].dropna().tolist()) - {"Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "RNA", "Targeted_Region", "Intron", "IGR", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame"} != set():
		logging.error("Please check the content in the column named Variant_Classification in the input MAF file! (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
		funcs.remove_whole_dir(output_path)
		sys.exit()
	df = df.dropna(subset = ["Tumor_Sample_Barcode", "Gene"])
	for col in ["ENSP", "Transcript_ID", "Gene"]:
		df[col] = df[col].apply(format_id)
	if df.shape[0] == 0:
		logging.error("Please check the columns 'Tumor_Sample_Barcode' and 'Gene'!")
		funcs.remove_whole_dir(output_path)
		sys.exit()
	
	# Annotate UniProt IDs
	if "UniProt" in df.columns:
		df = df.drop(columns = ["UniProt"])
	all_ids = set(df_annot[2])
	ensp_index = df["ENSP"].isin(all_ids)
	df_ensp = pd.merge(df, df_annot.rename(columns = {2: "ENSP"}), on = "ENSP")
	df = df[~ensp_index]
	enst_index = df["Transcript_ID"].isin(all_ids)
	df_enst = pd.merge(df, df_annot.rename(columns = {2: "Transcript_ID"}), on = "Transcript_ID")
	df = df[~enst_index]
	df_ensg = pd.merge(df, df_annot.rename(columns = {2: "Gene"}), on = "Gene")
	df = pd.concat([df_ensp, df_enst, df_ensg])
	df = df[df["UniProt"].isin(prolen_dict)]


	############################################################################################################################################
	if df.shape[0] == 0:
		logging.warning("Could not map the mutations to any protein!")
		funcs.remove_whole_dir(output_path)
		pd.DataFrame({"Subnetwork_UniProts": [], "Subnetwork_size": []}).to_csv(args.output_path + args.job_name + "_subnetworks.txt", sep = "\t", header = True, index = None)
		pd.DataFrame({"Signature_ID": [], "Type": [], "Uniprots": [], "Structure_source": [], "Mutation_frequency": [], "Raw_pvalue": [], "Adjusted_pvalue": []}).to_csv(args.output_path + args.job_name + "_signatures.txt", sep = "\t", header = True, index = None)
		pd.DataFrame({"Chromosome": [], "Start_Position": [], "End_Position": [], "Reference_Allele": [], "Tumor_Seq_Allele1": [], "Tumor_Seq_Allele2": [], "Variant_Classification": [], "Codons": [], "ENSP": [], "Protein_position": [], "Tumor_Sample_Barcode": [], "UniProt": [], "Signature_ID": []}).to_csv(args.output_path + args.job_name + "_drivers.txt", sep = "\t", header = True, index = None)
		sys.exit()
	############################################################################################################################################
	
	

	############################################################################################################################################
	# Build protein-protein interaction network
	# -----------------------------------------
	try:
		df_interactome = pd.read_csv(args.binary_interactome, sep = "\t", header = None, dtype = str)
		df_interactome[2] = df_interactome.apply(lambda x: tuple(sorted([x[0].split("-")[0], x[1].split("-")[0]])), axis = 1)
		G = nx.Graph(set(df_interactome[2]))
	except:
		funcs.remove_whole_dir(output_path)
		sys.exit("Please check the file with binary interactome! Each row must contain two UniProt IDs separated by tab. See example for details.")
	############################################################################################################################################

	
	
	############################################################################################################################################
	# Filter out mutations in unexpressed genes and remove unexpressed proteins from the network
	# ------------------------------------------------------------------------------------------
	if args.expressed_genes != None:
		try:
			expr_genes = {x.split(".")[0].split("-")[0] for x in pd.read_csv(args.expressed_genes, sep = "\t", header = None, dtype = str)[0]}
		except:
			funcs.remove_whole_dir(output_path)
			sys.exit("Please check the format of the file with expressed genes! Each row should only contain one Ensembl gene/transcript ID. See example for details.")
		df1 = df_annot[df_annot["UniProt"].isin(expr_genes)]
		df2 = df_annot[df_annot[2].isin(expr_genes)]
		
		# Proteins encoded by the 18 well-known cancer genes with low transcript detection levels
		# all isoforms
		expr_uniprots = {"P10275", "A0A024RBL6","A0A024RBQ4","A0A087WX99","A0A087WY21","A0A087X2D1","A0A0A0MSE1","A0A0C3SFZ7","A0A140VJJ0","A0A141AXF1","A0A1B0GWF4","A0A1B0GXI8","A0A1W2PQT9","A0A1W2PR07","A0A494C0I9","A0A5S6RJB7","A0A804HI10","A0A804HI76","A0A804HIH8","A0A8V8TPW8","B0ZTD4","B7Z2I3","C7FEN9","C9JU02","C9JXA2","C9JYS6","D2CGD1","D6RDX0","D6RG11","D6RIG5","D6RJH0","E7BSV0","E7EPY2","E7ER61","E7ERX0","E7EU48","E7EUL6","E7EVR7","E9PDR1","E9PGE9","F2YGG7","G3V4B9","H0Y3F0","H0Y7K5","H0YED9","H3BLT0","H7BXU9","H7C265","H9KVD4","J3KNN9","O00570","O15119","O15457","P00533","P16234","P19544","P21802","P29320","P36888","P48436","P55283","P55317","Q03112","Q15303","Q504U8","Q6P4R6","Q6PI38","Q9H6I2","Q9Y261","S4R381","S4R3B2"}
		# # canonical isoform
		# expr_uniprots = {'O00570', 'O15119', 'O15457', 'P00533', 'P10275', 'P16234', 'P19544', 'P21802', 'P29320', 'P36888', 'P48436', 'P55283', 'P55317', 'Q03112', 'Q15303', 'Q9H6I2', 'Q9Y261'}
		expr_uniprots = expr_uniprots.union(set(df1["UniProt"].tolist() + df2["UniProt"].tolist()))
		G = nx.Graph(G.subgraph(expr_uniprots))
		df = df[df["UniProt"].isin(expr_uniprots)]
		if df.shape[0] == 0:
			logging.warning("No mutation is in expressed genes!")
			funcs.remove_whole_dir(output_path)
			pd.DataFrame({"Subnetwork_UniProts": [], "Subnetwork_size": []}).to_csv(args.output_path + args.job_name + "_subnetworks.txt", sep = "\t", header = True, index = None)
			pd.DataFrame({"Signature_ID": [], "Type": [], "Uniprots": [], "Structure_source": [], "Mutation_frequency": [], "Raw_pvalue": [], "Adjusted_pvalue": []}).to_csv(args.output_path + args.job_name + "_signatures.txt", sep = "\t", header = True, index = None)
			pd.DataFrame({"Chromosome": [], "Start_Position": [], "End_Position": [], "Reference_Allele": [], "Tumor_Seq_Allele1": [], "Tumor_Seq_Allele2": [], "Variant_Classification": [], "Codons": [], "ENSP": [], "Protein_position": [], "Tumor_Sample_Barcode": [], "UniProt": [], "Signature_ID": []}).to_csv(args.output_path + args.job_name + "_drivers.txt", sep = "\t", header = True, index = None)
			sys.exit()
	############################################################################################################################################



	# ############################################################################################################################################
	# # Sanity check
	# df_inframe = df[df["Variant_Classification"].isin(inframe)].dropna(subset = ["Protein_position"])
	# df_lof = df[df["Variant_Classification"].isin(lof)]
	# function_input = []
	# for x in ["df_inframe", "df_lof"]:
	# 	for key, df_one in eval(x).groupby(["Chromosome", "Start_Position", "Tumor_Sample_Barcode", "UniProt"]):
	# 		function_input.append(df_one)
	# p = Pool(args.threads)
	# df = p.map(eliminate_redundancy, function_input)
	# p.close()
	# df = pd.concat(df)

	# # Write preprocessed file
	# df.to_csv(output_path + "Preprocessed.maf", sep = "\t", header = True, index = None)
	# ############################################################################################################################################



	############################################################################################################################################
	# Map genes to UniProts
	# ---------------------
	all_uniprots = set(df["UniProt"])
	n_lof = len(all_uniprots)
	n = len(df["Tumor_Sample_Barcode"].unique())
	uniprot2mis_mutcount = {}
	for x,df_one in mutrate_mis.groupby("UniProt"):
		uniprot2mis_mutcount[x] = df_one["Missense_Mutation_Count"].sum()
	uniprot2lof_mutcount = {}
	for x,df_one in mutrate_indel.groupby("UniProt"):
		uniprot2lof_mutcount[x] = df_one["LOF_Mutation_Count"].sum()
	############################################################################################################################################


	
	############################################################################################################################################
	# Identify selection signatures formed by LOF mutations
	# -----------------------------------------------------
	df_lof = df[df["Variant_Classification"].isin(lof)]
	df_lof_mut = copy.deepcopy(df_lof)
	if df_lof.shape[0] > 0:
		# protein-specific enrichment analysis
		logging.info(str(df_lof.shape[0]) + " available loss-of-function mutation(s).")
		# --------------------------
		# total_mut_expected = 0
		# for x in all_uniprots:
		# 	total_mut_expected = total_mut_expected + uniprot2mutability[x] * prolen_dict[x] * 3 * n * frac_lof
		# total_mut_observed = df_lof.shape[0]
		# --------------------------
		final_output_intra_lof = output_path + "All_intra_LoF_pvalue.txt"
		f = open(final_output_intra_lof, "w")
		f.write("\t".join(["Uniprots", "Mutated_sample", "Mutation_count", "Raw_pvalue"]) + "\n")
		for key, df_one in df_lof.groupby("UniProt"):
			if key in uniprot2lof_mutcount:
				expected = uniprot2lof_mutcount[key] * n
			else:
				expected = lof_factor_for_unknown_genes * n * prolen_dict[key] * 3 
			pval = poisson.sf(df_one.shape[0] - 1, expected)
			# enr, _, pval = funcs.calc_enr(df_one.shape[0], total_mut_observed, uniprot2mutability[key] * prolen_dict[key] * 3 * n * frac_lof, total_mut_expected)
			sample_list = []
			mutcount_list = []
			for sample, df_onesample in df_one.groupby("Tumor_Sample_Barcode"):
				sample_list.append(sample)
				mutcount_list.append(str(df_onesample.shape[0]))
			f.write("\t".join([key, ",".join(sample_list), ",".join(mutcount_list), str(pval)]) + "\n")
		f.close()

		# multiple test correction
		df_lof = pd.read_csv(final_output_intra_lof, sep = "\t")
		df_lof["Adjusted_pvalue"] = df_lof["Raw_pvalue"].apply(lambda x: min(x * n_lof, 1))
		df_lof["Type"] = "LoF_IntraProtein_Enriched"
		df_lof = df_lof[["Type", "Uniprots", "Mutated_sample", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
		lof_uniprots = df_lof[df_lof["Adjusted_pvalue"] < args.significance_level]["Uniprots"].tolist()
		logging.info(str(df_lof.shape[0]) + " proteins are significantly enriched of loss-of-function mutations.")
	else:
		logging.warning("No available loss-of-function mutation!")
		lof_uniprots = []
		df_lof = pd.DataFrame({"Type": [], "Uniprots": [], "Mutated_sample": [], "Mutation_count": [], "Raw_pvalue": [], "Adjusted_pvalue": []})
	############################################################################################################################################
	


	############################################################################################################################################
	# Identify selection signatures formed by in-frame mutations
	# ----------------------------------------------------------------
	df_mis = df[df["Variant_Classification"].isin(inframe)].dropna(subset = ["Protein_position", "Codons"])
	df_mis = df_mis[df_mis["Codons"].apply(lambda x: set(x.split("/")[0].upper()) - {"A", "T", "C", "G"} == set())]
	df_mis["Protein_position"] = df_mis["Protein_position"].apply(lambda x: x.split("/")[0])
	df_mis_mut = copy.deepcopy(df_mis)
	if df_mis.shape[0] == 0:
		logging.warning("No available in-frame mutation!")
		skip_flag = 1
	else:
		shortversion = funcs.get_res_annotation(df_mis, output_path)
		df = pd.read_csv(shortversion, sep = "\t", dtype = {"Uniprot_Position": int})
		if df.shape[0] == 0:
			logging.warning("No available in-frame mutation!")
			skip_flag = 1
		else:
			logging.info(str(df.shape[0]) + " available in-frame mutation(s)!")
			skip_flag = 0

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
			# identify spatial clusters and evaluate their significance
			uniprot2res = defaultdict(set)
			for x,df_one in df.groupby("Uniprot"):
				uniprot2res[x] = set(df_one["Uniprot_Position"])
			if args.no_structurome:
				clusters = output_path + "Intra_protein_clusters.txt"
				f = open(clusters, "w")
				for uniprot in uniprot2res:
					for res in uniprot2res[uniprot]:
						f.write("\t".join([uniprot + "_" + str(res), "[NA]"]) + "\n")
				f.close()
				final_output = funcs.get_cluster_pval(shortversion, clusters, n, output_path + "Intra_pvalue.txt", uniprot2mis_mutcount, prolen_dict, mis_factor_for_unknown_genes)
				df_structure = funcs.get_sig_cluster_structure(final_output, None, None, all_uniprots, prolen_dict)
				df_interface = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutation_count": [], "Raw_pvalue": [], "Adjusted_pvalue": []})
			else:
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
						final_output_intra_pdb = funcs.get_cluster_pval(shortversion, output_path + "Intra_protein_clusters_PDB.txt", n, output_path + "PDB_intra_pvalue.txt", uniprot2mis_mutcount, prolen_dict, mis_factor_for_unknown_genes)
					else:
						final_output_intra_pdb = None
					if len(df_inter) > 0:
						df_inter = pd.concat(df_inter)
						df_inter[1] = "PDB"
						df_inter.to_csv(output_path + "Inter_protein_clusters_PDB.txt", sep = "\t", header = None, index = None)
						final_output_inter_pdb = funcs.get_cluster_pval(shortversion, output_path + "Inter_protein_clusters_PDB.txt", n, output_path + "PDB_inter_pvalue.txt", uniprot2mis_mutcount, prolen_dict, mis_factor_for_unknown_genes)
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
				# p = Pool(80)
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
					# p = Pool(80)
					output = p.map(get_cluster_from_graph_multirun, function_input)
					p.close()
					df_intra = []
					for element in output:
						df_intra.append(pd.read_csv(element, sep = "\t", header = None))
					df_intra = pd.concat(df_intra)
					df_intra[1] = "AlphaFold2"
					df_intra.to_csv(output_path + "Intra_protein_clusters_AlphaFold2.txt", sep = "\t", header = None, index = None)
					final_output_intra_af2 = funcs.get_cluster_pval(shortversion, output_path + "Intra_protein_clusters_AlphaFold2.txt", n, output_path + "AlphaFold2_intra_pvalue_pLDDT0.txt", uniprot2mis_mutcount, prolen_dict, mis_factor_for_unknown_genes)
				else:
					final_output_intra_af2 = None

				# PIONEER
				INSIDER_clusters = funcs.get_cluster_from_interface(uniprot2res, args.ires_file, [",".join(sorted([u,v])) for u,v in G.edges()], output_path + "PIONEER_clusters.txt")
				if INSIDER_clusters != None:
					final_output_inter_pioneer = funcs.get_cluster_pval(shortversion, INSIDER_clusters, n, output_path + "PIONEER_inter_pvalue.txt", uniprot2mis_mutcount, prolen_dict, mis_factor_for_unknown_genes)
				else:
					final_output_inter_pioneer = None

				df_structure = funcs.get_sig_cluster_structure(final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2, all_uniprots, prolen_dict)
				df_interface = funcs.get_sig_cluster_interface(final_output_inter_pioneer)
	if skip_flag == 1:
		df_structure = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutation_count": [], "Raw_pvalue": [], "Adjusted_pvalue": []})
		df_interface = pd.DataFrame({"Type": [], "Structure_source": [], "Uniprots": [], "Residues": [], "Mutation_count": [], "Raw_pvalue": [], "Adjusted_pvalue": []})
	df_all = pd.concat([df_structure, df_interface, df_lof])[["Type", "Structure_source", "Uniprots", "Residues", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
	if df_all.shape[0] == 0:
		logging.warning("No selection signature can be found!")
		df_all["Signature_ID"] = np.nan
		df_all["Mutation_frequency"] = np.nan
	else:
		df_all["Structure_source"] = df_all["Structure_source"].fillna("[NA]")
		df_all["Residues"] = df_all["Residues"].fillna("[NA]")
		df1 = df_all[df_all["Type"] == "LoF_IntraProtein_Enriched"]
		if df1.shape[0] > 0:
			df1["Mutation_frequency"] = df1.apply(lambda x: x["Uniprots"] + ":" + str(sum([int(y) for y in x["Mutation_count"].split(",")])), axis = 1)
			df1["Signature_ID"] = df1["Uniprots"].apply(lambda x: x + "_LoF")
			df1 = df1.sort_values(by = "Raw_pvalue")
		else:
			df1["Mutation_frequency"] = np.nan
			df1["Signature_ID"] = np.nan
		df2 = df_all[df_all["Type"] != "LoF_IntraProtein_Enriched"]
		if df2.shape[0] > 0:
			df2["Mutation_frequency"] = df2.apply(get_mis_mutfreq, axis = 1)
			df2_final = []
			for key,df_one in df2.groupby("Uniprots"):
				df_one = df_one.sort_values(by = "Raw_pvalue").reset_index().drop(columns = ["index"]).reset_index().rename(columns = {"index": "Signature_ID"})
				df_one["Signature_ID"] = df_one.apply(lambda x: "_".join(x["Uniprots"].split(",") + ["InFrame" + str(x["Signature_ID"] + 1)]), axis = 1)
				df2_final.append(df_one)
			df2 = pd.concat(df2_final).sort_values(by = "Raw_pvalue")
		else:
			df2["Mutation_frequency"] = np.nan
			df2["Signature_ID"] = np.nan
		df_all = pd.concat([df2, df1])
	df_all["Raw_pvalue"] = df_all["Raw_pvalue"].apply(lambda x: '{:0.2e}'.format(x))
	df_all["Adjusted_pvalue"] = df_all["Adjusted_pvalue"].apply(lambda x: '{:0.2e}'.format(x))
	canonical_isoforms = funcs.extract_uniprot_ids(args.canonical_isoform)
	df_all["Canonical_isoform"] = df_all["Uniprots"].apply(lambda x: whether_canonical(x, canonical_isoforms))
	df_all[["Signature_ID", "Type", "Uniprots", "Canonical_isoform", "Structure_source", "Mutation_frequency", "Raw_pvalue", "Adjusted_pvalue"]].to_csv(args.output_path + args.job_name + "_signatures.txt", sep = "\t", header = True, index = None)
	############################################################################################################################################
	
		

	############################################################################################################################################
	# Network diffusion
	# ----------------------------------------------	
	if args.no_network == False:
		if all_uniprots.intersection(set(G.nodes())) == set():
			logging.warning("Skip the diffusion process!")
			pd.DataFrame({"Subnetwork_UniProts": [], "Subnetwork_size": []}).to_csv(args.output_path + args.job_name + "_subnetworks.txt", sep = "\t", header = True, index = None)
		else:
			final_output_intra_pdb = output_path + "PDB_intra_pvalue.txt"
			final_output_inter_pdb = output_path + "PDB_inter_pvalue.txt"
			final_output_intra_af2 = output_path + "AlphaFold2_intra_pvalue_pLDDT0.txt"
			final_output_inter_pioneer = output_path + "PIONEER_inter_pvalue.txt"
			final_output_intra_lof = output_path + "All_intra_LoF_pvalue.txt"
			if args.alternative_heat_score == "no_cutoff":
				tag = "intercept" + str(args.intercept) + "_" + args.resolution + "res_edgeweight" + str(not args.no_edge_weight)
				G_heat_diffusion = funcs.get_initial_distribution_alternate(G, final_output_intra_pdb, final_output_intra_af2, final_output_inter_pdb, final_output_inter_pioneer, final_output_intra_lof, output_path + "initial_state_" + tag + ".graphml.gz", args.intercept, not args.no_edge_weight)
			elif args.alternative_heat_score == "with_cutoff":
				tag = "intercept" + str(args.intercept) + "_" + args.resolution + "res_sig" + str(args.significance_level)
				df_intra = df_intra[df_intra["Adjusted_pvalue"] < args.significance_level]
				df_inter = df_inter[df_inter["Adjusted_pvalue"] < args.significance_level]
				df_lof = df_lof[df_lof["Adjusted_pvalue"] < args.significance_level]
				G_heat_diffusion = funcs.get_initial_distribution_alternate_sig(G, df_intra, df_inter, df_lof, args.intercept, output_path + "initial_state_" + tag + ".graphml.gz")
			else:
				tag = "weight" + str(args.enhanced_conductivity) + "_" + args.resolution + "res_sig" + str(args.significance_level)
				# get heat sources
				df_intra = df_intra[df_intra["Adjusted_pvalue"] < args.significance_level]
				df_inter = df_inter[df_inter["Adjusted_pvalue"] < args.significance_level]
				pairs = [tuple(item.split(",")) for item in df_inter["Uniprots"].unique()]
				ones = df_intra["Uniprots"].unique().tolist()
				# get heat diffusion network
				G_heat_diffusion = funcs.get_initial_distribution(G, pairs, ones, lof_uniprots, 1.0, args.enhanced_conductivity, output_path + "initial_state_" + tag + ".graphml.gz")
		
			# find appropriate delta
			if len(G_heat_diffusion.edges()) < args.max_subnetwork_size:
				delta = -np.inf
			else:
				f = open(output_path + "choose_delta_" + tag + ".txt", "w")
				f.write("\t".join(["delta", "top_edge_cutoff", "subnetwork_sizes"]) + "\n")
				f.close()
				function_input = []
				for i in range(args.delta_trial):
					function_input.append([G_heat_diffusion, output_path + "choose_delta_" + tag + ".txt", args.restart_prob, args.max_subnetwork_size, 2062 + i])
				p = Pool(args.threads)
				output = p.map(get_one_delta_multirun, function_input)
				p.close()
				if args.resolution == "low":
					delta = pd.read_csv(output_path + "choose_delta_" + tag + ".txt", sep = "\t")["delta"].min()
				else:
					delta = pd.read_csv(output_path + "choose_delta_" + tag + ".txt", sep = "\t")["delta"].median()

			# identify subnetworks
			all_subnetworks = funcs.identify_hot_modules(G_heat_diffusion, output_path + "final_state_" + tag + ".graphml.gz", beta = args.restart_prob, delta = delta)
			f = open(args.output_path + args.job_name + "_subnetworks_" + tag + ".txt", "w")
			f.write("\t".join(["Subnetwork_UniProts", "Subnetwork_size"]) + "\n")
			for x in sorted(all_subnetworks, key = lambda z: len(z), reverse = True):
				f.write("\t".join([",".join(sorted(x)), str(len(x))]) + "\n")
			f.close()
			sizes = [len(x) for x in all_subnetworks]
			logging.info(str(len(all_subnetworks)) + " subnetworks are identified, whose sizes range from " + str(min(sizes)) + " to " + str(max(sizes)) + ".")
	############################################################################################################################################
	


	# ############################################################################################################################################
	# # remove intermediate files
	# # -------------------------
	# funcs.remove_whole_dir(output_path)
	# ############################################################################################################################################


