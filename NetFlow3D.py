description = '''
- NetFlow3D pipeline.
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
import sys
from scipy.stats import poisson
file_path = os.path.dirname(os.path.abspath(__file__))


def get_bmr_uniprot(background_mutability_file, id_mapping_file, output, indel_inframe_frac = 0.09, mutrate_lower_quantile = 0.01):
	"""
	Purpose: Calculate the per-basepair in-frame and loss-of-function (LOF) mutation rates for each protein and write the results to a file.

	Parameters
	----------
	background_mutability_file : str
		Path to the file containing background mutability data.
	id_mapping_file : str
		Path to the file containing ID mapping data.
	mutrate_lower_quantile : float
		Quantile threshold for setting minimum mutation rates.
	indel_inframe_frac : float
		Fraction of in-frame indels among all indels.
	output : str
		Path to the output file where the calculated mutation rates will be saved.
	"""
	# Load background mutability data
	mutrate = pd.read_csv(background_mutability_file, sep="\t")
	
	# Calculate per-basepair in-frame and LOF mutation rates for each gene
	mutrate["BMR_Inframe"] = mutrate.apply(
		lambda x: x["BMR_SNV"] * x["Missense_Fraction"] + x["BMR_Indel"] * indel_inframe_frac, axis=1)
	mutrate["BMR_LOF"] = mutrate.apply(
		lambda x: x["BMR_SNV"] * (x["Nonsense_Fraction"] + x["SpliceSite_Fraction"]) + x["BMR_Indel"] * (1 - indel_inframe_frac), axis=1)
	
	# Ensure mutation rates are not below the specified quantile
	for x in ["Inframe", "LOF"]:
		quantile = mutrate["BMR_" + x].quantile(q=mutrate_lower_quantile)
		mutrate.loc[mutrate["BMR_" + x] < quantile, "BMR_" + x] = quantile

	# Load ID mapping data
	df_annot = pd.read_csv(id_mapping_file, sep="\t", dtype=str, header=None)
	df_annot[2] = df_annot[2].apply(lambda x: x.split(".")[0])
	df_annot[0] = df_annot[0].apply(lambda x: x.split("-")[0])
	df_annot = df_annot[[0, 2]].drop_duplicates().rename(columns={0: "UniProt"})

	# Merge mutation rates with UniProt ID mapping
	mutrate = pd.merge(mutrate, df_annot.rename(columns={2: "Hugo_Symbol"}), on="Hugo_Symbol")
	
	# Write the final mutation rates to the output file
	mutrate.to_csv(output, sep="\t", header=True, index=None)
	
	return


def get_expr_uniprot(id_mapping_file, expr_whitelist_file, expr_input_file, output):
	"""
	Processes input expression data to map gene IDs to UniProt IDs and writes the results to a file.

	Parameters
	----------
	id_mapping_file : str
		Path to the file containing ID mapping data.
	output : str
		Path to the output file where the mapped UniProt IDs will be saved.
	expr_whitelist_file : str
		Path to the file containing a whitelist of genes always considered expressed.
	expr_input_file : str, optional
		Path to the file containing the list of expressed genes.
	"""
	try:
		# Load the whitelist of genes always considered expressed
		df_expr_whitelist = pd.read_csv(expr_whitelist_file, sep="\t", header=None, dtype=str)
		
		# Load the list of expressed genes from input file
		df_expr_input = pd.read_csv(expr_input_file, sep="\t", header=None, dtype=str)
		
		# Combine the whitelist and input genes, extracting the base gene ID (before '.' and '-')
		expr_genes = {x.split(".")[0].split("-")[0] for x in pd.concat([df_expr_whitelist, df_expr_input])[0]}
	except Exception:
		sys.exit("Please check the format of input expression data! Each row should contain only one ID. Exiting process.")
	
	# Load ID mapping data
	df_annot = pd.read_csv(id_mapping_file, sep="\t", dtype=str, header=None)
	
	# Process the mapping data to extract base IDs (before '.' and '-')
	df_annot[2] = df_annot[2].apply(lambda x: x.split(".")[0])
	df_annot[0] = df_annot[0].apply(lambda x: x.split("-")[0])
	df_annot = df_annot[[0, 2]].drop_duplicates()

	# Map input gene IDs to UniProt IDs using the processed ID mapping data
	df1 = df_annot[df_annot[0].isin(expr_genes)]
	df2 = df_annot[df_annot[2].isin(expr_genes)]
	df_output = pd.concat([df1, df2])[[0]].drop_duplicates()
	
	# Write the mapped UniProt IDs to the output file
	df_output.to_csv(output, sep="\t", header=None, index=None)
	
	return


def mutation_preprocessing(input_maf, id_mapping_file, output, expr_uniprot_file = None):
	"""
	Preprocesses input MAF file.

	Parameters
	----------
	input_maf : str
		Path to the input MAF file.
	id_mapping_file : str
		Path to the file containing ID mapping data.
	output : str
		Path to the output file where the preprocessed data will be saved.
	expr_uniprot_file : str, optional
		Path to the file containing the list of expressed UniProt IDs.
	"""
	def format_id(ensembl_id):
		"""
		Formats Ensembl ID by removing version numbers.
		"""
		if isinstance(ensembl_id, str):
			return ensembl_id.split(".")[0]
		else:
			return ensembl_id

	# Load MAF file
	try:
		df = pd.read_csv(input_maf, sep="\t", dtype=str)
	except Exception:
		sys.exit("Please check the input MAF file! Exiting process. (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")

	# Check MAF file format
	required_cols = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Tumor_Sample_Barcode", "Protein_position", "ENSP", "Transcript_ID", "Gene"]
	if set(required_cols) - set(df.columns.tolist()) != set():
		logging.error("Please check the input MAF file! The following columns must be included: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Protein_position, Variant_Classification, ENSP, Transcript_ID, Gene. (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
		sys.exit("Exiting process.")
	if set(df["Variant_Classification"].dropna().tolist()) - {"Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "RNA", "Targeted_Region", "Intron", "IGR", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame"} != set():
		logging.error("Please check the content in the column named Variant_Classification in the input MAF file! Exiting process. (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)")
		sys.exit()

	# Load ID mapping data
	df_annot = pd.read_csv(id_mapping_file, sep="\t", dtype=str, header=None)
	df_annot[2] = df_annot[2].apply(lambda x: x.split(".")[0])
	df_annot[0] = df_annot[0].apply(lambda x: x.split("-")[0])
	df_annot = df_annot[[0, 2]].drop_duplicates().rename(columns={0: "UniProt"})
	all_ids = set(df_annot[2])

	# Annotate UniProt IDs
	if "UniProt" in df.columns:
		df = df.drop(columns=["UniProt"])[required_cols].drop_duplicates()  # If "UniProt" column already exists, drop it
	df = df.dropna(subset=["Gene", "Tumor_Sample_Barcode"])
	for col in ["ENSP", "Transcript_ID", "Gene"]:
		df[col] = df[col].apply(format_id)  # Remove version number in Ensembl IDs (e.g. drop ".8" from "ENST00000380152.8")
	ensp_index = df["ENSP"].isin(all_ids)  # If an Ensembl protein ID is available, use it to map to the corresponding UniProt ID.
	df_ensp = pd.merge(df, df_annot.rename(columns={2: "ENSP"}), on="ENSP")
	df = df[~ensp_index]  # If an Ensembl protein ID is not available, try using Ensembl transcript ID to map to the corresponding UniProt ID.
	enst_index = df["Transcript_ID"].isin(all_ids)
	df_enst = pd.merge(df, df_annot.rename(columns={2: "Transcript_ID"}), on="Transcript_ID")
	df = df[~enst_index]  # If an Ensembl transcript ID is not available, try using Ensembl gene ID to map to the corresponding UniProt ID.
	df_ensg = pd.merge(df, df_annot.rename(columns={2: "Gene"}), on="Gene")
	df = pd.concat([df_ensp, df_enst, df_ensg])

	# Exclude mutations mapped to unexpressed UniProt IDs
	if expr_uniprot_file is not None:
		expr_uniprots = set(pd.read_csv(expr_uniprot_file, sep="\t", header=None)[0])
		df = df[df["UniProt"].isin(expr_uniprots)]

	# Write preprocessed file
	if df.shape[0] == 0:
		logging.warning("No mutations remain after preprocessing.")
	df.to_csv(output, sep="\t", header=True, index=None)

	return


def lof_analysis(preprocessed_maf, mutrate_file, prolen_dict_file, output, lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}):
	"""
	Performs statistical tests to compare observed vs. expected loss-of-function (LOF) mutations per protein.

	Parameters
	----------
	preprocessed_maf : str
		Path to the preprocessed MAF file containing mutation data.
	mutrate_file : str
		Path to the file containing mutation rates.
	prolen_dict_file : str
		Path to the JSON file containing protein lengths.
	output : str
		Path to the output file where results will be saved.
	lof : set, optional
		Set of variant classifications considered as LOF mutations.
	"""
	
	# Load protein length (in amino acid residues)
	with open(prolen_dict_file) as f:
		prolen_dict = json.load(f)

	# Load observed LOF mutation data from the preprocessed MAF file
	df_nonsyn = pd.read_csv(preprocessed_maf, sep="\t", dtype=str)
	sample_size = len(df_nonsyn["Tumor_Sample_Barcode"].unique())
	df_lof = df_nonsyn[df_nonsyn["Variant_Classification"].isin(lof) & df_nonsyn["UniProt"].isin(prolen_dict)]
	all_uniprots = set(df_lof["UniProt"])

	# Calculate background LOF mutation rate for each protein
	df_mutrate = pd.read_csv(mutrate_file, sep="\t")
	median = df_mutrate[["Hugo_Symbol", "BMR_LOF"]].drop_duplicates()["BMR_LOF"].median()
	df_mutrate = df_mutrate.groupby("UniProt")
	uniprot2lof_mutrate = {}
	for x in all_uniprots.intersection(set(df_mutrate.groups.keys())):
		# Sum the LOF mutation rates for all genes encoding the same protein
		uniprot2lof_mutrate[x] = df_mutrate.get_group(x)["BMR_LOF"].sum()
	for x in all_uniprots - set(df_mutrate.groups.keys()):
		# For proteins without a known LOF mutation rate, use the median rate of all genes
		uniprot2lof_mutrate[x] = median
	
	if df_lof.shape[0] > 0:
		# Calculate LOF enrichment p-value for each protein
		with open(output, "w") as f:
			f.write("\t".join(["Uniprots", "Mutation_count", "Raw_pvalue"]) + "\n")
			for key, df_one in df_lof.groupby("UniProt"):
				# Expected LOF mutation count = per-basepair LOF mutation rate * coding length * number of WES samples
				expected = uniprot2lof_mutrate[key] * prolen_dict[key] * 3 * sample_size
				# Calculate p-value assuming LOF mutation count follows a Poisson distribution
				pval = poisson.sf(df_one.shape[0] - 1, expected)
				# Log the number of LOF mutations per protein per sample
				mutcount_list = [str(df_onesample.shape[0]) for sample, df_onesample in df_one.groupby("Tumor_Sample_Barcode")]
				f.write("\t".join([key, ",".join(mutcount_list), str(pval)]) + "\n")

		# Multiple test correction
		df = pd.read_csv(output, sep="\t")
		df["Adjusted_pvalue"] = df["Raw_pvalue"].apply(lambda x: min(x * df.shape[0], 1))
		df["Type"] = "LoF_IntraProtein_Enriched"
		df = df[["Type", "Uniprots", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
	else:
		logging.warning("No loss-of-function mutations available after preprocessing!")
		df = pd.DataFrame({x: [] for x in ["Type", "Uniprots", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]})
	
	# Write the final results to the output file
	df.to_csv(output, sep="\t", header=True, index=None)
	
	return


def inframe_analysis(preprocessed_maf, mutrate_file, prolen_file, PDB_intra_resource, PDB_inter_resource, AF2_intra_resource, PIONEER_inter_resource, binary_interactome, output_path, intra_dist_upperlimit = 6, inter_dist_upperlimit = 9, threads = 10, inframe = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}):
	"""
	Performs statistical tests to compare observed vs. expected in-frame mutations per protein, per residue, and per 3D cluster.

	Parameters
	----------
	preprocessed_maf : str
		Path to the preprocessed MAF file containing mutation data.
	mutrate_file : str
		Path to the file containing mutation rates.
	output_path : str
		Path to the directory where output files will be saved.
	inframe : set, optional
		Set of variant classifications considered as in-frame mutations.
	"""

	# Load protein length (in amino acid residues)
	with open(prolen_file) as f:
		prolen_dict = json.load(f)

	# Load observed in-frame mutation data
	df_nonsyn = pd.read_csv(preprocessed_maf, sep="\t", dtype=str)
	sample_size = len(df_nonsyn["Tumor_Sample_Barcode"].unique())
	df_mis = df_nonsyn[df_nonsyn["Variant_Classification"].isin(inframe)].dropna(subset=["Protein_position"])
	df_mis["Protein_position"] = df_mis["Protein_position"].apply(lambda x: x.split("/")[0])
	all_uniprots = set(df_mis["UniProt"])

	# Calculate background in-frame mutation rate for each protein
	df_mutrate = pd.read_csv(mutrate_file, sep="\t")
	median = df_mutrate[["Hugo_Symbol", "BMR_Inframe"]].drop_duplicates()["BMR_Inframe"].median()
	df_mutrate = df_mutrate.groupby("UniProt")
	uniprot2inframe_mutrate = {}
	for x in all_uniprots.intersection(set(df_mutrate.groups.keys())):
		# Sum the in-frame mutation rates for all genes encoding the same protein
		uniprot2inframe_mutrate[x] = df_mutrate.get_group(x)["BMR_Inframe"].sum()
	for x in all_uniprots - set(df_mutrate.groups.keys()):
		# For proteins without a known in-frame mutation rate, use the median rate of all genes
		uniprot2inframe_mutrate[x] = median

	# Define output file paths for different analysis levels
	final_output_intra_res = output_path + "Residue_intra_pvalue.txt"
	final_output_intra_uniprot = output_path + "UniProt_intra_pvalue.txt"
	final_output_intra_pdb = output_path + "PDB_intra_pvalue.txt"
	final_output_inter_pdb = output_path + "PDB_inter_pvalue.txt"
	final_output_intra_af2 = output_path + "AlphaFold2_intra_pvalue_pLDDT0.txt"
	final_output_inter_pioneer = output_path + "PIONEER_inter_pvalue.txt"

	if df_mis.shape[0] == 0:
		logging.warning("No in-frame mutations available after preprocessing!")
	else:
		# Statistical tests per protein
		df_byprot = {"Uniprots": [], "Observed_mut": []}  # Collect the number of observed in-frame mutations per protein
		for prot, df_one in df_mis.groupby("UniProt"):
			df_byprot["Uniprots"].append(prot)
			df_byprot["Observed_mut"].append(df_one[["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"]].drop_duplicates().shape[0])
		df_byprot = pd.DataFrame(df_byprot).sort_values(by = "Uniprots")
		
		# Calculate expected in-frame mutation count
		df_byprot["Expected_mut"] = df_byprot["Uniprots"].apply(lambda x: uniprot2inframe_mutrate[x] * prolen_dict[x] * 3 * sample_size)
		
		# Perform statistical tests assuming a Poisson distribution for mutation counts
		df_byprot["Raw_pvalue"] = df_byprot.apply(lambda row: poisson.sf(row["Observed_mut"] - 1, row["Expected_mut"]), axis=1)
		
		# Adjust p-values for multiple testing
		df_byprot["Adjusted_pvalue"] = df_byprot["Raw_pvalue"].apply(lambda x: min(x * df_byprot.shape[0], 1))
		
		# Log the number of in-frame mutations per protein per sample
		df_byprot["Residues"] = df_byprot.apply(lambda x: x["Uniprots"] + "_NA:" + str(x["Observed_mut"]), axis=1)
		
		# Save per-protein results
		df_byprot.to_csv(final_output_intra_uniprot, sep="\t", header=True, index=None)
		

		# Statistical tests per residue
		# Collect in-frame mutation information for each residue
		funcs.get_res_annotation(df_mis, output_path + "Per_residue_info.txt")
		
		# Create a mapping of UniProt IDs to mutated residue positions
		df = pd.read_csv(output_path + "Per_residue_info.txt", sep="\t")
		uniprot2res = defaultdict(set)
		for uniprot, df_one in df.groupby("UniProt"):
			uniprot2res[uniprot] = set(df_one["UniProt_Position"])
		
		# Add a "Residues" column with a combined UniProt ID and residue position
		df["Residues"] = df.apply(lambda x: f"{x['UniProt']}_{x['UniProt_Position']}", axis=1)
		df["Structure_source"] = "[NA]"
		
		# Save the per-residue mutation data
		df.sort_values(by = "Residues").rename(columns={"UniProt": "Uniprots"})[["Residues", "Structure_source", "Uniprots"]].to_csv(final_output_intra_res, sep="\t", header=True, index=None)
		
		# Calculate p-values for each residue cluster
		funcs.get_cluster_pval(output_path + "Per_residue_info.txt", final_output_intra_res, sample_size, uniprot2inframe_mutrate)
		
		# Load the calculated p-values
		df = pd.read_csv(final_output_intra_res, sep="\t")
		
		# Adjust p-values for multiple testing
		df["Adjusted_pvalue"] = df["Raw_pvalue"].apply(lambda x: min(x * df.shape[0], 1))
		
		# Save the adjusted p-values
		df.to_csv(final_output_intra_res, sep="\t", header=True, index=None)
		

		# Perform statistical tests for 3D clusters from PDB structures
		os.mkdir(output_path + "PDB_graph")

		# Prepare input for multiprocessing
		function_input = []
		for key in uniprot2res:
			function_input.append([
				key, uniprot2res, output_path + "PDB_graph/", PDB_intra_resource, PDB_inter_resource, 
				intra_dist_upperlimit, inter_dist_upperlimit
			])

		# Generate MR-MR (MR: mutated residue) contact graphs using multiprocessing
		with Pool(threads) as p:
			output = p.map(get_graph_by_uniprot_binary_multirun, function_input)

		# Flatten the output list
		function_output = [item for sublist in output for item in sublist]

		if function_output:
			os.mkdir(output_path + "PDB_cluster")

			# Prepare input for cluster generation
			function_input = [
				[x, output_path + "PDB_cluster/" + os.path.split(x)[1].split(".")[0] + ".txt"] 
				for x in set(function_output)
			]

			# Generate clusters using multiprocessing
			with Pool(threads) as p:
				output = p.map(get_cluster_from_graph_multirun, function_input)

			# Separate intra- and inter-protein clusters
			df_intra = []
			df_inter = []
			for x in output:
				df = pd.read_csv(x, sep="\t", header=None)
				uniprot_id = os.path.split(x)[1].split(".")[0].split("_")
				if len(uniprot_id) == 1:
					df_intra.append(df)
				else:
					df = df[df[1] == "inter"]
					df_inter.append(df)

			# Process intra-protein clusters
			if df_intra:
				df_intra = pd.concat(df_intra).drop(columns={1}).rename(columns={0: "Residues"})
				df_intra["Structure_source"] = "PDB"
				df_intra.sort_values(by = "Residues").to_csv(final_output_intra_pdb, sep="\t", header=True, index=None)
				funcs.get_cluster_pval(output_path + "Per_residue_info.txt", final_output_intra_pdb, sample_size, uniprot2inframe_mutrate)

			# Process inter-protein clusters
			if df_inter:
				df_inter = pd.concat(df_inter).drop(columns={1}).rename(columns={0: "Residues"})
				df_inter["Structure_source"] = "PDB"
				df_inter.sort_values(by = "Residues").to_csv(final_output_inter_pdb, sep="\t", header=True, index=None)
				funcs.get_cluster_pval(output_path + "Per_residue_info.txt", final_output_inter_pdb, sample_size, uniprot2inframe_mutrate)
				df_inter = pd.read_csv(final_output_inter_pdb, sep="\t")
				df_inter["Uniprots"] = df_inter["Uniprots"].apply(funcs.binary_interaction)
				df_inter.to_csv(final_output_inter_pdb, sep="\t", header=True, index=None)

			# Clean up temporary directory
			funcs.remove_whole_dir(output_path + "PDB_cluster")


		# 3D clustering algorithm for AlphaFold2 structures
		os.mkdir(output_path + "AlphaFold2_graph_pLDDT0")
		
		# Prepare input for multiprocessing
		function_input = []
		for key in uniprot2res:
			function_input.append([
				key, uniprot2res, output_path + "AlphaFold2_graph_pLDDT0/", 
				AF2_intra_resource, None, 
				intra_dist_upperlimit, None
			])

		# Generate MR-MR contact graphs using multiprocessing
		with Pool(threads) as p:
			output = p.map(get_graph_by_uniprot_binary_multirun, function_input)

		# Flatten the output list
		function_output = [item for sublist in output for item in sublist]

		if function_output:
			# Create directory for AlphaFold2 clusters
			os.mkdir(output_path + "AlphaFold2_cluster")

			# Prepare input for cluster generation
			function_input = [
				[x, output_path + "AlphaFold2_cluster/" + os.path.split(x)[1].split(".")[0] + ".txt"] 
				for x in set(function_output)
			]

			# Generate clusters using multiprocessing
			with Pool(threads) as p:
				output = p.map(get_cluster_from_graph_multirun, function_input)

			# Collect intra-protein clusters
			df_intra = []
			for element in output:
				df_intra.append(pd.read_csv(element, sep="\t", header=None))
			df_intra = pd.concat(df_intra).drop(columns={1}).rename(columns={0: "Residues"})
			df_intra["Structure_source"] = "AlphaFold2"

			# Save intra-protein cluster data
			df_intra.sort_values(by = "Residues").to_csv(final_output_intra_af2, sep="\t", header=True, index=None)
			
			# Calculate p-values for each cluster
			funcs.get_cluster_pval(output_path + "Per_residue_info.txt", final_output_intra_af2, sample_size, uniprot2inframe_mutrate)
			
			# Clean up temporary directory
			funcs.remove_whole_dir(output_path + "AlphaFold2_cluster")

		
		# 3D clustering algorithm for PIONEER PPI interface data
		try:
			# Load the binary interactome data
			df_interactome = pd.read_csv(binary_interactome, sep="\t", header=None, dtype=str)
			
			# Process interactome data to create a list of binary interactions
			binary_interactome = df_interactome.apply(
				lambda x: ",".join(sorted([x[0].split("-")[0], x[1].split("-")[0]])), axis=1
			).tolist()
		except Exception:
			# Exit the process if there is an issue with the interactome file
			sys.exit("Please check the binary interactome file! Each row must contain two UniProt IDs separated by a tab. See example for details. Exiting process.")
		
		# Identify clusters from the PIONEER interface data
		funcs.get_cluster_from_interface(uniprot2res, PIONEER_inter_resource, binary_interactome, final_output_inter_pioneer)
		
		# Calculate p-values for each cluster
		funcs.get_cluster_pval(output_path + "Per_residue_info.txt", final_output_inter_pioneer, sample_size, uniprot2inframe_mutrate)

	# Ensure all output files exist and initialize with empty data if they do not
	for file in [final_output_intra_res, final_output_intra_uniprot, final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2, final_output_inter_pioneer]:
		if not os.path.exists(file):
			# Create an empty DataFrame with the required columns
			df = pd.DataFrame({x: [] for x in ["Residues", "Structure_source", "Uniprots", "Mutation_count", "Observed_mut", "Expected_mut", "Raw_pvalue"]})
			
			# Save the empty DataFrame to the specified file
			df.to_csv(file, sep="\t", header=True, index=None)

	return [final_output_intra_res, final_output_intra_uniprot, final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2, final_output_inter_pioneer]


def generate_result_table(final_output_intra_lof, final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2, final_output_inter_pioneer, canonical_isoform_file, output):
	"""
	Generates a summary table combining data from various sources. The results are saved to a specified output file.
	"""
	
	def get_mis_mutfreq(row):
		"""Generates a string with residues and their corresponding mutation counts."""
		res = row["Residues"].split(",")
		mutcount = row["Mutation_count"].split(",")
		output = [f"{res[i]}:{mutcount[i]}" for i in range(len(res))]
		return ",".join(output)
	
	# Load data from various sources
	df_lof = pd.read_csv(final_output_intra_lof, sep="\t")
	df_structure = funcs.get_sig_cluster_structure(final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2)
	df_interface = pd.read_csv(final_output_inter_pioneer, sep="\t")
	df_interface["Type"] = "InFrame_InterProtein_Cluster"
	df_interface["Adjusted_pvalue"] = df_interface["Raw_pvalue"].apply(lambda x: min(x * df_interface.shape[0], 1)) # Adjust p-values for PIONEER inter-chain clusters using Bonferroni correction
	
	# Combine all data into a single DataFrame
	df_all = pd.concat([df_structure, df_interface, df_lof])[["Type", "Structure_source", "Uniprots", "Residues", "Mutation_count", "Raw_pvalue", "Adjusted_pvalue"]]
	df_all["Structure_source"] = df_all["Structure_source"].fillna("[NA]")
	
	# Process LoF enriched data
	df1 = df_all[df_all["Type"] == "LoF_IntraProtein_Enriched"]
	if not df1.empty:
		df1["Mutation_frequency"] = df1.apply(lambda x: f"{x['Uniprots']}:{sum([int(y) for y in x['Mutation_count'].split(',')])}", axis=1)
		df1["Signature_ID"] = df1["Uniprots"].apply(lambda x: f"{x}_LoF")
		df1 = df1.sort_values(by="Raw_pvalue")
	
	# Process in-frame clusters data
	df2 = df_all[df_all["Type"] != "LoF_IntraProtein_Enriched"]
	if not df2.empty:
		df2["Mutation_frequency"] = df2.apply(get_mis_mutfreq, axis=1)
		df2_final = []
		for key, df_one in df2.groupby("Uniprots"):
			df_one = df_one.sort_values(by="Raw_pvalue").reset_index(drop=True).reset_index().rename(columns={"index": "Signature_ID"})
			df_one["Signature_ID"] = df_one.apply(lambda x: "_".join(x["Uniprots"].split(",") + [f"InFrame{x['Signature_ID'] + 1}"]), axis=1)
			df2_final.append(df_one)
		df2 = pd.concat(df2_final).sort_values(by="Raw_pvalue")
	else:
		df2["Mutation_frequency"] = pd.Series(dtype='object')
		df2["Signature_ID"] = pd.Series(dtype='object')
	
	# Combine LoF and in-frame data
	df_all = pd.concat([df2, df1])
	df_all["Raw_pvalue"] = df_all["Raw_pvalue"].apply(lambda x: '{:0.2e}'.format(x))
	df_all["Adjusted_pvalue"] = df_all["Adjusted_pvalue"].apply(lambda x: '{:0.2e}'.format(x))
	
	# Annotate canonical isoforms
	canonical_isoforms = funcs.extract_uniprot_ids(canonical_isoform_file)
	df_all["Canonical_isoform"] = df_all["Uniprots"].apply(lambda x: funcs.whether_canonical(x, canonical_isoforms))
	
	# Select and order columns for the final output
	df_all = df_all[["Signature_ID", "Type", "Uniprots", "Canonical_isoform", "Structure_source", "Mutation_frequency", "Raw_pvalue", "Adjusted_pvalue"]]
	
	# Save the final summary table to the specified output file
	df_all.to_csv(output, sep="\t", header=True, index=None)

	return


def network_propagation(binary_interactome, output_file, output_path, expr_uniprot_file=None, final_output_intra_lof=None, final_output_intra_pdb=None, final_output_inter_pdb=None, final_output_intra_af2=None, final_output_inter_pioneer=None, intercept = 1.0, restart_prob = 0.5, max_subnetwork_size = 5, random_seed = 2062, delta_trial = 20, threads = 5):
	"""
	Perform 3D-structurally-informed PPI network propagation to identify significantly interconnected modules.

	Parameters
	----------
	output_file : str
		Path to the file where the results will be saved.
	output_path : str
		Directory path for intermediate files and results.
	expr_uniprot_file : str, optional
		Path to the file containing expressed UniProt IDs.
	"""

	# Build PPI network
	try:
		df_interactome = pd.read_csv(binary_interactome, sep="\t", header=None, dtype=str)
		binary_interactome = set(df_interactome.apply(lambda x: (x[0].split("-")[0], x[1].split("-")[0]), axis=1))
		G = nx.Graph(binary_interactome)
	except Exception:
		sys.exit("Please check the binary interactome file! Each row must contain two UniProt IDs separated by a tab. See example for details. Exiting process.")
	
	# Exclude unexpressed proteins from the PPI network
	if expr_uniprot_file:
		expr_uniprots = set(pd.read_csv(expr_uniprot_file, sep="\t", header=None)[0])
		G = nx.Graph(G.subgraph(expr_uniprots))
	
	# Network propagation: Get initial heat distribution
	G_heat_diffusion = funcs.get_initial_distribution(
		G, final_output_intra_pdb, final_output_intra_af2, final_output_inter_pdb, 
		final_output_inter_pioneer, final_output_intra_lof, output_path + "initial_state.graphml.gz", intercept
	)
	
	if sum(nx.get_node_attributes(G_heat_diffusion, 'heat_score').values()) > 0:
		# Determine edge threshold
		if len(G_heat_diffusion.edges()) < max_subnetwork_size:
			delta = -np.inf
		else:
			with open(output_path + "choose_delta.txt", "w") as f:
				f.write("\t".join(["delta", "top_edge_cutoff", "subnetwork_sizes"]) + "\n")
			
			function_input = [
				[G_heat_diffusion, output_path + "choose_delta.txt", restart_prob, max_subnetwork_size, random_seed + i]
				for i in range(delta_trial)
			]
			with Pool(threads) as p:
				output = p.map(get_one_delta_multirun, function_input)
			
			delta = pd.read_csv(output_path + "choose_delta.txt", sep="\t")["delta"].min()

		# Identify interconnected modules
		all_subnetworks = funcs.identify_hot_modules(
			G_heat_diffusion, output_path + "final_state.graphml.gz", beta=restart_prob, delta=delta
		)
		
		# Save the results
		with open(output_file, "w") as f:
			f.write("\t".join(["Subnetwork_UniProts", "Subnetwork_size"]) + "\n")
			for subnetwork in sorted(all_subnetworks, key=lambda x: len(x), reverse=True):
				f.write("\t".join([",".join(sorted(subnetwork)), str(len(subnetwork))]) + "\n")
	else:
		# If no heat sources are found, create empty files
		pd.DataFrame(columns=["delta", "top_edge_cutoff", "subnetwork_sizes"]).to_csv(output_path + "choose_delta.txt", sep="\t", index=False)
		pd.DataFrame(columns=["Subnetwork_UniProts", "Subnetwork_size"]).to_csv(output_file, sep="\t", index=False)

	return

def get_graph_by_uniprot_binary_multirun(args):
	return funcs.get_graph_by_uniprot_binary(*args)

def get_cluster_from_graph_multirun(args):
	return funcs.get_cluster_from_graph(*args)

def get_one_delta_multirun(args):
	return funcs.get_one_delta(*args)


if __name__ == "__main__":
	# User input
	parser = argparse.ArgumentParser(description=description)
	
	# Required arguments
	parser.add_argument('-m', '--input_maf', required=True, type=str, 
						help='Mutation data in MAF format (MAF file format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)')
	parser.add_argument('-I', '--job_name', required=True, type=str, 
						help='Specify a name for your job. The output files will be stored in a folder with this name')
	
	# Optional arguments
	parser.add_argument('-t', '--threads', type=int, default=5, 
						help='Number of threads to use (default: 5)')
	parser.add_argument('-X', '--expr_input_file', type=str, default=None, 
						help='Text file with genes expressed in your context. Each row should contain a gene (HUGO/ENSP/UniProt IDs accepted). Default: all genes considered expressed')
	parser.add_argument('-G', '--expr_whitelist_file', type=str, default=file_path + "/metadata/expr_white_list.txt", 
						help='Text file with genes always considered expressed regardless of RNA-seq data (default: list of 18 well-known cancer genes with low transcript detection levels)')
	parser.add_argument('-n', '--binary_interactome', type=str, default=file_path + "/metadata/HomoSapiens_binary_HINThq.txt", 
						help='Text file with protein-protein interactions. Each row should contain an interaction (IDs separated by tab; HUGO/ENSP/UniProt IDs accepted). Default: human binary interactome from HINT')
	parser.add_argument('-N', '--no_network', action='store_true', 
						help='Do not perform network propagation')
	parser.add_argument('-o', '--output_path', type=str, default=file_path + "/output/", 
						help='Path to the job folder. Default: output folder')
	
	# Native files
	parser.add_argument('-B', '--background_mutability_file', type=str, default=file_path + "/metadata/background.txt", 
						help='Background per-basepair mutation rate of each gene (native file)')
	parser.add_argument('-l', '--prolen_file', type=str, default=file_path + "/metadata/uniprot2prolen.json", 
						help='Sequence length map of human proteins in amino acid residues (native file)')
	parser.add_argument('-a', '--PDB_intra_resource', type=str, default=file_path + "/graph/PDB_intra/", 
						help='Intra-chain residue-residue distances from PDB structures (native resource)')
	parser.add_argument('-e', '--PDB_inter_resource', type=str, default=file_path + "/graph/PDB_inter/", 
						help='Inter-chain residue-residue distances from PDB structures (native resource)')
	parser.add_argument('-d', '--AF2_intra_resource', type=str, default=file_path + "/graph/AF2_pLDDT0/", 
						help='Intra-chain residue-residue distances from AlphaFold structures (native resource)')
	parser.add_argument('-i', '--PIONEER_inter_resource', type=str, default=file_path + "/metadata/HomoSapiens_interfaces_PIONEER_veryhigh.txt", 
						help='Protein-protein interaction interfaces predicted by PIONEER (native file)')
	parser.add_argument('-P', '--id_mapping_file', type=str, default=file_path + "/metadata/HUMAN_9606_idmapping.dat.gz", 
						help='ID conversion file (native file)')
	parser.add_argument('-c', '--canonical_isoform', type=str, default=file_path + "/metadata/UP000005640_9606.fasta", 
						help='Canonical isoforms (native file)')
	
	# Default settings
	parser.add_argument('-q', '--mutrate_lower_quantile', type=float, default=0.01, 
						help='Lower quantile threshold for gene-specific per-basepair mutation rate (default: 0.01)')
	parser.add_argument('-F', '--indel_inframe_frac', type=float, default=0.09, 
						help='Fraction of in-frame indels among all coding indels (default: 0.09)')
	parser.add_argument('-A', '--intra_dist_upperlimit', type=float, default=6.0, 
						help='Intra-chain distance cutoff (default: 6 angstroms)')
	parser.add_argument('-E', '--inter_dist_upperlimit', type=float, default=9.0, 
						help='Inter-chain distance cutoff (default: 9 angstroms)')
	parser.add_argument('-r', '--restart_prob', type=float, default=0.5, 
						help='Restart probability in the diffusion model (default: 0.5)')
	parser.add_argument('-x', '--max_subnetwork_size', type=int, default=5, 
						help='Interconnected module size cutoff for determining appropriate delta (default: 5)')
	parser.add_argument('-D', '--delta_trial', type=int, default=20, 
						help='Number of permutations for determining the value of delta (default: 20)')
	parser.add_argument('-C', '--intercept', type=float, default=1.0, 
						help='Baseline weight for edges (default: 1.0)')
	parser.add_argument('-s', '--random_seed', type=int, default=2062, 
						help='Seed for randomizing heat diffusion process (default: 2062)')
	
	args = parser.parse_args()

	# Check input files
	if not path.exists(args.input_maf):
		sys.exit(f"No such file: {args.input_maf}")
	if args.expr_input_file and not path.exists(args.expr_input_file):
		sys.exit(f"No such file: {args.expr_input_file}")
	if not path.exists(args.expr_whitelist_file):
		sys.exit(f"No such file: {args.expr_whitelist_file}")
	if not path.exists(args.binary_interactome):
		sys.exit(f"No such file: {args.binary_interactome}")
	if not os.path.isdir(args.output_path):
		sys.exit(f"No such directory: {args.output_path}")
	
	# Create job directory
	output_path = os.path.join(args.output_path, args.job_name) + "/"
	try:
		os.mkdir(output_path)
	except FileExistsError:
		sys.exit("The job name is already in use! Please specify another name. Exiting process.")
	
	# Log input parameters
	log_file = os.path.join(output_path, "Input_parameters.txt")
	with open(log_file, 'w') as f:
		f.write("parameter\tvalue\n")
		for arg in vars(args):
			f.write(f"{arg}\t{getattr(args, arg)}\n")
	
	# Run NetFlow3D analysis
	# Calculate per-basepair background mutation rates for each UniProt ID
	get_bmr_uniprot(args.background_mutability_file, args.id_mapping_file, os.path.join(output_path, "mutrate.txt"), args.indel_inframe_frac, args.mutrate_lower_quantile)

	# Generate a list of expressed UniProt IDs if an expression input file is provided; otherwise, consider all UniProt IDs as expressed.
	if args.expr_input_file:
		expr_uniprots = os.path.join(output_path, "Expr_uniprots.txt")
		get_expr_uniprot(args.id_mapping_file, args.expr_whitelist_file, args.expr_input_file, expr_uniprots)
	else:
		expr_uniprots = None

	# Preprocess input mutation data in MAF format
	mutation_preprocessing(args.input_maf, args.id_mapping_file, os.path.join(output_path, "Preprocessed.maf"), expr_uniprots)

	# Perform loss-of-function (LoF) mutation analysis and calculate p-values
	final_output_intra_lof = os.path.join(output_path, "All_intra_LoF_pvalue.txt")
	lof_analysis(os.path.join(output_path, "Preprocessed.maf"), os.path.join(output_path, "mutrate.txt"), args.prolen_file, final_output_intra_lof)

	# Perform in-frame mutation analysis and calculate p-values
	final_output_intra_res, final_output_intra_uniprot, final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2, final_output_inter_pioneer = inframe_analysis(
		os.path.join(output_path, "Preprocessed.maf"), os.path.join(output_path, "mutrate.txt"), args.prolen_file, 
		args.PDB_intra_resource, args.PDB_inter_resource, args.AF2_intra_resource, args.PIONEER_inter_resource, 
		args.binary_interactome, output_path, args.intra_dist_upperlimit, args.inter_dist_upperlimit, args.threads
	)

	# Generate a summary table of mutation signatures
	generate_result_table(final_output_intra_lof, final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2, final_output_inter_pioneer, args.canonical_isoform, os.path.join(args.output_path, f"{args.job_name}_signatures.txt"))

	if args.no_network == False:
		# Perform 3D-structurally-informed PPI network propagation to identify significantly interconnected modules.
		network_propagation(args.binary_interactome, os.path.join(args.output_path, f"{args.job_name}_subnetworks.txt"), 
			output_path, os.path.join(output_path, "Expr_uniprots.txt"), final_output_intra_lof, final_output_intra_pdb, 
			final_output_inter_pdb, final_output_intra_af2, final_output_inter_pioneer, args.intercept, args.restart_prob,
			args.max_subnetwork_size, args.random_seed, args.delta_trial, args.threads)
