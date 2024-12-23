import unittest
import filecmp
from io import StringIO
import tempfile
import os
import NetFlow3D


def check_output_file_correctness(source_file, target_string):
	with open(source_file, 'r') as file:
		data = file.read()
	return target_string == data


class TestNetFlow3D(unittest.TestCase):
	
	def test_get_bmr_uniprot(self):
		 # background_mutability_file, id_mapping_file, output, indel_inframe_frac, mutrate_lower_quantile
		 for i, o in [
			 [["metadata/background.txt", "metadata/HUMAN_9606_idmapping.dat.gz", StringIO(), 0.09, 0.01], "example/output/test/mutrate.txt"],
		 ]:
			 with self.subTest(i=i, o=o):
				 NetFlow3D.get_bmr_uniprot(*i)
				 generated_output = i[2].getvalue()
				 self.assertTrue(check_output_file_correctness(o, generated_output))
				
	def test_get_expr_uniprot(self):
		 # id_mapping_file, expr_whitelist_file, expr_input_file, output
		 for i, o in [
			 [["metadata/HUMAN_9606_idmapping.dat.gz", "metadata/expr_white_list.txt", "example/input/expressed_genes.txt", StringIO()], "example/output/test/Expr_uniprots.txt"]
		 ]:
			 with self.subTest(i=i, o=o):
				 NetFlow3D.get_expr_uniprot(*i)
				 generated_output = i[3].getvalue()
				 self.assertTrue(check_output_file_correctness(o, generated_output))
				
	def test_mutation_preprocessing(self):
		 # input_maf, id_mapping_file, output, expr_uniprot_file=None
		 for i, o in [
			 [["example/input/mutations.maf", "metadata/HUMAN_9606_idmapping.dat.gz", StringIO(), "example/output/test/Expr_uniprots.txt"], "example/output/test/Preprocessed.maf"]
		 ]:
			 NetFlow3D.mutation_preprocessing(*i)
			 generated_output = i[2].getvalue()
			 self.assertTrue(check_output_file_correctness(o, generated_output))
			
	def test_lof_analysis(self):
		# preprocessed_maf, mutrate_file, prolen_dict_file, output, lof={"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
		for i, o in [
			[["example/output/test/Preprocessed.maf", "example/output/test/mutrate.txt", "metadata/uniprot2prolen.json", None], "example/output/test/All_intra_LoF_pvalue.txt"]
		]:
			with self.subTest(i=i, o=o):
				i = list(i)
				tempfilename = tempfile.NamedTemporaryFile(delete=False).name
				i[3] = tempfilename
				NetFlow3D.lof_analysis(*i)
				result = filecmp.cmp(tempfilename, o)
				os.unlink(tempfilename)
				self.assertTrue(result)
				
	def test_inframe_analysis(self):
		# preprocessed_maf, mutrate_file, prolen_file, PDB_intra_resource, PDB_inter_resource, AF2_intra_resource, PIONEER_inter_resource, binary_interactome, output_path, intra_dist_upperlimit = 6, inter_dist_upperlimit = 9, threads = 10, inframe = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}
		for i, o in [
			[["example/output/test/Preprocessed.maf", "example/output/test/mutrate.txt", "metadata/uniprot2prolen.json", "graph/PDB_intra/", "graph/PDB_inter/", "graph/AF2_pLDDT0/", "metadata/HomoSapiens_interfaces_PIONEER_veryhigh.txt", "metadata/HomoSapiens_binary_HINThq.txt", None, 6, 9, 1], "example/output/test/"]
		]:
			with self.subTest(i=i, o=o):
				output_dir = tempfile.TemporaryDirectory()
				i = list(i)
				i[8] = output_dir.name + "/"
				NetFlow3D.inframe_analysis(*i)
				result = all(filecmp.cmp(o + f, output_dir.name + "/" + f) for f in ["Residue_intra_pvalue.txt", "UniProt_intra_pvalue.txt", "PDB_intra_pvalue.txt", "PDB_inter_pvalue.txt", "AlphaFold2_intra_pvalue_pLDDT0.txt", "PIONEER_inter_pvalue.txt"])
				output_dir.cleanup()
				self.assertTrue(result)
				
	def test_generate_result_table(self):
		# final_output_intra_lof, final_output_intra_pdb, final_output_inter_pdb, final_output_intra_af2, final_output_inter_pioneer, canonical_isoform_file, output
		for i, o in [
			[["example/output/test/All_intra_LoF_pvalue.txt", "example/output/test/PDB_intra_pvalue.txt", "example/output/test/PDB_inter_pvalue.txt", "example/output/test/AlphaFold2_intra_pvalue_pLDDT0.txt", "example/output/test/PIONEER_inter_pvalue.txt", "metadata/UP000005640_9606.fasta", StringIO()], "example/output/test_signatures.txt"]
		]:
			with self.subTest(i=i, o=o):
				NetFlow3D.generate_result_table(*i)
				generated_output = i[6].getvalue()
				self.assertTrue(check_output_file_correctness(o, generated_output))

	def test_network_propagation(self):
		# binary_interactome, output_file, output_path, expr_uniprot_file=None, final_output_intra_lof=None, final_output_intra_pdb=None, final_output_inter_pdb=None, final_output_intra_af2=None, final_output_inter_pioneer=None, intercept = 1.0, restart_prob = 0.5, max_subnetwork_size = 5, random_seed = 2062, delta_trial = 20, threads = 5
		for i, o in [
			[["metadata/HomoSapiens_binary_HINThq.txt", StringIO(), "example/output/test/", "example/output/test/Expr_uniprots.txt", "example/output/test/All_intra_LoF_pvalue.txt", "example/output/test/PDB_intra_pvalue.txt", "example/output/test/PDB_inter_pvalue.txt", "example/output/test/AlphaFold2_intra_pvalue_pLDDT0.txt", "example/output/test/PIONEER_inter_pvalue.txt", 1.0, 0.5, 5, 2062, 20, 5], "example/output/test_subnetworks.txt"]
		]:
			# with self.subTest(i=i, o=o):
			# 	NetFlow3D.network_propagation(*i)
			# 	generated_output = i[1].getvalue()
			# 	self.assertTrue(check_output_file_correctness(o, generated_output))

			with self.subTest(i=i, o=o):
				i = list(i)
				tempfilename = tempfile.NamedTemporaryFile(delete=False).name
				i[1] = tempfilename
				NetFlow3D.network_propagation(*i)
				result = filecmp.cmp(tempfilename, o)
				os.unlink(tempfilename)
				self.assertTrue(result)
			
unittest.main(argv=[''], verbosity=2, exit=False)
