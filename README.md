# NetFlow3D
NetFlow3D is a computational tool aiming at mapping how somatic mutations act across scales in cancer. If you find NetFlow3D helpful, please cite https://doi.org/10.1101/2023.03.06.531441. You can also upload your data to our web server (http://netflow3d.yulab.org) and run NetFlow3D there.

Python 3.9.10
## Versions the software has been tested on
- Python: 3.9.10
- Linux Distribution: Rocky Linux 9.0 (Blue Onyx)

## Prerequisites
The Python Standard Library and the following packages:
- scipy (version 1.9.3)
- numpy (version 1.23.5)
- networkx (version 2.8.8)
- pandas (version 1.5.2)

## Installation (~5min)

	git clone https://github.com/haiyuan-yu-lab/NetFlow3D.git
	cd NetFlow3D
	python NetFlow3D.py -h

## Run the Unit Tests
	python test_netflow3d.py
	
## Usage
Your command should be in the following format (the contents in `[]` are optional):

	python NetFlow3D.py -m <input_maf> -I <job_name> [-X <expressed_genes>] [-n <binary_interactome>] [-o <output_path>] [-t <threads>]

### Required arguments
- `-m <input_maf>`: replace `<input_maf>` with the path to your MAF file.
- `-I <job_name>`: replace `<job_name>` with a name you preferred for the current job.

### Optional arguments
- `-X <expressed_genes>`: replace `<expressed_genes>` with the path to your file which stores a complete list of expressed genes/proteins (see [Optional input](#optional-input) for how to generate the file). If not specified, all genes/proteins will be considered expressed.
- `-n <binary_interactome>`: replace `<binary_interactome>` with the path to your file which stores a complete list of existing protein-protein interactions (see [Optional input](#optional-input) for how to generate the file). If not specified, NetFlow3D will use the high quality binary interactome of Homo sapiens curated by HINT (http://hint.yulab.org/).
- `-o <output_path>`: replace `<output_path>` with a directory where the output files will be stored. If not specified, the output files will be stored in `./output/`.
- `-t <threads>`: replace `<threads>` with a postive integer. This argument specifies the number of threads to use. If not specified, NetFlow3D will use 5 threads.
	
We provide example input files in `./example/input/`. Here is an example of your command (please run the following command to see if NetFlow3D is working properly, taking ~1min):

	python NetFlow3D.py -m example/input/mutations.maf -I test -X example/input/expressed_genes.txt -n example/input/interactome.txt -t 10
	
If you run the above command, the output should be found in `./output/`, including `test_signatures.txt`, `test_subnetworks_intercept1.0_lowres_edgeweightTrue.txt`, and a folder `test/`. To get an idea of what the output files should look like, please see example output files in `./example/output/`.


## Prepare input files
### Required input
- A Mutation Annotation Format (MAF) file (https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format)

	Required columns:
	<ul>
	<li>Chromosome</li>
	<li>Start_Position</li>
	<li>Variant_Classification</li>
	<li>ENSP</li>
	<li>Transcript_ID</li>
	<li>Gene</li>
	<li>Protein_position</li>
	<li>Tumor_Sample_Barcode</li>
	</ul>

	Other columns can also be present in the MAF file but they will not be used. 

### Optional input
- A text file containing a complete list of genes/transcripts expressed in the contexts where the mutations occur. One ID per line. Ensembl gene ID and Ensembl transcript ID are accepted. Example:

	>ENSG00000163166<br>
	>ENSG00000110422<br>
	>ENSG00000077312<br>
	>ENSG00000180660<br>
	>ENSG00000186635<br>

- A text file containing a complete list of protein-protein interactions existing in the contexts where the mutations occur. One interaction per line. Protein IDs should be separated by tab. Only UniProt ID is accepted. Example:

	>Q9H4A3&emsp;Q9HBL0<br>
	>Q15654&emsp;Q15797<br>
	>P63279&emsp;Q13643<br>
	>O43236&emsp;O43236<br>
	>P01112&emsp;P04049<br>


## Output files
NetFlow3D will output the following and files and a folder. `{job_name}` will be replaced by the job name you specified before. If you run the example command, `{job_name}` will be replaced by `test`. 
- `{job_name}`_signatures.txt

	This a tab-separated file containing the significant 3D clusters and LOF enrichment signals identified by NetFlow3D. The first line is a header. Eight columns are present:
	1. Signature_ID
	2. Type
	3. Uniprots
	4. Canonical_isoform
	6. Structure_source (`[NA]` means not applicable)
	7. Mutation_frequency 

		The content format in this column depends on the content in "Type":
		- If the content in "Type" is “LoF_IntraProtein_Enriched”, the format of this column is `{UniProt ID}:{number of LoF mutations in all samples}`
		- Otherwise, the format of this column is `{residue1}:{number of mutated samples},{residue2}:{number of mutated samples},...`
	8. LoF_enrichment (`[NA]` means not applicable)
	9. Raw_pvalue
	10. Adjusted_pvalue

- `{job_name}`_subnetworks_intercept1.0_lowres_edgeweightTrue.txt

	This is a tab-separated file containing the interconnected modules identified by NetFlow3D. Two columns are present:
	1. Subnetwork_UniProts	
	2. Subnetwork_size (i.e. number of proteins in the interconnected module)

- `{job_name}/`

	This is a folder containing intermediate files by NetFlow3D:
	1. `initial_state_intercept1.0_lowres_edgeweightTrue.graphml.gz`: input to the network propagation model in NetFlow3D.
	2. `final_state_intercept1.0_lowres_edgeweightTrue.graphml.gz`: output from the network propagation model in NetFlow3D.
	3. `choose_delta_intercept1.0_lowres_edgeweightTrue.txt`: δ's from randomized input
	4. `ShortVersion_mutation_data.txt`: mutation information summarized to each residue
	5. `PIONEER_inter_pvalue.txt`, `PDB_intra_pvalue.txt`, `PDB_inter_pvalue.txt`, `AlphaFold2_intra_pvalue_pLDDT0.txt`: 3D cluster information.
	6.  `All_intra_LoF_pvalue.txt`: LOF enrichment information.
	7. `PDB_graph`, `AlphaFold2_graph_pLDDT0`: residue-residue contact map.


