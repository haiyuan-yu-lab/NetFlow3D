# NetFlow3D
NetFlow3D is a computational tool aiming at mapping how somatic mutations act across scales in cancer. If you find NetFlow3D helpful, please cite https://doi.org/10.1101/2023.03.06.531441. You can also upload your data to our web server (http://netflow3d.yulab.org) and run NetFlow3D there.

## Prerequisites
The Python Standard Library and the following packages:
- scipy (version 1.9.3)
- numpy (version 1.23.5)
- networkx (version 2.8.8)
- pandas (version 1.5.2)
- statsmodels (version 0.13.5)

## Installation

	git clone https://github.com/hyulab/NetFlow3D.git
	cd NetFlow3D
	python NetFlow3D.py -h
	
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
	<li>Strand</li>
	<li>Reference_Allele</li>
	<li>Tumor_Seq_Allele1</li>
 	<li>Tumor_Seq_Allele2</li>
	<li>Variant_Classification</li>
	<li>ENSP</li>
	<li>Transcript_ID</li>
	<li>Gene</li>
	<li>Protein_position</li>
	<li>Codons</li>
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
NetFlow3D will output the following files. `{job_name}` will be replaced by the job name you specified before. If you run the example command, `{job_name}` will be replaced by `test`. 
- `{job_name}`_signatures.txt

	This a tab-separated file containing the significant 3D clusters and LOF enrichment signals identified by NetFlow3D. The first line is a header. Eight columns are present:
	1. Signature_ID
	2. Type
	3. Affected_genes
	4. Structure_source (`[NA]` means not applicable)
	5. Mutation_frequency 

		The content format in this column depends on the content in "Type":
		- If the content in "Type" is “LoF_IntraProtein_Enriched”, the format of this column is `{gene}:{number of LoF mutations in all samples}`
		- Otherwise, the format of this column is `{residue1}:{number of mutated samples},{residue2}:{number of mutated samples},...`
	6. LoF_enrichment (`[NA]` means not applicable)
	7. Raw_pvalue
	8. Adjusted_pvalue
	9. Subnetwork_ID 
	
		- If a significant 3D cluster or LOF enrichment signal can be mapped to a subnetwork, this field will contain the information about the subnetwork
		- If a significant 3D cluster or LOF enrichment signal can not be mapped to a subnetwork, this field will be `[NA]`
- `{job_name}`_drivers.txt

	This is a tab-separated file containing the candidate driver mutations identified by NetFlow3D. The first line is a header. The columns include:
	1. All columns in the input MAF file
	2. UniProt
	3. Signature_ID (This column contains the significant 3D cluster(s) or LOF enrichment signal by which a mutation is incorporated. If a mutation is incorporated by multiple significant 3D clusters, their IDs will be separated by comma)

- `{job_name}`_subnetworks.txt

	This is a tab-separated file containing the subnetworks with strong internal heat exchanges identified by NetFlow3D. Subnetworks including >=2 proteins will be reported. The first line is a header. Two columns are present:
	1. Subnetwork_genes	
	2. Subnetwork_size (i.e. number of proteins in the subnetwork)


