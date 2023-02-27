# NetFlow3D
NetFlow3D is a computational tool aiming at mapping how somatic mutations act across scales in cancer. If you find our tool helpful, please cite. You can also upload your data to our web server (http://netflow3d.yulab.org) and run NetFlow3D there.

## Prerequisites
The Python Standard Library and the following packages:
- scipy
- numpy
- networkx
- pandas
- statsmodels

## Installation

	git clone https://github.com/hyulab/NetFlow3D.git
	cd NetFlow3D
	python NetFlow3D.py -h
	
## Usage
Your command should be in the following format (the contents in `[]` are optional):

	python NetFlow3D.py -m <input_maf> -R <resolution> -I <job_name> [-X <expressed_genes>] [-n <binary_interactome>] [-o <output_path>] [-L <logfile_path>] [-t <threads>]

### Required arguments
- `-m <input_maf>`: replace `<input_maf>` with the path to your MAF file.
- `-R <resolution>`: replace `<resolution>` with `low` or `high`. This argument specifies the resolution for identifying subnetworks with strong internal heat exchanges. By default, `low` is recommended. If set to `high`, large subnetworks may be split into smaller subnetworks.
- `-I <job_name>`: replace `<job_name>` with a name you preferred for the current job.

### Optional arguments
- `-X <expressed_genes>`: replace `<expressed_genes>` with the path to your file which stores a complete list of expressed genes/proteins (see [Optional input](#optional-input) for how to generate the file). If not specified, all genes/proteins will be considered expressed.
- `-n <binary_interactome>`: replace `<binary_interactome>` with the path to your file which stores a complete list of existing protein-protein interactions (see [Optional input](#optional-input) for how to generate the file). If not specified, NetFlow3D will use the high quality binary interactome of Homo sapiens curated by HINT (http://hint.yulab.org/).
- `-o <output_path>`: replace `<output_path>` with a directory where the output files will be stored. If not specified, the output files will be stored in `./output/`.
- `-L <logfile_path>`: replace `<logfile_path>` with a directory where the log file will be stored. If not specified, the log file will be stored in `./log/`.
- `-t <threads>`: replace `<threads>` with a postive integer. This argument specifies the number of threads to use. If not specified, NetFlow3D will use 5 threads.
	
We provide example input files in `./example/input/`. Here is an example of your command (please run the following command to see if NetFlow3D is working normally):

	python NetFlow3D.py -m example/input/mutations.maf -R low -I test -X example/input/expressed_genes.txt -n example/input/interactome.txt -t 10
	
If you run the above command correctly, the output files should be found in `./output/` with prefix `test_`. To get an idea of what the output files should look like, please see example output files in `./example/output/`.


## Prepare input files
### Required input
- A Mutation Annotation Format (MAF) file (https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format)

	Required columns:
	<ul>
	<li>Hugo_Symbol</li>
	<li>Variant_Classification</li>
	<li>ENSP</li>
	<li>Protein_position</li>
	<li>Codons</li>
	<li>Tumor_Sample_Barcode</li>
	</ul>
	
	Optional columns:
	<ul>
	<li>Transcript_ID</li>
	<li>Gene</li>
	</ul>

	Other columns can also be present in the MAF file but they will not be used. 

### Optional input
- A text file containing a complete list of genes/proteins expressed in the cells where the mutations occur. One ID per line. Gene name, Ensembl gene ID, Ensembl transcript ID, Ensembl protein ID, and UniProt ID are accepted. Example:

	>ENSG00000163166<br>
	>ENSG00000110422<br>
	>ENSG00000077312<br>
	>ENSG00000180660<br>
	>ENSG00000186635<br>

- A text file containing a complete list of protein-protein interactions existing in the cells where the mutations occur. One interaction per line. Protein IDs should be separated by tab. Ensembl protein ID and UniProt ID are accepted. Example:

	>Q9H4A3&emsp;Q9HBL0<br>
	>Q15654&emsp;Q15797<br>
	>P63279&emsp;Q13643<br>
	>O43236&emsp;O43236<br>
	>P01112&emsp;P04049<br>


## Output files
NetFlow3D will output the following files. `{job_name}` will be replaced by the name you specified. If you run the example command, `{job_name}` will be replaced by `test`. 
- `{job_name}`_signatures.txt

	This a tab-separated file containing the selection signatures identified by NetFlow3D. The first line is a header. Eight columns are present:
	1. Signature_ID
	2. Type
	3. Affected_genes
	4. Structure_source (`[NA]` means not applicable)
	5. Mutation_frequency 

		The content format in this column depends on the content in "Type":
		- If the content in "Type" is “LoF_IntraProtein_Enriched”, the format of this column is `{gene}:{# of LoF mutations in all samples}`
		- Otherwise, the format of this column is `{residue1}:{# of mutated samples},{residue2}:{# of mutated samples},...`
	6. LoF_enrichment (`[NA]` means not applicable)
	7. Raw_pvalue
	8. Adjusted_pvalue
	9. Subnetwork_ID 
	
		- If a selection signature can be mapped to a subnetwork, this field will contain the information about the subnetwork
		- If a selection signature can not be mapped to a subnetwork, this field will be `[NA]`
- `{job_name}`_drivers.txt

	This is a tab-separated file containing the potentionally functional mutations identified by NetFlow3D. The first line is a header. The columns include:
	1. All columns in the input MAF file
	2. UniProt
	3. Signature_ID (This column indicates the selection signature where a mutation is involved. If a mutation is involved in multiple selection signatures, their IDs will be separated by comma)

- `{job_name}`_subnetworks.txt

	This is a tab-separated file containing the subnetworks with strong internal interconnectivity identified by NetFlow3D. Subnetworks including >=2 proteins will be reported. The first line is a header. Two columns are present:
	1. Subnetwork_genes	
	2. Subnetwork_size


