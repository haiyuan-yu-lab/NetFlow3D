# Net3D
Net3D 是用来decode somatic mutations干嘛的。 我们还有一个web server也可以用来run Net3D on your data. 如果你用它请cite paper link。

## Installation
Net3D is available on PyPI, which means you can install it with the following command:

	pip install Net3D

Alternatively, you can clone this repository to a local directory:

	git clone https://github.com/zzyingying753/Net3D.git
	cd Net3D

## Prerequisites
The Python Standard Library and the following packages:
- scipy
- numpy
- networkx (2.5 or later)
- pandas
- statsmodels

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

## Usage
Type your command in the following format. The content in `[]` are optional.

	python Net3D.py -m <input_maf> -R <resolution> -I <job_name> [-X <expressed_genes>] [-n <binary_interactome>] [-o <output_path>] [-L <logfile_path>] [-t <threads>]

### Required arguments
- `-m <input_maf>`: replace `<input_maf>` with the path to your MAF file.
- `-R <resolution>`: replace `<resolution>` with `low` or `high`. This argument specifies the resolution when identifying strongly connected subnetworks. If set to `high`, smaller subnetworks are favored. If set to `low`, larger subnetworks are favored.
- `-I <job_name>`: replace `<job_name>` with a name you preferred for the current job.

### Optional arguments
- `-X <expressed_genes>`: replace `<expressed_genes>` with the path to your file which stores a complete list of expressed genes/proteins (see [Optional input](#optional-input) for how to generate the file). If not specified, all genes/proteins will be considered to be expressed.
- `-n <binary_interactome>`: replace `<binary_interactome>` with the path to your file which stores a complete list of existing protein-protein interactions (see [Optional input](#optional-input) for how to generate the file). If not specified, Net3D will use the high quality binary interactome of Homo sapiens curated by HINT (http://hint.yulab.org/).
- `-o <output_path>`: replace `<output_path>` with a directory where the output files will be stored. If not specified, the output files will be stored in `./output/`.
- `-L <logfile_path>`: replace `<logfile_path>` with a directory where the log file will be stored. If not specified, the log file will be stored in `./log/`.
- `-t <threads>`: replace `<threads>` with a postive integer. This argument specifies the number of threads to use. If not specified, Net3D will use five threads.
	
An example of your command (please run the following command to see if Net3D is working normally):

	python Net3D.py -m example/input/mutations.maf -R low -I test -X example/input/expressed_genes.txt -n example/input/interactome.txt -t 10

## Output files
Net3D will output the following files. `{job_name}` will be replaced by the name you specified. If you run the example command, `{job_name}` will be replaced by `test`. 
- `{job_name}`_signatures.txt

	This a tab-separated file containing the selection signatures identified by Net3D. The first line is a header. Eight columns are present:
	1. Signature_ID
	2. Type
	4. Affected_genes
	5. Structure_source (`[NA]` means not applicable)
	6. Mutation_frequency 

		The content format in this column depends on the content in "Type":
		- If the content in "Type" is “LoF_IntraProtein_Enriched”, the format of this column is `{gene}:{# of LoF mutations in all samples}`
		- Otherwise, the format of this column is `{residue1}:{# of mutated samples},{residue2}:{# of mutated samples},...`
	9. LoF_enrichment (`[NA]` means not applicable)
	10. Raw_pvalue
	11. Adjusted_pvalue

- `{job_name}`_drivers.txt

	This is a tab-separated file containing the potentionally functional mutations identified by Net3D. The first line is a header. The columns include:
	1. All columns in the input MAF file
	2. UniProt
	3. Signature_ID (This column indicates the selection signature where a mutation is involved. If a mutation is involved in multiple selection signatures, their IDs will be separated by comma)

- `{job_name}`_subnetworks.txt

	This is a tab-separated file containing the subnetworks with strong internal interconnectivity identified by Net3D. The first line is a header. Two columns are present:
	1. Subnetwork_genes	
	2. Subnetwork_size


