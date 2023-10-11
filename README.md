# PhenVen - Patient to HPO gene/disease phenotype overlap analysis.

### Description

PhenVen is a tool to assist with a review of the phenotype overlap
between a patient and genes/diseases anotated by the [Human Phenotype Ontology](https://hpo.jax.org/app/).  PhenVen was primarily 
designed as a tool to provide the detail needed to explore the details of phenotype 
overlap between the phenotype terms of a patient and one or more 
candidate genes e.g. you know/suspect that there is good phenotype overlap, but you
need to review a detailed list of phenotypes are driving that signal.

PhenVen can also 
prioritize candidate genes and/or diseases based on a patient's
phenotype in a diagnostic style analysis, similar to how you might use
[Phenomizer](https://pubmed.ncbi.nlm.nih.gov/19800049/) in a phenotype first
analysis.  PhenVen might also be run with a list of candidate genes 
that were identified by WGS/WES gene/variant prioritization methods 
to add a phenotype review component to that analysis - similar to how you 
might use [Phevor](https://pubmed.ncbi.nlm.nih.gov/24702956/), but with ability to explore the details of the phenotype overlap.

### Synopsis

```
phenven.py -terms hpo_terms.tsv --genes genes.tsv \
           -phen2gene phenotype_to_genes.txt -json hp.json -n 16
```

### Documentation

```
 Arguments:
  -h, --help            show this help message and exit
  --genes GENES, -g GENES
                        A comma-separated list of target/candidate genes. (default: None)
  --gene_file GENE_FILE, -f GENE_FILE
                        A file of target/candidate genes - one per row first column. (default: None)
  --terms PROBAND_TERMS_FILE, -t PROBAND_TERMS_FILE
                        A file containing the list of HPO IDs (first column) associated with the proband. (default: None)
  --phen2gene PHEN2GENE_FILE, -p PHEN2GENE_FILE
                        The gzipped phenotype_to_genes.txt.gz file from HPO (default: None)
  --json JSON_FILE, -j JSON_FILE
                        A gzipped hp.json.gz file from HPO (default: None)
  --jobs JOBS, -n JOBS  The number of jobs to run in parallel (default: 1)
```

### Installation

Downloading phenven.py Python script or cloning the PhenVen repo from GitHub is all you need as code from PhenVen.  The phenven.py script imports several Python libraries that are required dependencies and will need to be installed (they are not core Python libraries) if you don't already have them.  It is recommended that you use conda or pip to install the following libraries:

* numpy
* pandas
* tqdm
* joblib
* networkx

```
conda create env phenven
conda activate phenven
conda install numpy pandas tqdm joblib networkx
```

### Getting Help

Sorry, no fancy ReadtheDocs pages or Discord server for this repo.  You can post an issue to this GitHub repo if you need to report a bug.

### License

PhenVen is license under the GNU General Public License Version 3

### Author

PhenVen was written by Barry Moore barry.moore@genetics.utah.edu
