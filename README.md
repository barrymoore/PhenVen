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
phenven.py -terms hpo_terms_all.tsv --genes candidate_genes.txt \
           -phen2gene phenotype_to_genes.txt -json hp.json -n 4
```

### Documentation


### Installation


### Docker


### Getting Help


### License

PhenVen is license under the GNU General Public License Version 3

### Author

PhenVen was written by Barry Moore barry.moore@genetics.utah.edu
