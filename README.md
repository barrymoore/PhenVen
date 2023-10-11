# PhenVen - Proband/gene/disease phenotype overlap analysis.

### Description

PhenVen is a tool to assist with a review of the phenotype overlap
between a patient and a set of candidate genes.  PhenVen can
prioritize candidate genes and/or diseases based on a patient's
phenotype in a diagnostic style analysis, similar to how you might use
Phevor or Phenomizer (phenotype first analysis).  PhenVen might also
be run with a list of candidate genes that were identified by WGS/WES
gene/variant prioritization analysis to add phenotype review of those
candidate genes.  However, PhenVen was primarily designed as a tool
to provide the detail needed to explore the phenotype overlap between the
phenotype terms of a patient and one or more candidate genes that have
already been phenotype-prioritized by a tool like Phevor or Exomiser -
e.g. you know/suspect that there is good phenotype overlap, but you
need to review a detailed list of phenotypes are driving that signal.

### Synopsis

```
phenven.py -terms hpo_terms_all.tsv --genes candidate_genes.txt \
           -phen2gene phenotype_to_genes.txt -json hp.json -n 4
```


patient_phenotype.tsv

../data/:
README.md
genes_to_phenotype.txt
hp.json
phenotype.hpoa
phenotype_to_genes.txt

### Documentation


### Installation


### Docker


### Getting Help


### License

PhenVen is license under the GNU General Public License Version 3

### Author

PhenVen was written by Barry Moore barry.moore@genetics.utah.edu
