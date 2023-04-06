#!/usr/bin/env python

import re
import sys
import json
import argparse
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from joblib import Parallel, delayed

description_text = (
    """
    Synopsis:
    
    phenven.py --terms patient_phenotype.tsv         \
               --gene_file candidate_genes.tsv       \
               --phen2gene phenotype_to_genes.txt.gz \
               --json hp.json.gz                     \
               --jobs 4                              \
               > phenotype_overlap.tsv
    
    Description:

    PhenVen is a tool to assist with a review of the phenotype overlap
    between a patient and a set of candidate genes.  PhenVen can
    prioritize candidate genes and/or diseases based on a patient's
    phenotype in a diagnosis style analysis, similar to how you might
    use Phevor or Phenomizer (phenotype first analysis).  PhenVen
    might also be run with a list of candidate genes that were
    identified by WGS/WES gene/variant prioritization analysis to aid
    phenotype review of those candidate genes.  Finally, PhenVen could
    be used to provide detail and clarity regarding the extent of
    phenotype overlap between the phenotype terms of a patient and one
    or more candidate genes that have already been
    phenotype-prioritized by a tool like Phevor or Exomiser - e.g. you
    know/suspect that there is good phenotype overlap, but you want to
    review which phenotypes are driving that signal.

    """
)

def main():
    parser = argparse.ArgumentParser(
        description=description_text,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--genes', '-g',
                        help='A comma-separated list of candidate genes.')
    parser.add_argument('--gene_file', '-f',
                        help='A file of candidate genes - one per row first column.')
    parser.add_argument('--terms', '-t', dest='proband_terms_file',
                        help='A file containing the list of HPO IDs (first column) associated with the proband.')
    parser.add_argument('--phen2gene', '-p', dest='phen2gene_file',
                        help='The gzipped phenotype_to_genes.txt.gz file from HPO')
    parser.add_argument('--json', '-j', dest='json_file',
                        help='A gzipped hp.json.gz file from HPO')
    parser.add_argument('--jobs', '-n', type=int, default=1,
                        help='The number of jobs to run in parallel')
    args = parser.parse_args()

    # Parse Text Files to Datafames
    sys.stderr.write('INFO : loading_data_file : ' + args.gene_file + '\n')
    df_cnd = object()
    if args.genes is not None:
        gene_list = args.genes.split(',')
        df_cnd = pd.DataFrame(gene_list, columns=['gene'])
    if args.gene_file is not None:
        df_cnd = pd.read_table(args.gene_file, names=['gene'])

    sys.stderr.write('INFO : loading_data_file : ' + args.proband_terms_file + '\n')
    df_prb = pd.read_table(args.proband_terms_file)
    df_prb.drop_duplicates(subset='id', inplace=True)
    df_prb.set_index('id', inplace=True)
    prb_ids = set(df_prb.index.tolist())
    
    # Parse HPO Phenotype_to_Gene File to Dataframe
    # 'id' 'term' 'entrez-gene-id' 'entrez-gene-symbol' 'Additional Info from G-D source' 'G-D source' 'disease-ID for link'
    sys.stderr.write('INFO : loading_data_file : ' + args.phen2gene_file + '\n')
    df_p2g = pd.read_table(args.phen2gene_file,
                           skiprows=1,
                           header=None,
                           index_col=False,
                           names=['id', 'term', 'gene_id', 'gene', 'source_info', 'source_id', 'disease_id'])
    df_p2g.drop_duplicates(subset=['id', 'gene'], inplace=True)
    df_p2g = df_p2g.set_index('id')
        
    # Get term frequency and information content (ic)
    term_counts = pd.DataFrame(df_p2g.reset_index()['id'].value_counts())
    term_counts.rename(columns={'id': 'count'}, inplace=True)
    term_counts.index.rename('id', inplace=True)
    term_counts['freq'] = term_counts['count'] / term_counts['count'].sum()
    term_counts['ic'] = -1 * np.log10(term_counts['freq'])
    term_counts.loc['HP:0000118'] = [term_counts['count'].sum(), 1, 0]
    df_p2g = df_p2g.join(term_counts)
    
    # Parse HPO JSON
    sys.stderr.write('INFO : loading_data_file : ' + args.json_file + '\n')
    with open(args.json_file) as json_fh:
        data = json.load(json_fh)
        
    # Parse Nodes from HPO JSON Data
    sys.stderr.write('INFO : processing_hpo_graph\n')
    nodes = list()
    for node in data['graphs'][0]['nodes']:
        hp_id = node['id']
        hp_id = re.sub('.*\/', '', hp_id)
        hp_id = re.sub('_', ':', hp_id)
        if not re.match('^HP:', hp_id):
            continue
        hp_lbl = node['lbl']
        hp_def = 'No definition'
        if 'meta' in node.keys():
            if 'definition' in node['meta'].keys():
                hp_def = node['meta']['definition']['val']
                nodes.append((hp_id, {'term': hp_lbl, 'definition': hp_def}))
                
    # Parse Edges from HPO JSON Data
    edges = []
    for edge in data['graphs'][0]['edges']:
        sub = edge['sub']
        sub = re.sub('.*\/', '', sub)
        sub = re.sub('_', ':', sub)
        
        obj = edge['obj']
        obj = re.sub('.*\/', '', obj)
        obj = re.sub('_', ':', obj)
        
        edges.append((obj, sub))

    # Create NetworkX Graph of HPO
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    # Trim HPO to subgraph of 'Phenotypic abnormality'
    root = 'HP:0000118' # Phenotypic abnormality
    sub_ids = nx.descendants(G, 'HP:0000118')
    sub_ids.add('HP:0000118')
    pabG = G.subgraph(sub_ids)

    # Get subgraph/leaves of proband HPO ancestors
    prb_id_ancestors = set()
    for id in prb_ids:
        prb_id_ancestors.add(id)
        if id not in pabG.nodes:
            # WARN
            continue
        prb_id_ancestors.update(nx.ancestors(pabG, id))
        prb_id_ancestors.add('HP:0000118')
        prbG = pabG.subgraph(prb_id_ancestors)
        prb_leaves = {node for node in prbG.nodes() if prbG.out_degree(node)==0}



    sys.stderr.write('INFO : identify_lcas_and_shpl \n')
    all_df_lcas = []
    all_gene_lcas = []
    # Check for gene_ids in genG_ids
    all_df_lcas = Parallel(n_jobs=args.jobs)(delayed(get_all_lcas)(prb_id_ancestors, prb_leaves, gene, df_p2g, pabG)
                                             for gene in tqdm(df_cnd['gene'].tolist()))

    sys.stderr.write('INFO : processing_term_pairs \n')
    df_lcas = pd.concat(all_df_lcas)
    df_lcas.set_index('lca_id', inplace=True)
    df_lcas = df_lcas.join(term_counts, how='left')
    df_lcas.sort_values(by='ic', ascending=False, inplace=True)
    df_lcas.drop_duplicates(subset=['prb_id', 'gene'], keep='first', inplace=True)

    sys.stderr.write('INFO : merging_gene_data \n')
    df_lcas = df_lcas.merge(df_p2g['term'], left_on='prb_id', right_index=True).drop_duplicates()
    df_lcas.rename(columns={'term': 'prb_term'}, inplace=True)

    df_lcas = df_lcas.merge(df_p2g['term'], left_on='gene_id', right_index=True).drop_duplicates()
    df_lcas.rename(columns={'term': 'gene_term'}, inplace=True)

    df_lcas = df_lcas.merge(df_p2g['term'], left_on='anc_id', right_index=True).drop_duplicates()
    df_lcas.rename(columns={'term': 'anc_term'}, inplace=True)

    df_lcas = df_lcas.merge(df_p2g['term'], left_index=True, right_index=True).drop_duplicates()
    df_lcas.rename(columns={'term': 'lca_term'}, inplace=True)

    sys.stderr.write('INFO : sorting_data \n')
    df_lcas.sort_values(by='ic', ascending=False, inplace=True)
    df_lcas.reset_index(inplace=True, drop=False)
    df_lcas.rename(columns={'index': 'lca_id'}, inplace=True)

    sys.stderr.write('INFO : print_data \n')
    # Fix ordering here
    df_lcas = df_lcas.loc[:,['gene', 'ic', 'shpl', 'prb_term', 'gene_term',
                             'lca_term', 'anc_term', 'prb_id', 'gene_id',
                             'lca_id', 'anc_id', 'count', 'freq']]
    print(df_lcas.to_csv(sep='\t', index=False))

def get_lcas(prb_id, genes, pgG, pgUG):

    lcas = []
    for gene_id in genes:
        lca = nx.lowest_common_ancestor(pgG, gene_id, prb_id)
        paths = nx.all_simple_paths(pgG, 'HP:0000118', prb_id)
        anc = lca
        for path in paths:
            if (1 < len(path)):
                anc = path[1]
                break
        shpl = nx.shortest_path_length(pgUG, gene_id, prb_id)
        lcas.append((anc, prb_id, gene_id, lca, shpl))
    return lcas

def get_all_lcas(prb_id_ancestors, prb_id_leaves, gene, df_p2g, pabG):
    # Get subgraph/leaves of gene HPO ancestors
    # gene = df_cnd.iloc[0]['gene']
    gene_ids = set(df_p2g.query('gene == @gene').index.tolist())
    gene_id_ancestors = set()
    for id in gene_ids:
        gene_id_ancestors.add(id)
        if id not in pabG.nodes:
            # WARN
            continue
        gene_id_ancestors.update(nx.ancestors(pabG, id))
    gene_id_ancestors.add('HP:0000118')
    genG = pabG.subgraph(gene_id_ancestors)
    gen_leaves = [node for node in genG.nodes() if genG.out_degree(node)==0]

    # Get subgraph of proband & gene IDs
    prb_gene_ids = prb_id_ancestors.union(gene_id_ancestors)
    prb_gene_ids.add('HP:0000118')
    pgG = pabG.subgraph(prb_gene_ids)
    pgG_ids = set(pgG.nodes)
    pgUG = nx.Graph(pgG)
        
    lcas = []
    for prb_id in prb_id_leaves:
        for gene_id in gen_leaves:
            lca = nx.lowest_common_ancestor(pgG, gene_id, prb_id)
            paths = nx.all_simple_paths(pgG, 'HP:0000118', prb_id)
            anc = lca
            for path in paths:
                if (1 < len(path)):
                    anc = path[1]
                    break
            shpl = nx.shortest_path_length(pgUG, gene_id, prb_id)
            lcas.append((anc, prb_id, gene_id, lca, shpl))
    
    # all_gene_lcas = [];
    # for lca in lcas:
    #     all_gene_lcas.extend(lca)
    
    df_lcas = pd.DataFrame(lcas, columns=['anc_id', 'prb_id', 'gene_id', 'lca_id', 'shpl'])
    df_lcas['gene'] = gene
    return df_lcas



if __name__ == "__main__":
    main()
