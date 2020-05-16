import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r, Formula
from rpy2.robjects.packages import importr
r.source('limma_edge.R')

pandas2ri.activate()


def differential_expression_tcga_gender(cancer='laml', genenames=None):
    """A function to perform differential expression on TCGA RNA-seq data.
    
    Parameters
    ----------
    cancer : string {'brca', 'kirc', 'luad', 'ucec', 'thca', 'lusc', 'prad', 
                     'hnsc', 'lgg', 'coad', 'skcm', 'blca', 'lihc', 'stad', 
                     'ov', 'kirp', 'cesc', 'sarc', 'pcpg', 'paad', 'read', 
                     'gbm', 'esca', 'tgct', 'laml', 'thym', 'kich', 'meso', 
                     'uvm', 'acc', 'ucs', 'dlbc', 'chol'}
        Type of cancer to perform differential expression on.

    genenames : pandas Dataframe 
        DataFrame of gene names to join to limma output. Index column must be
        Ensembl ID, and the other column desired gene names.

    Returns
    -------
    limma_object : pandas DataFrame
        Limma object containing results of differential expression.
    """
    if genenames is None:
        gencode_names = pd.read_csv('gencode.v22.genes.txt', sep='\t')[['gene_id', 'gene_name', 'seqname']]
        gencode_names.set_index('gene_id', inplace=True)

    # Remove Y genes from gencode
    gencode_names = gencode_names[~(gencode_names['seqname'] == 'chrY')]

    expression_data = pd.read_csv("data/merged_rna/merged_rna_" + cancer + ".tsv", sep='\t').set_index('gene_name')
    
    # Remove Y genes from expression data
    expression_data = expression_data[expression_data.index.isin(gencode_names.index)]

    metadata = pd.read_csv("data/metadata/metadata.files_" + cancer + ".txt", sep='\t')
    metadata = metadata.set_index('cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id')

    metadata = metadata[~(metadata['cases.0.samples.0.sample_type'] == 'Solid Tissue Normal')]

    # Subset metadata by removing normal samples
    print(metadata['cases.0.samples.0.sample_type'].value_counts())

    male_sample_list = list(metadata[metadata['cases.0.demographic.gender'] == 'male'].index)
    female_sample_list = list(metadata[metadata['cases.0.demographic.gender'] == 'female'].index)

    design_matrix = pd.DataFrame([{'index': x, 
                                   'male': int(x in male_sample_list), 
                                   'female': int(x in female_sample_list)} 
                                   for x in expression_data.columns]).set_index('index')

    # Dictionary to hold expression data and design matrix
    processed_data = {'expression': expression_data, 'design': design_matrix}
    
    # Run through limma_edge function written in R
    limma_object = r.limma_edge(pandas2ri.py2rpy(processed_data['expression']), pandas2ri.py2rpy(processed_data['design']))

    # Load gencode names to join to limma_object
    # gencode_names = pd.read_csv('gencode.v22.genes.txt', sep='\t')[['gene_id', 'gene_name']]
    # gencode_names.set_index('gene_id', inplace=True)

    limma_object = limma_object.join(gencode_names)

    print(cancer, len(male_sample_list), len(female_sample_list),
            sum(limma_object['adj.P.Val'] < 0.05))

    print(limma_object[((limma_object['adj.P.Val'] <= 0.05) & (abs(limma_object['logFC']) >= 1.5))])

    return limma_object


if __name__=='__main__':
    # ['brca', 'kirc', 'luad', 'ucec', 'thca', 'lusc', 'prad', 'hnsc', 'lgg', 'coad', 'skcm', 'blca', 'lihc', 
    #            'stad', 'ov', 'kirp', 'cesc', 'sarc', 'pcpg', 'paad', 'read', 'gbm', 'esca', 'tgct', 'laml', 'thym', 
    #            'kich', 'meso', 'uvm', 'acc', 'ucs', 'dlbc', 'chol']

    deg_df = pd.DataFrame()

    for cancer in ['esca','kirc', 'luad', 'thca', 'lusc', 'hnsc', 'lgg', 'coad', 'skcm', 'blca', 'lihc', 
                   'stad', 'kirp', 'sarc', 'pcpg', 'paad', 'read', 'gbm', 'laml', 'thym', 
                   'kich', 'meso', 'uvm', 'acc', 'dlbc', 'chol']:

        deg = differential_expression_tcga_gender(cancer=cancer)
        deg = deg[((deg['adj.P.Val'] <= 0.05) & (abs(deg['logFC']) >= 1.5))]
        deg['cancer'] = cancer
        deg = deg[['gene_name', 'cancer']]
        deg_df = deg_df.append(deg)

        # gencode_names = pd.read_csv('gencode.v22.genes.txt', sep='\t')[['gene_id', 'gene_name']]
        # gencode_names.set_index('gene_id', inplace=True)

        # expression_data = pd.read_csv("data/merged_rna/merged_rna_" + cancer + ".tsv", sep='\t').set_index('gene_name')
        # metadata = pd.read_csv("data/metadata/metadata.files_" + cancer + ".txt", sep='\t')
        # metadata = metadata.set_index('cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id')

        # print(cancer)
        # print(metadata['cases.0.samples.0.sample_type'].value_counts())
    
    count_df = deg_df.groupby('gene_name').count()
    count_df.to_csv('gene_count_df.csv')
    deg_df.to_csv('def_df.csv')
