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
        gencode_names = pd.read_csv('gencode.v22.genes.txt', sep='\t')[['gene_id', 'gene_name']]
        gencode_names.set_index('gene_id', inplace=True)

    expression_data = pd.read_csv("data/merged_rna/merged_rna_" + cancer + ".tsv", sep='\t').set_index('gene_name')
    metadata = pd.read_csv("data/metadata/metadata.files_" + cancer + ".txt", sep='\t')
    metadata = metadata.set_index('cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id')

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
    gencode_names = pd.read_csv('gencode.v22.genes.txt', sep='\t')[['gene_id', 'gene_name']]
    gencode_names.set_index('gene_id', inplace=True)

    limma_object = limma_object.join(gencode_names)

    print(cancer, len(male_sample_list), len(female_sample_list),
            sum(limma_object['adj.P.Val'] < 0.01))

    return limma_object


if __name__=='__main__':
    # ['brca', 'kirc', 'luad', 'ucec', 'thca', 'lusc', 'prad', 'hnsc', 'lgg', 'coad', 'skcm', 'blca', 'lihc', 
    #            'stad', 'ov', 'kirp', 'cesc', 'sarc', 'pcpg', 'paad', 'read', 'gbm', 'esca', 'tgct', 'laml', 'thym', 
    #            'kich', 'meso', 'uvm', 'acc', 'ucs', 'dlbc', 'chol']
    for cancer in ['chol', 'kirc', 'luad', 'thca', 'lusc', 'hnsc', 'lgg', 'coad', 'skcm', 'blca', 'lihc', 
               'stad', 'kirp', 'sarc', 'pcpg', 'paad', 'read', 'gbm', 'esca', 'laml', 'thym', 
               'kich', 'meso', 'uvm', 'acc', 'dlbc', 'chol']:

        print(differential_expression_tcga_gender(cancer=cancer).head(20))
        from bokeh.plotting import figure, output_file, show, ColumnDataSource
        from bokeh.io import output_notebook
        from bokeh.models import HoverTool, LogColorMapper
        import numpy as np

        p = figure(plot_width=800, plot_height=400, title='TCGA-LAML Men vs Women Gene Signature Volcano Plot')
        Limma_Data = differential_expression_tcga_gender(cancer=cancer)
        # show the results

        source = ColumnDataSource(data=dict(
            x=Limma_Data['logFC'],
            y=-1*np.log10(Limma_Data['P.Value']),
            gene=Limma_Data['gene_name'],
            gene_id=Limma_Data['gene_symbol']
        ))

        TOOLTIPS = [
            ("Gene Name", "@gene"),
            ("Ensembl ID", "@gene_id"),
            ("-log10P", "$y")
        ]

        h = HoverTool(tooltips = TOOLTIPS)

        p = figure(plot_width=900, plot_height=550, title="TCGA-KIRC Normal vs Tumor Gene Expression Signature Volcano Plot")

        p.circle('x', 'y', size=4, source=source)

        p.xaxis.axis_label = "log(Fold-Change)"
        p.yaxis.axis_label = "-log10(P-Value)"

        p.add_tools(h)

        show(p)