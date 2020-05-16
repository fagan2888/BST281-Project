import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()


def expression_violinplot(gene='ALB', cancers=['gbm', 'lihc', 'acc', 'laml', 'thym']):
    """Take gene and make expression plot from subset of cancers."""

    gencode_names = pd.read_csv('gencode.v22.genes.txt', sep='\t')[['gene_id', 'gene_name']]
    gencode_names.set_index('gene_id', inplace=True)

    gene_ensembl = gencode_names[gencode_names['gene_name'] == gene]

    violin_df = pd.DataFrame()
    for cancer in cancers:
        expression_data = pd.read_csv("data/merged_rna/merged_rna_" + cancer + ".tsv", sep='\t').set_index('gene_name')
        expression_data = expression_data.loc[gene_ensembl.index]
        
        metadata = pd.read_csv("data/metadata/metadata.files_" + cancer + ".txt", sep='\t')
        metadata = metadata.set_index('cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id')

        male_sample_list = list(metadata[metadata['cases.0.demographic.gender'] == 'male'].index)
        female_sample_list = list(metadata[metadata['cases.0.demographic.gender'] == 'female'].index)

        design_matrix = pd.DataFrame([{'index': x, 
                                    'male': int(x in male_sample_list), 
                                    'female': int(x in female_sample_list)} 
                                    for x in expression_data.columns]).set_index('index')

        merged_df = pd.merge(expression_data.T, design_matrix,
                             left_index=True, right_index=True)

        merged_df['cancer'] = cancer

        violin_df = violin_df.append(merged_df)
    
    violin_df['gender'] = ['male' if x==1 else 'female' for x in violin_df['male'].values]

    print(violin_df)


    gene_ens = gene_ensembl.index.values[0]

    violin_df[gene_ens] = violin_df[gene_ens] / violin_df[gene_ens].std()
    sns.boxplot(x='cancer', y=gene_ens, hue='gender', data=violin_df)
    plt.title('TCGA Expression Data for ' + gene)
    plt.show()


if __name__=='__main__':
    expression_violinplot(gene='PI3', cancers=sorted(['kirc', 'luad', 'thca', 'lusc', 'hnsc', 'lgg', 'coad', 'skcm', 'blca', 'lihc', 
               'stad', 'kirp', 'sarc', 'pcpg', 'paad', 'read', 'gbm', 'esca', 'laml', 'thym', 
               'kich', 'meso', 'uvm', 'acc', 'dlbc', 'chol']))