import umap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from sklearn.preprocessing import StandardScaler

gencode_names = pd.read_csv('gencode.v22.genes.txt', sep='\t')
gencode_names.set_index('gene_id', inplace=True)
y_gene_list = gencode_names[gencode_names['seqname'] == 'chrY']

for cancer in ['brca', 'kirc', 'luad', 'ucec', 'thca', 'lusc', 'prad', 'hnsc', 'lgg', 'coad', 'skcm', 'blca', 'lihc', 
               'stad', 'ov', 'kirp', 'cesc', 'sarc', 'pcpg', 'paad', 'read', 'gbm', 'esca', 'tgct', 'laml', 'thym', 
               'kich', 'meso', 'uvm', 'acc', 'ucs', 'dlbc', 'chol']:

    expression_data = pd.read_csv("data/merged_rna/merged_rna_" + cancer + ".tsv", sep='\t').set_index('gene_name')
    metadata = pd.read_csv("data/metadata/metadata.files_" + cancer + ".txt", sep='\t')
    metadata = metadata.set_index('cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id')

    male_sample_list = list(metadata[metadata['cases.0.demographic.gender'] == 'male'].index)
    design_matrix = pd.DataFrame([{'index': x, 
                                   'male': int(x in male_sample_list)} 
                                   for x in expression_data.columns]).set_index('index')
    
    
    full_umap_data_lab = expression_data.T.join(design_matrix)
    # Keep only Y genes
    full_umap_data = full_umap_data_lab.T[full_umap_data_lab.columns.isin(y_gene_list.index)].T
    
    test_df = full_umap_data.join(design_matrix)
    test_df = test_df[test_df['male'] == 0]
    test_df = test_df.loc[(test_df > 100).any(axis=1)]
    print(cancer, test_df)

    # Drop zero columns and normalize
    full_umap_data = full_umap_data.loc[:, (full_umap_data != 0).any(axis=0)]
    full_umap_data = (full_umap_data - full_umap_data.mean()) / full_umap_data.std()

    reducer = umap.UMAP()
    exp_umap = reducer.fit_transform(full_umap_data)

    # Plot
    plt.scatter(exp_umap[:, 0], exp_umap[:, 1], c=full_umap_data_lab['male'], 
                cmap=plt.cm.get_cmap('coolwarm', 2)) 
    
    plt.title('umap projection: ' +  cancer)
    plt.legend()
    plt.colorbar()
    plt.show()