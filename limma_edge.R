# R function for performing differential expression between genders for RNA-seq
# data.

require(edgeR)


limma_edge <- function(expression_df, design_df){
	# Create DGEList object
	design_df <- as.matrix(design_df)
	dge <- DGEList(counts=expression_df)
	
	# Calculate normalization factors to scale raw library sizes
	dge <- calcNormFactors(dge)
	
	# Filter counts for genes which are 'expressed enough', at least a cpm count of 1
	cpm_counts <- cpm(dge)
	keep <- rowSums(cpm_counts > 1) >= 1
	dge <- dge[keep,, keep.lib.sizes=FALSE]  # Update lib size too
	
    # Transform RNA-Seq for linear modelling and fit model - use quantile norm?
	voom_obj <- voom(dge, design=design_df, normalize="quantile")
	fit <- lmFit(voom_obj, design_df)

	# Make contrast matrix and fit prior to differential expression
	cont_matrix <- makeContrasts(expression_diff=female-male, levels=c('male', 'female'))
	fit <- contrasts.fit(fit, contrasts=cont_matrix)

	# Compute differential expression stats using empirical Bayes statistics
	fit <- eBayes(fit)

	# Get dataframe of top ranked genes 
	limma_dataframe <- topTable(fit, number=nrow(expression_df))
	limma_dataframe$gene_symbol <- rownames(limma_dataframe)

	return(limma_dataframe)
}
