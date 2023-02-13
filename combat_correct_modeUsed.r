library(sva)

combat_correct <- function(exp_df, studies_df, groups_df, cutoff, outFile){
	exp_matrix = as.matrix(exp_df)
	idx <- rowMeans(exp_matrix) > cutoff
	exp_matrix <- exp_matrix[idx,]
	mod = model.matrix(~as.factor(group), data=groups_df)
	exp_df_corr = ComBat(exp_matrix, batch=studies_df$studyID, mod=mod)
	write.table(exp_df_corr, file=outFile)
	return(exp_df_corr)
}