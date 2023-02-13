library(sva)

combat_correct <- function(exp_df, studies_df, cutoff, outFile){
	exp_matrix = as.matrix(exp_df)
	idx <- rowMeans(exp_matrix) > cutoff
	exp_matrix <- exp_matrix[idx,]
	exp_df_corr = ComBat(exp_matrix, batch=studies_df$studyID)
	write.table(exp_df_corr, file=outFile)
	return(exp_df_corr)
}