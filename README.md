# Comparative human/mouse pan-cancer analyses uncover regulations of cancer and cancer immune divergence by clade-specific lncRNAs

This directory contains data and scripts related to our work: "Comparative human/mouse pan-cancer analyses uncover regulations of cancer and cancer immune divergence by clade-specific lncRNAs".

In this work, we integrated 9058 samples from 13 human cancer types, their mouse counterparts, and the corresponding normal tissues to study human/mouse transcriptomic divergences in the cancer context. We also studied the relationships between divergent genes and clade-specific lncRNAs across cancers, and identified human-specific and cancer-specific modules of interactions between divergent genes and clade-specific lncRNAs.


# Data
The datasets were too big to be uploaded via "Upload Files", so they were hosted by "Release" in this repository, which can be accessed from https://github.com/LinjieCodes/ComparativePanCancer/releases/tag/Data .

1. **TPM**  --  the gene expression matrixes normalized by TPM.

2. **NX**  --  the gene expression matrixes normalized by z-score which were ultimately used for most subsequent analyses.

3. **LongTarget_output**  --  the raw output files of LongTarget running which predicts lncRNA/DNA bindings.

4. **TTS**  --  the predicted DNA binding sites of clade-specific lncRNAs, and their related cis-regulatory elements annotated by ENCODE project

# Scripts
1. **download_ENA_RNASeq-er_data.py**  --  to extract raw count data from the RNASeq-er API which has re-assembled RNA-seq reads deposited to the ENA website.

2. **extract_compara_orthologs.py**  --  to obtain human-mouse one-to-one orthologues from the Compara API on the EMBL-EBI website.

3. **tmm.r**  --  the TMM (trimmed means of M values) algorithm implemented in the edgeR package.

4. **tmm_tpm.py**  --  to normalize raw counts using the TMM and TPM methods to remove the bias of library size on expression quantification.

5. **combat_correct.r and combat_correct_modeUsed.r**  --  the batch effect-removing algorithm implemented in the ComBat package.

6. **combat_for_mouse.py and combat_for_human_modeUsed.py**  --  to adjust gene expression matrixes using ComBat to remove batch effects between different datasets.

7. **scale_geneExp.py**  --  to normalize the gene expression matrixes by the z-score method using the preprocessing.StandardScaler function in the scikit-learn package.

8. **anosim.py**  --  the ANOSIM analysis to test whether there are confounding factors would bias the cross-species comparisons.

9. **tsne.py**  --  the t-distributed Stochastic Neighbor Embedding (t-SNE) method to visualize high-dimensional NX profiles.

10. **diffExp_cancerNormal.py**  --  the differential expression analysis between cancer and normal samples.

11. **diffExp_interspecies.py**  --  the differential expression analysis between human and mouse samples.

12. **barplot_divergentGene_percentage.py**  --  to display the percentages of divergent genes across cancers.

13. **cancerExpressedLnc.py**  --  to identify clade-specific lncRNAs that are robustly expressed in each cancer type.

14. **cladeLnc_divergentGene_corr.py**  --  to calculate the expression correlation between clade-specific lncRNAs and divergent genes in each cancer.

15. **divergentGeneCorrelation.py**  --  to calculate the expression correlation between divergent genes in each cancer.

16. **coExpCluster.py**  --  to cluster clade-specific lncRNAs and divergent genes in each cancer into gene modules based on the expression correlation distances.

17. **coExpModule_immuneInfiltration.py**  --  to examine the relationships between gene modules and immune cell infiltration.

18. **compare_immune_infiltration2.py**  --  the human/mouse comparison of immune cell infiltration. 

19. **immuneDivergentModule_spearman_infiltration_eachGene.py**  --  to calculate the spearman correlation between genes in each gene module and immune cell infiltration. 

20. **markerEnrichment.py**  --  the marker gene enrichment analysis to reveal the relationships between gene modules and immune cells.

21. **percentage_module_enrichedWithLnc.py**  --  to calculate the percentages of gene modules that are enriched with clade-specific lncRNAs.

22. **plot_geneExample_immuneCellInfiltration_correlation.py**  --  to display examples in gene modules, and their correlation with immune cell infiltration.

23. **prioritize_immuneLncs.py**  --  to prioritize top lncRNA regulators from all immune divergence-associated modules across cancers.

24. **tts_cre_enrichment.py**  --  the CRE enrichment analysis to study the relationships between predicted lncRNA/DNA binding sites and ENCODE-annotated cis-regulatory elements.

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com.

# Related website
To obtain more details about lncRNA/DNA binding prediction using LongTarget, please go to our website http://lncRNA.smu.edu.cn .

# Citation
Lin et al. Comparative human/mouse pan-cancer analyses uncover regulations of cancer and cancer immune divergence by clade-specific lncRNAs.
