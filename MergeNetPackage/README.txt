MergeNetPackageVer 1.0, released on 06/01/2012.
This archive contains the excutable programs which have been tested on several Linux machines using the standard g++ compiler.

Email netlab@osumc.edu with questions or comments.

This README describes how to use MergeNetPackage to mine biomedical network data. 
Reference: Yang Xiang, David Fuhry, Kamer Kaya, Ruoming Jin, Umit Catalyurek, Kun Huang, "Merging Network Patterns: A General Framework to Summarize Biomedical Network Data", Network Modeling Analysis in Health Informatics and Bioinformatics, In press, DOI:10.1007/s13721-012-0009-3
In addition, this package uses available software "mafia" (http://himalaya-tools.sourceforge.net/Mafia/). 

You may use our program for research, publication and education purposes, provided that you properly cite our paper in your publications.
In the following we illustrate how to use this package to summarize unweighted bipartite graph patterns, and weighted graph patterns:

Part I: Summarize unweighted bipartite graph patterns
Dataset: genes_to_phenotype_hpo-Revision2250, downloaded from http://human-phenotype-ontology.org/ on Feb 23, 2012
The dataset is preprocessed into three parts:
(1) genes_to_phenotype_hpo-Revision2250-Feb23-2012.dat, in folder "MAFIA-MBC/datasets". This data is in transactional data format, where each line corresponds to a gene and each number in the line corresponds to a phenotype which the gene is related to.
(2) genes_to_phenotype_hpo-Revision2250-Feb23-2012_GeneList.txt, in folder "SumNetwork/datasets". Each line is a gene name of the corresponding line in the file "genes_to_phenotype_hpo-Revision2250-Feb23-2012.dat".  
(3) genes_to_phenotype_hpo-Revision2250-Feb23-2012_PhenotypeList.txt, in folder "SumNetwork/datasets". Each line is a phenotype name. The line number, starting from 1, is the corresponding phenotype ID in the file "genes_to_phenotype_hpo-Revision2250-Feb23-2012.dat".

Step 1: Get Closed Frequent Itemset:
In directory "MAFIA-MBC/", run:
qsub script_file_mafia_g2f_005

Step 2: Get maximal bicliques based on the closed frequent itemset:
In directory "MAFIA-MBC/", run:
qsub script_file_MBC_g2f_005

Step 3: (First round) Merge maximal bicliques using SingleMerge strategy with density threshold being 0.7:
In directory "SumNetwork/", run:
qsub script_file_GraphMining_g2p_S_07

Step 4: (Second round) Merge maximal bicliques using MultiMerge strategy with density threshold being 0.7:
In directory "SumNetwork/", run:
qsub script_file_GraphMining_g2p_S_07_2nd

Step 5: Post processing the merging results into gene-phenotype list:
In directory "SumNetwork/", using Matlab run:
postprocessing_biclique.m 

Part II: Summarize weighted graph patterns
NCBI GEO datset GSE2034
The dataset is preprocessed into parts:
(1) GeneCoexpWang.matrix, not included in this package due to its size. The is a gene coexpression data built on GSE2034.
(2) GeneCoexpWang_MAFIA_06.dat, in folder "MAFIA-MBC/datasets". This dataset is a result of converting "GeneCoexpWang.matrix" to a transacation data for MAFIA to process. The edge cutoff threshold is 0.6. 
(3) GeneCoexpWang_genelist.txt, in folder "SumNetwork/datasets". Each line is a gene name of the corresponding row (or column) in the file "GeneCoexpWang.matrix". The line number, starting from 1, is the corresponding gene ID in the file "GeneCoexpWang_MAFIA_06.dat".

Step 1: Get Closed Frequent Itemset:
In directory "MAFIA-MBC/", run:
qsub script_file_mafia_Wang_06_01

Step 2: Get cliques based on the closed frequent itemset:
In directory "MAFIA-MBC/", run:
qsub script_file_MBC_Wang_06_clique_01

Step 3: (First round) Merge maximal bicliques using MultiMerge strategy with density threshold being 0.7:
In directory "SumNetwork/", run:
qsub script_file_GraphMining_Wang_0601_07

Step 4: (Second round) Merge maximal bicliques using SingleMerge strategy with density threshold being 0.7:
In directory "SumNetwork/", run:
qsub script_file_GraphMining_Wang_0601_07_2nd_S

Step 5: Post processing the merging results into gene list:
In directory "SumNetwork/", using Matlab run:
postprocessing_clique.m

Below are optional steps for additional mining/analysis as described in our paper:

Step 6: Performance survival tests on the discovered gene list
This requires a large amount of external data therefore it is not included.
To continue the following analysis, we put the result file "Wang_0601_07_2nd_S.survival" in directory "SumNetwork/Survival". Note that the float number at the end of each line is the p-value of survival test on the GSE1456 dataset.

Step 7: Preprocess (encode) the survival test result into network patterns:
In directory "SumNetwork/", using Matlab run:
preprocessing_clique.m

Step 8: Merging network patterns into macro level networks and build hierachical tree with density threshold being 0.7:
In directory "SumNetwork/", run:
script_file_GraphMining_Wang_0601_07_tree_d04

Step 9: Visualize macro level networks by prepare a visualization file for Gephi:
In directory "SumNetwork/", run:
script_file_GraphVis_Wang_0601_07_2nd_S_w02_d04
or
script_file_GraphVis_Wang_0601_07_2nd_S_w02_d04

Part III: Additional Hints:
(1) The best way to understand the above procedue is to read our paper.
(2) If your machine does not have a queuing system, please open the script file and use the command in the file.
(3) Simply running the executable file, MBC, GraphMining, GraphVis, and you will see the command line formats.
(4) Our program also supports summarizing weighted bipartite network patterns, and unweighted graph patterns, if properly setting up the workflow.


