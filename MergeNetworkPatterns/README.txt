SumNetwork Version 1.2. Released July 04, 2012. Copyright(C) The Network Laboratory of Biomedical Informatics, The Ohio State University (netlab.bmi.osumc.edu).

This archive contains the excutable programs which have been tested on several Linux machines using the standard g++ compiler.

Email netlab@osumc.edu with questions or comments.

This README describes how to use MergeNetPackage to mine biomedical network data. 
Reference: Yang Xiang, David Fuhry, Kamer Kaya, Ruoming Jin, Umit Catalyurek, Kun Huang, "Merging Network Patterns: A General Framework to Summarize Biomedical Network Data", Network Modeling Analysis in Health Informatics and Bioinformatics, DOI:10.1007/s13721-012-0009-3
In addition, this package uses available software "mafia" (http://himalaya-tools.sourceforge.net/Mafia/). 

You may use our program for research, publication and education purposes, provided that you properly cite our paper in your publications.
In the following we illustrate how to use this package to summarize unweighted bipartite graph patterns, and weighted graph patterns:

Part I: Summarize unweighted bipartite graph patterns

Preprocessing 1: Get genes_to_phenotype.txt from Human Phenotype Ontology Consortium website: http://human-phenotype-ontology.org, and save it in directory "/MAFIA-MBC/Matlab-Preprocessing/". 
Preprocessing 2: Run: ./MAFIA-MBC/Matlab-Preprocessing/preprocessingg2f.m

The dataset is preprocessed into three parts:
(1) genes_to_phenotype.dat, in folder "MAFIA-MBC/datasets". This data is in transactional data format, where each line corresponds to a gene and each number in the line corresponds to a phenotype which the gene is related to.
(2) genes_to_phenotype_GeneList.txt, in folder "SumNetwork/datasets". Each line is a gene name of the corresponding line in the file "genes_to_phenotype.dat".  
(3) genes_to_phenotype_PhenotypeList.txt, in folder "SumNetwork/datasets". Each line is a phenotype name. The line number, starting from 1, is the corresponding phenotype ID in the file "genes_to_phenotype.dat".

!!!Important Notes:
(1) The genes_to_phenotype.txt file, including its content and format, has been updated from time to time by Human Phenotype Ontology Consortium. 
(2) The current release version has passed the test on the genes_to_phenotype.txt downloaded on 2/23/2012.  
	The file format may change such that in the future it is possible the provided Matlab-Preprocessing code cannot handle the file. We are not responsible for maintaining and updating our program.


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
Please download a latest version of GEO GeneCoexpression dataset for the following procedure. 

The dataset shall be preprocessed into the following parts:
(1) GeneCoexp.matrix. This gene coexpression data can be obtained by Pearson correlation, Spearman correlation or other methods.
(2) GeneCoexp_MAFIA_06.dat, saved in folder "MAFIA-MBC/datasets". This dataset is a result of converting "GeneCoexp.matrix" to a transacation data for MAFIA to process. As an example, here the edge cutoff threshold is set to be 0.6. 
(3) GeneCoexp_genelist.txt, saved in folder "SumNetwork/datasets". Each line is a gene name of the corresponding row (or column) in the file "GeneCoexp.matrix". The line number, starting from 1, is the corresponding gene ID in the file "GeneCoexp_MAFIA_06.dat".

Step 1: Get Closed Frequent Itemset:
In directory "MAFIA-MBC/", run:
qsub script_file_mafia_GeneCoexp_06_01

Step 2: Get cliques based on the closed frequent itemset:
In directory "MAFIA-MBC/", run:
qsub script_file_MBC_GeneCoexp_06_clique_01

Step 3: (First round) Merge maximal bicliques using MultiMerge strategy with density threshold being 0.7:
In directory "SumNetwork/", run:
qsub script_file_GraphMining_g2p_S_07

Step 4: (Second round, if needed) Merge maximal bicliques using SingleMerge strategy with density threshold being 0.7:
In directory "SumNetwork/", run:
qsub script_file_GraphMining_g2p_S_07_2nd

Step 5: Post processing the merging results into gene list:


Part III: Additional Hints:
(1) The best way to understand the above procedue is to read our paper.
(2) If your machine does not have a queuing system, please open the script file and use the command in the file.
(3) Simply running the executable file, MBC, GraphMining, and you will see the command line formats.
(4) Our program also supports summarizing weighted bipartite network patterns, and unweighted graph patterns, if properly setting up the workflow.


