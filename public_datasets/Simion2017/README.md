SuppData_Metazoa_2017
---

This repository contains supplemental materials from :

Paul Simion, Hervé Philippe, Denis Baurain, Muriel Jager, Daniel J. Richter, Arnaud Di Franco, Béatrice Roure, Nori Satoh, Éric Quéinnec, Alexander Ereskovsky, Pascal Lapébie, Erwan Corre, Frédéric Delsuc, Nicole King, Gert Wörheide, Michaël Manuel, A Large and Consistent Phylogenomic Dataset Supports Sponges as the Sister Group to All Other Animals, Current Biology, Available online 16 March 2017, ISSN 0960-9822, http://dx.doi.org/10.1016/j.cub.2017.02.031.


## Dataset

**alignements_1719genes.tgz**  
This archive contains the 1719 alignments as obtained after our complete dataset building protocol.

**supermatrix_97sp_401632pos_1719genes.fasta**  
this file corresponds to the concatenation of the entire dataset (1719 gene alignements).

**supermatrix_90sp_136618pos_heterop60.puz**  
Dataset after the removal of the 60% most heteropecileous sites.

**supermatrix_90sp_102464pos_heterop70.puz**  
Dataset after the removal of the 70% most heteropecileous sites.

**supermatrix_90sp_268032pos_dayhoff6.phy**  
Dayhoff6 recoding of our dataset (constant sites were removed).
You probably want to remove the first line of this file when using other softwares than PhyloBayes.

**partition_401632pos_1719genes.part**  
Partitionning scheme corresponding to the boundaries of the 1719 genes.

**supermatrix_whelan2015_81006pos_NoDemo.phy**  
Dataset from Whelan et al. 2015 in which demosponges have been removed.

**supermatrix_whelan2015_81006pos_NoDemoCalcHomo.phy**  
Dataset from Whelan et al. 2015 in which all sponges except hexactinellides have been removed.

**partition_whelan2015_81006pos.part**  
Partitionning scheme corresponding to the boundaries of the 251 genes from Whelan et al. 2015.


## Trees

### complete dataset analyses

**tree_97sp_CAT.tre**  
10 jackknife replicates of 100,000 position each - CAT+G4 (Phylobayes)

**tree_90sp_CAT.tre**  
100 jackknife replicates of 100,000 position each - CAT+G4 (Phylobayes)

**tree_97sp_LGF-PARTITION.tre**  
100 bootstraps - partitionning by gene - LG+G4+F (RAxML)

**tree_90sp_LGF-PARTITION.tre**  
100 bootstraps - partitionning by gene - LG+G4+F (RAxML)

### removal of heteropecileous sites

**tree_90sp_CAT_heterop60.tre**  
2 independant MCMC chains - CAT+G4 (Phylobayes)

**tree_90sp_CAT_heterop70.tre**  
2 independant MCMC chains - CAT+G4 (Phylobayes)

### model comparison when reducing taxonomic sampling

**tree_NoDemo_CAT.tre**  
10 jackknife replicates of 100,000 position each - CAT+G4 (Phylobayes)

**tree_NoDemo_LGF-PARTITION.tre**  
100 bootstraps - partitionning by gene - LG+G4+F (RAxML)

**tree_NoDemoCalcHomo_CAT.tre**  
10 jackknife replicates of 100,000 position each - CAT+G4 (Phylobayes)

**tree_NoDemoCalcHomo_LGF-PARTITION.tre**  
100 bootstraps - partitionning by gene - LG+G4+F (RAxML)

**tree_Whelan2015_NoDemo_CAT.tre**  
2 independant MCMC chains - CAT+G4 (Phylobayes)

**tree_Whelan2015_NoDemo_LGF-PARTITION.tre**  
100 bootstraps - partitionning by gene - LG+G4+F (RAxML)

**tree_Whelan2015_NoDemoCalcHomo_CAT.tre**  
2 independant MCMC chains - CAT+G4 (Phylobayes)

**tree_Whelan2015_NoDemoCalcHomo_LGF-PARTITION.tre**  
100 bootstraps - partitionning by gene - LG+G4+F (RAxML)

### Pipeline cleaning example : the rpl2 gene

**tree_rpl2_FIFO.pdf**  
rpl2 tree when using initial datasets (right after the Filter Focus step - FIFO).

**tree_rpl2_DCC.pdf**  
rpl2 tree when using datasets after the De-Cross-Contamination step (DCC).

**tree_rpl2_DC1.pdf**  
rpl2 tree when using datasets after the De-Contamination step 1 (DC1).

**tree_rpl2_DC2.pdf**  
rpl2 tree when using datasets after the De-Contamination step 2 (DC2).

**tree_rpl2_DC3.pdf**  
rpl2 tree when using datasets after the De-Contamination step 3 (DC3).

**tree_rpl2_FINAL.pdf**  
Final rpl2 tree after all cleaning steps of our pipeline, as used for supermatrix concatenation.


## Softwares

**utilities_src.tgz**  
These are the sources of the C programs used in our dataset assembly procedure.
Please note that the "De-Cross-Contamination" (DCC) step of our procedure has been re-worked
as a dedicated software (named "CroCo") currently under final development and that will be published elsewhere.


## Information

**taxon_names.dict**  
This dictionnary indicates the corresponds between full taxon names and their short version found in phylip files.

**choanoflagellate names**  
Choanoflagellate names used in our study corresponds to the valid names recently published by Carr et al. ([2016](http://www.sciencedirect.com/science/article/pii/S1055790316302743)).


## Post-Publication Materials
The *Post-Publication_Materials_following_arising_questions* folder contains additional files and analyses details that we made available following arising technical questions.

**Post-publication_Note.pdf**  
Additional note on the MCMC convergence approximation strategy used and on the further evaluation of chain convergence.

**100_jackknife_replicates_90sp_CAT.trees**  
All jackknife trees corresponding to the 100 replicates used to build *tree_90sp_CAT.tre* - CAT+G4 (Phylobayes)

**10_jackknife_replicates_97sp_CAT.trees**  
All jackknife trees corresponding to the 10 replicates used to build *tree_97sp_CAT.tre* - CAT+G4 (Phylobayes)

**10_jackknife_replicates_80sp_CAT.trees**  
All jackknife trees corresponding to the 10 replicates used to build *tree_NoDemo_CAT.tre* - CAT+G4 (Phylobayes)

**10_jackknife_replicates_70sp_CAT.trees**  
All jackknife trees corresponding to the 10 replicates used to build *tree_NoDemoCalcHomo_CAT.tre* - CAT+G4 (Phylobayes)
