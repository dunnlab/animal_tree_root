
Rooting the animal tree of life
===============================

Introduction
------------

There has long been uncertainty about the relationship of Ctenophora to other animals \[1\].

Defining Porifera-sister and Ctenophora-sister... XXX Explain that it is about support for calde of all animals except Ctenophora or Porifera...

Variation across studies
------------------------

### Models of molecular evolution

Models of molecular evolution have several components that each consider different aspects of the evolutionary process. The models that have been used to model protein evolution in studies of the animal root have largely differed according to three components: the exchangability matrix *E*, the rate of evolution, and the state equilibrium frequencies *Π*.

The Exchangeability matrix *E* describes the rate at which one amino acid changes to another. Exchangability matrices have been used in the studies under consideration here include:

-   F81 \[2\] corresponds to equal rates between all states. The F81 matrix is also sometimes referred to as the Poisson matrix. It has no free parameters to estimate since all off-diagonal elements are set to 1.

-   WAG \[3\] is an emprically derived exchangability matrix based on a dataset of 182 globular protein families. It has no free parameters to estimate since all off-diagonal elements are set according to values estimated from this particular sample dataset.

-   LG \[4\], like WAG, is an emprically derived exchangability matrix. It is based on a much larger set of genes, and variation in rates across sites was taken into consideration when it was calculated. It has no free parameters to estimate since all off-diagonal elements are set according to values estimated from this particular sample dataset.

-   GTR, the General Time Reversible exchangability matrix, has free parameters for all off-diagonal elements that describe the exchangability of different amino acids. It is consrained so that changes are reversible, *i.e.* the rates above the diagonal are the same as those below the diagonal. This leaves 190 parameters that must be estimated from the data long with the other model parameters and the phylogenetic tree topology. This estimation requires a considerable amount of data and computational power, but if successful has the advantage of being based on the dataset at hand rather than a different dataset (as for LG and WAG).

While the exchancheability matrix describes the relative rate of different changes between amino acids, the actual rate can be further scaled. There are couple approaches that have been used in the studies considered here:

-   Site homogenous rates. The rates of evolution are assumed to be the same at all sites in the amino acid alignment.

-   Gamma rate heterogeneity. Each site is assigned to a different rate class with its own rate value. This accommodates different rates of evolution across different sites. Gamma is used so commonly that sometimes it isn't even specified, making it difficult at times to know if a study uses Gamma or not.

The vector of equilibrium frequencies *Π* describes the stationary frequency of amino acids. There are a few approaches that have been used across the studies considered here:

-   Empirical site homogenous. The frequency of each amino acid is observed from the matrix under consideration and applied to homogenously to all sites in the matrix.

-   Estimated site homogenous. The frequency of each amino acid is inferred along with other model parameters, under the assumption that it is the same at all sites.

-   CAT site heterogeneous \[5\]. Each site is assigned to a class with its own equilibrium frequencies. The number of classes, assignment of sites to classes, and equilibrium frequencies within the data are all estimated in a Bayesian framework.

Models can be assembled by selecting different options for all these different components. The models that are applied in practice area heavily influenced by engineering and computational costs, as well as convention. For example, on the questions considered here F81 and GTR exchangability matrices have only been used in combination with CAT site heterogeneous models of equilibrium frequency. LG and WAG exchangability matrices have only been used with site homogenous estimates of equilibrium frequency. This is further confused by the abreviations that are used for models. Papers often discuss CAT and WAG models as if they are exclusive, but these particular terms apply to non-exclusive model components-- CAT refers to variation across sites and WAG a particular exchangability matrix. CAT is generally shorthand for F81+CAT and WAG is shorthand for WAG+homogenous equilibrium frequency estimation. One could, though, run a WAG+CAT model.

To avoid confusion on this point, we always specify the exchangability matrix first, followed by modifiers that describe accomodation of heterogeneity in equilibrium frequencies (*e.g.*, CAT) or rate (*e.g.*, Gamma). If there are no modifiers, then it is implied that site homogenous models are used.

### Gene sampling

### Outgroup taxon sampling

Choanimalia, Holozoa, Opisthokonta

### Ingroup taxon sampling

Sensitivity to ingroup sampling has recieved less attention than sensitivity to outgroup sampling. This may be because results have tended to be more sensitive to outgroup sampling.

Overview of published analyses
------------------------------

![](manuscript_files/figure-markdown_github/support-1.png)

A total of 134 analyses were transcribed from the literature.

Details of published analyses
-----------------------------

### Dunn *et al.* 2008

Dunn *et al.* \[6\] added Expressed Sequence Tag (EST) data for 29 animals. It was the first phylogenomic analysis that included ctenophores, and therefore that could test the relationships of both Ctenophora and Porifera to the rest of animals. It was the first phylogenetic analysis to recover Ctenophora as the sister group to all other animals.

The data matrix was constructed using a semi-automated approach. Genes were translated into proteins, promiscuous domains were masked, all gene sequences from all species were compared to each other with blastp, genes were clustered based on this similarity with TribeMCL \[7\], and these clusters were filtered to remove those with poor taxon sampling and high rates of lineage-specific duplciations. Gene trees were then constructed, and in clades of sequences all from the same species all but one sequence were removed (these groups are often due to assembly errors). The remaining gene trees with more than one sequence for any taxon were then manually inspected. If strongly supported deep nodes indicative of paralogy were found, the entire gene was discarded. If the duplications for a a small number of taxa were unresolved, all genes from those taxa were excluded. Genes were then realigned and sites were filtered with Gblocks \[8\], resulting in a 77 taxon matrix. Some taxa in this matrix were quite unstable, which obscured other strongly-supported relationships. Unstable taxa were identified with leaf stability indeces \[9\], as implemented in phyutility \[10\], and removed from the matrix. This resulted in the 64-taxon matrix that is the focus of most of their analyses. Phylogenetic analyses were conducted under the F81+CAT model in phylobayes \[5\], and under the WAG model in MrBayes \[11\] and RAxML \[12\].

Regarding the recovery of Ctenophora-sister, the authors concluded:

> The placement of ctenophores (comb jellies) as the sister group to all other sampled metazoans is strongly supported in all our analyses. This result, which has not been postulated before, should be viewed as provisional until more data are considered from placozoans and additional sponges.

Note that there was, in fact, an exception to strong support. An analysis of the 40 ribosomal proteins in the matrix recovered Ctenophora-sister with only 69% support. This study did not include *Trichoplax*.

### Philippe *et al.* 2009

Philippe *et al.* 2009 \[Philippe:2009hh\]...

### Hejnol et al. 2009

### Pick *et al.* 2010

Pick *et al.* \[13\] sought to test wether Ctenophora-sister was an artefact of insufficient taxon sampling. They added new and additional published sequence data to the 64-taxon matrix of Dunn *et al.* \[6\]. The new taxa included 12 sponges, 1 ctenophore, 5 cnidarians, and *Trichoplax*. They firther modified the matrix by removing 2,150 sites that were poorly sampled or aligned. They considered two different sets of outgroups: Choanoflagellata (resulting in Choanimalia) and the same sampling as Dunn *et al.* (resulting in Opisthokonta).

All their analyses were conducted under the F81+CAT+Gamma model in phylobayes \[5\], in both a bayesian framework and with bootstrapping. All analyses have the same ingroup sampling and site removal so it isn't possible to independently assess the impact of these factors. Analyses with Choanimalia sampling recovered Porifera-sister with 72% posterior probability (PP) and 91% bootstrap support (BS). With broader Opisthokonta sampling, support for Porifera-sister is 84% PP. This is an inveresting case where increased outgroup sampling leads to increased support for Porifera-sister.

The authors argue that previous results supporting Ctenophora-sister "are artifacts stemming from insufficient taxon sampling and long-branch attraction (LBA)" and that "this hypothesis should be rejected". Thout the posterior probabilities supporting Porifera-sister are not strong, they conclude:

> Results of our analyses indicate that sponges are the sister group to the remaining Metazoa, and Placozoa are sister group to the Bilateria

They also investigated saturation, and coclude that Dunn *et al.* \[6\] is more saturated than Philippe *et al.* 2009 \[Philippe:2009hh\]. Note that the Pick *et al.* \[13\] dataset is not reanalyzed here because partition data are not available, and due to site filtering the partition file from Dunn *et al.* \[6\] cannot be applied to this matrix.

References
----------

1. Wallberg A, Thollesson M, Farris J, Jondelius U. The phylogenetic position of the comb jellies (Ctenophora) and the importance of taxonomic sampling. Cladistics. 2004;20: 558–578. Available: <http://onlinelibrary.wiley.com/doi/10.1111/j.1096-0031.2004.00041.x/full>

2. Felsenstein J. Evolutionary trees from DNA sequences: a maximum likelihood approach. Journal of Molecular Evolution. 1981;17: 368–376. Available: <http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=7288891&retmode=ref&cmd=prlinks>

3. Whelan S, Goldman N. A General Empirical Model of Protein Evolution Derived from Multiple Protein Families Using a Maximum-Likelihood Approach. Molecular Biology and Evolution. 2001;18: 691–699. doi:[10.1093/oxfordjournals.molbev.a003851](https://doi.org/10.1093/oxfordjournals.molbev.a003851)

4. Le SQ, Gascuel O. An improved general amino acid replacement matrix. Molecular Biology and Evolution. 2008;25: 1307–1320. doi:[10.1093/molbev/msn067](https://doi.org/10.1093/molbev/msn067)

5. Lartillot N. A Bayesian Mixture Model for Across-Site Heterogeneities in the Amino-Acid Replacement Process. Molecular Biology and Evolution. 2004;21: 1095–1109. doi:[10.1093/molbev/msh112](https://doi.org/10.1093/molbev/msh112)

6. Dunn CW, Hejnol A, Matus DQ, Pang K, Browne WE, Smith SA, et al. Broad phylogenomic sampling improves resolution of the animal tree of life. Nature. 2008;452: 745–749. doi:[10.1038/nature06614](https://doi.org/10.1038/nature06614)

7. Enright A, Van Dongen S, Ouzounis C. An efficient algorithm for large-scale detection of protein families. Nucleic Acids Research. Oxford University Press; 2002;30: 1575–1584. doi:[10.1093/nar/30.7.1575](https://doi.org/10.1093/nar/30.7.1575)

8. Castresana J. Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular Biology and Evolution. 2000;17: 540–552. Available: [http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Citation&list\_uids=10742046 ](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Citation&list_uids=10742046 )

9. Thorley J, Wilkinson M. Testing the phylogenetic stability of early tetrapods. Journal of Theoretical Biology. 1999;200: 343–344. doi:[10.1006/jtbi.1999.0999](https://doi.org/10.1006/jtbi.1999.0999)

10. Smith SA, Dunn CW. Phyutility: a phyloinformatics tool for trees, alignments and molecular data. Bioinformatics. Oxford University Press; 2008;24: 715–716. doi:[10.1093/bioinformatics/btm619](https://doi.org/10.1093/bioinformatics/btm619)

11. Ronquist F, Huelsenbeck JP. MrBayes 3: Bayesian phylogenetic inference under mixed models. Bioinformatics. 2003;19: 1572–1574. doi:[10.1093/bioinformatics/btg180](https://doi.org/10.1093/bioinformatics/btg180)

12. Stamatakis A. RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics. 2006;22: 2688–2690. doi:[10.1093/bioinformatics/btl446](https://doi.org/10.1093/bioinformatics/btl446)

13. Pick KS, Philippe H, Schreiber F, Erpenbeck D, Jackson DJ, Wrede P, et al. Improved phylogenomic taxon sampling noticeably affects nonbilaterian relationships. Molecular Biology and Evolution. 2010;27: 1983–1987. doi:[10.1093/molbev/msq089](https://doi.org/10.1093/molbev/msq089)
