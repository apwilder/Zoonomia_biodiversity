# Zoonomia_biodiversity
Scripts for data generation and analysis of genomic correlates of extinction risk across 240 mammal species
Scripts for data generation and analysis of genomic correlates of extinction risk across 240 mammal species in the Zoonomia multispecies genome alignment.

Below is a brief description of the scripts used for data generation and analysis.

1.	CallBranchMutations.sh
For a given species in the HAL alignment, runs the halBranchMutations and halLiftover functions from the Comparative Genomics Toolkit to call derived substitutions relative to the closest reconstructed ancestral node, then assigns phylop scores for substitution positions lifted to human autosomes using bedmap. PhyloP scores were estimated using Phast v 1.5 as described by Christmas et al. (accepted) Science.

2.	IDandFilterHypervariableWindows.R
For each species, estimates the distribution of counts of substitutions in 1KB windows across the genome, then identifies outlier windows with more substitutions than expected from a Poisson distribution using the R package alphaOutlier.

3.	Apply_HypervariableFilter.sh
Removes identified outlier windows from bed files of substitutions.

4.	SummarizePhyloP_fromBranchMutations.R
Compiles phyloP stats (e.g. mean phyloP score, proportion of substitutions at conserved sites, kurtosis of phylop distribution) across species into a matrix

5.	AnnotateCodingMutations.sh
Builds databases from gtf annotations for each genome for annotating variants in SnpEff. Annotates heterozygous variants from short read data and substitutions from HAL multispecies alignment. Combines info on phyloP of coding sites.

6.	50KBwindows_Phylop.sh
Takes 50kb windows with mean heterozygosity and ROH status (ROH=0, non-ROH=1) lifted to human genome, adds phylop scores and summarizes all heterozygosity, ROH and phylop stats within windows.

7.	50KBwindows_Coding.sh 
Summarizes coding variants and phylop in 50kb windows lifted to the human genome (by Ross Swofford @ Broad Inst.)

8.	GetSingleCopyBUSCOgenes.sh
Gets a list of mammalian BUSCO genes for restricting analysis to single copy BUSCO genes.

9.	FilterSingleCopyBUSCOgenes.R
With GetBUSCOgenenames.sh, compiles list of single copy busco genes.

10.	Summarize_HeterozygVars_snpeff_buscogenes_IMPC.R
Combine and summarize coding annotations for heterozygous sites across IMPC-annotated BUSCO genes for each species.

11.	Summarize_HomozygVars_snpeff_buscogenes_IMPC.R
Combine and summarize coding annotations for substitutions across IMPC-annotated BUSCO genes for each species.

12.	AdjustHomLoadStats4Mutations.R
Compiles datatables with genome-wide heterozygosity, ROH, Ne and homozygous load statistics and adjusts homozygous load statistics for branch lengths.

13.	CombineHomHetCodingMutations.R
Summarizes homozygous and heterozygous mutations in protein-coding genes with associated phenotypes categorized by IMPC.

14.	Make50KBmatrices4MLmodels.R
Compiles 50KB window-based stats for input to machine learning models

15.	Regressions_Main_text.R
Regression analyses reported in main text and supplement. 

