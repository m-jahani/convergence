# convergence

## TOPCAN.R

Identification of top candidate windows and testing window_by_window signatures of convergence (null-W test) with considering recombination effect.

## LDCLUST.R

to cluster convergent windows based on their LD and genetic distance

The clustering procedure conducts a stepwise process as follows: at the initial step, it takes the first convergent window from each chromosome and estimates the Linkage Disequilibrium (LD) with the second convergent window. At the following step, the procedure moves a window forward and calculates LD between the second and third windows. The process continues moving forward step by step as long as these conditions are met in each step; 1) the calculated LD is significantly different (P < .05) from the null distribution, and 2) the distance between the two windows is not more than 5 cM (editable). Immediately after a step has not met either of the two conditions, a cluster forms with all previous steps’ windows. These windows are then subsetted from the data, and the clustering algorithm starts from the first step with the remaining windows. This process continues until all the convergent windows have been clustered.
In each step of the procedure, PLINK (reference) was used to calculate LD (r2) among all pairwise combinations of the SNPs within the two windows. Later on, all the calculated r2 values were summarized by averaging to estimate the LD level among the two windows of the step. 
The null distribution in each step was constructed with LD estimations among 10,000 random, non-convergent window pairs.
