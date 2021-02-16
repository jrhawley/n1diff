# Using _a priori_ information to refine search space

In our estimates of small sample sizes and how our estimates compare to the fully-balanced dataset of 48 biological replicates, we found a decreased mean-square error (MSE) in the fold change estimates compared to the OLS estimates.
However, the gains were marginal, only a 1.24% decrease in MSE.

This was due to the large shrinkage coefficient (i.e. it was often close to 1).
This provided minimal shrinkage of fold change estimates, and resulted in marginal gains.

In this investigation, we want to see whether using _a priori_ information to reduce the search space of the transcripts being considered for differential expression, can lead to greater reductions in MSE.

## Materials

We use the same RNA-seq data as before, from [Gierlinksi _et al._, Bioinformatics, 2015](https://academic.oup.com/bioinformatics/article/31/22/3625/240923).

Instead of investigating the entire transcriptome, we can use information about the mutant itself.
The _Snf2_ protein is part of the SNI/SNF complex in yeast.
Given the protein's function and the importance of chromatin remodelling in cellular processes, we can focus on the other members of the SWI/SNF complex in yeast.
They are as follows [^1]:

- ARP4
- RSC1
- SWI1
- SWI2
- SWI3
- SNF2
- SNF5
- SNF12

## Methods

To investigate the effect of smaller search spaces we perform the same differntial analyses as before.
We propose two experiments.
The first only considers the transcripts in the SWI/SNF complex, to assess how _a priori_ biological information can help improve fold-change estimation.

The second considers smaller random sets of transcripts of various sizes.
We consider transcript sets of sizes {3, 10, 50, 100, 500, 1000, and 5817}.
The smallest size of 3 is chosen since it is the smallest set required for the James-Stein estimates to theoretical achieve a lower MSE than the OLS estimate.
The remaining sizes are to assess the benefits across a variety of applications.
The largest size of 5817 transcripts is the entire transcriptome.

Samples are randomly selected 30 times in both experiments.
To match realistic RNA-seq experiments, 6 samples in total are selected in each iteration.
Half of the iterations contain 1 WT and 5 mutant samples, and the other half contain 5 WT and 1 mutant samples.
For the second experiment, transcripts are randomly selected from the filtered transcript set of the Sleuth objects.
This sampling is performed once per sample-level selection, so that each iteration contains a randomly selected set of samples and a random, indepedent, set of transcripts.

Code for these experiments can be found in [`swi-snf-sleuth.R`](swi-snf-sleuth.R) and [`random-sleuth.R`](random-sleuth.R), respectively.

## Results

[^1]: Members of the SWI/SNF complex in yeast are obtained from [Table 1](https://www.nature.com/articles/onc20094/tables/1) in [Reisman, Glaros, and Thompson, Oncogene, 2009](https://doi.org/10.1038/onc.2009.4).
