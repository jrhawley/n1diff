# Smaller sample sizes

This folder compares the James-Stein estimates to smaller samples sizes of a well-controlled differential expression study.
We compare the James-Stein estimates to the OLS estimates of an unbalanced experimental design, as well as the OLS estimates of a similarly sized, but balanced, experimental design.

## Materials

We use the RNA-seq data from Gierlinksi et al, Bioinformatics, 2015, which has 48 biological replicates of _Saccharomyces cerevisiae_ cells.
The cells are either mutated to remove the _snf2_ gene on chrXV, or wild type.
The full, 48 replicate, differential analysis results can be found in [`../../data/Gierlinski_2015/`](../../data/Gierlinski_2015/).

## Methods

We select $$N$$ total samples in an experiment, where $$N \in {4, 6, ..., 22, 24}$$ and compare the differential expression results of an evenly balanced experimental design (i.e. $$\frac{N}{2}$$ wild type and $$\frac{N}{2}$$ mutated samples), an unbalanced experimental design (i.e. 1 mutated and $$N - 1$$ wild type samples, or vice versa), and an unbalanced experimental design with the James-Stein estimates.

A key parameter in the James-Stein estimates is the selection of which transcripts to consider, whereas differential expression analysis is typically concerned with all transcripts.
For these results, we use all transcripts.
Determining the optimal number of transcripts is the subject of another analysis.

To assess the stability of the results for each method, we randomly select $$N$$ samples, according to the experimental design, 30 times, and calculate the true positives, false positives, true negatives, and false negatives for each method by comparing them to the results from the full 48 replicate experiment.

Methods are compared with a one-way ANOVA test, with a significance threshold of 0.05.
The methods are compared with the following model:

$$
S \sim \text{Method} + N
$$

where $$S$$ is the statistic of interest, such as TPR or ACC, and $$N$$ is the total number of samples in the comparison.

## Results

Largely, we find that the James-Stein estimates produce fewer false negatives and more true positives than the OLS estimators.
This translates to slight but significantly higher true positive rate, negative predictive value, higher accuracy, and higher Matthews correlation coefficient.
We find on average, the following changes to statistics in the James-Stein estimates compared to the OLS esimates:

| Statistic | $$\beta_{JS}$$       | pval                 | qval                 |
| --------- | -------------------- | -------------------- | -------------------- |
| TPR       | 0.030354838205541    | 0.000192480705728893 | 0.000451528670412287 |
| TNR       | -0.00575908608351453 | 0.148108697610896    | 0.148108697610896    |
| PPV       | -0.0078619449606338  | 0.0953925069121369   | 0.114471008294564    |
| NPV       | 0.00753139140269705  | 0.000301019113608192 | 0.000451528670412287 |
| ACC       | 0.016912601960141    | 0.000143852644329807 | 0.000451528670412287 |
| MCC       | 0.0232626079808969   | 0.00026665173474846  | 0.000451528670412287 |

When comparing all 5817 transcripts in the yeast genome, we find an average increase of 3% in the TPR in the James-Stein estimates, compared to the standard OLS estimates.

## Conclusions

The James-Stein estimates have a higher TPR across sample sizes than OLS estimates, demonstrating its utility in unbalanced experimental designs.
