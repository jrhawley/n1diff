# Canonical differential expression analysis

Consider the statistical model employed in the `sleuth` R package for differential gene expression analysis.
It is described as follows:

Consider a set of transcripts, $S$, measured in $N$ samples with an experimental design matrix $X \in \mathbb{R}^{N \times p}$.
Let $Y_{si}$ be the natural log of the abundance of transcript $s$ in sample $i$.
Given the design matrix, $X = [x_1^T; x_2^T; ... x_n^T], x_i \in \mathbb{R}^p$, we model the abundance of transcripts as

$$
Y_{si} = x_i^T B_s + \epsilon_{si}
$$

where $\epsilon_{si} \sim \mathcal{N}(0, \sigma_s^2)$ is the biological noise of transcript $s$ in sample $i$ and $B_s \in \mathbb{R}^p$ is the fixed effect of the covariates on the expression of transcript $s$.
Due to inferential noise from sequencing, the $Y_{si}$'s are not observed directly, but indirectly through the observed perturbations $D_{si}$.
This can be modelled as

$$
D_{si} | Y_{si} = Y_{si} + \zeta_{si}
$$

where $\zeta_{si} \sim \mathcal{N}(0, \tau_s^2)$ is the inferential noise of transcript $s$ in sample $i$.
Both biological and inferential noise are each IID and independent of each other, namely

$$
\mathbb{Cov}[\epsilon_{si}, \epsilon_{rj}] = \sigma_s^2\delta_{i,j}\delta_{s,r} \\
\mathbb{Cov}[\zeta_{si}, \zeta_{rj}] = \tau_s^2\delta_{i,j}\delta_{s,r} \\
\mathbb{Cov}[\epsilon_{si}, \zeta_{rj}] = 0 \\
\forall s,r \forall i,j
$$

Thus, our transcript abundances are modelled as a multivariate normal distribution

$$D_{s} | Y_{s} \sim \mathcal{N}_N(X B_s, (\sigma_t^2 + \tau_t^2)I_N)$$

where $I_N \in \mathbb{R}^{N \times N}$ is the identity matrix.
The goal of the differential analysis is to estimate the $|S| \cdot p$ coefficients in $B_s \forall s \in S$, and to determine which coefficients differ significantly from 0.
This is achieved through a Wald test or likelihood ratio test after estimating the inferential variance, $\tau_s^2$, through bootstrapping and the biological variance, $\sigma_s^2$, through dispersion estimation and shrinkage.

The estimator for the differential effect is the ordinary least squares estimate:

$$\hat{B}_s = (X^TX)^{-1}X^T d_s$$

where $d_s$ is the observed abundances given by

$$ d_{si} = \ln \left(\frac{k_{si}}{\hat{f_i}} + 0.5 \right) $$

$$ \hat{f_i} = \underset{s \in S^*}{\mathrm{median}} \frac{k_{si}}{\sqrt[N]{\underset{j = 1}{\overset{N}{\prod }}k_{sj}}} $$

where $k_{si}$ is the estimated read count from `kallisto` for transcript $s$ in sample $i$ and $\hat{f_i}$ is the scaling factor for sample $i$, calculated from the set of all transcripts that pass initial filtering, $S^*$.

## Bias and variance of the OLS estimator

As shown in Supplementary Note 2 of [the `sleuth` paper](https://doi.org/10.1038/nmeth.4324), the estimator is unbiased $\left( \mathbb{E} \left[ \hat{B}_s^{(OLS)} \right] = B_s \right)$.
It can also be shown that $\mathbb{V} \left[ \hat{B}_s^{(OLS)} \right] = (X^TX)^{-1} X^T \Sigma X (X^TX)^{-1}$ for a general covariance matrix $\Sigma$.
In the case where $\Sigma = (\sigma_t^2 + \tau_t^2)I_n$ (such as in the `sleuth` model), this reduces to $\text{Var} \left[ \hat{\beta_t}^{(OLS)} \right] = (\sigma_t^2 + \tau_t^2)(X^TX)^{-1}$.

Consider a simple experimental design where our only covariate of interest is the presence of a mutation.
Then our design matrix, with the first column being the intercept and the second being the mutation status, looks like so:

$$
X =
\begin{bmatrix}
1 & 1 \\
\vdots_{n_{mut}} & \vdots_{n_{mut}} \\
1 & 1 \\
1 & 0 \\
\vdots_{n_{nonmut}} & \vdots_{n_{nonmut}} \\
1 & 0
\end{bmatrix}
$$

The variance of the OLS estimator is then

$$
\text{Var} \left[ \hat{\beta_t}^{(OLS)} \right] = \frac{(\sigma_t^2 + \tau_t^2)}{n_{mut} n_{nonmut}}
\begin{bmatrix}
n_{mut} & -n_{mut} \\
- n_{mut} & n_{mut} + n_{nonmut} \\
\end{bmatrix}
$$

Importantly, the estimate for the coefficient measuring the effect that the presence of the mutation has variance $\frac{(\sigma_t^2 + \tau_t^2)(n_{mut} + n_{nonmut})}{n_{mut} n_{nonmut}}$.
When we only have 1 mutated sample, as per the motivation of this work, this reduces to $\frac{(\sigma_t^2 + \tau_t^2)(1 + n_{nonmut})}{n_{nonmut}}$.

## Desired properties of an alternative estimator

Our goal is to produce an alternative estimator whose theoretical variance is less than the value from the previous section.
Given the unstable and uncertain nature of estimation using a single observation, reducing variance by biasing the estimators towards 0 is a conservative approach that is appropriate for this setting.
Some shrinkage approaches like ridge and lasso regression approach this by placing a bound on the $L^1$ and $L^2$ norms of the coefficients.
In differential gene expression, the effect of a covariate on a single transcript is not necessarily related to that of another transcript, so fixing an upper bound on some $L^p$ norm of the effects is not desired since the magnitude of effects is not known ahead of time.

Moreover, this type of shrinkage can increase the estimated effect of some covariates while reducing others closer to zero.
If the design matrix $X$ is not orthogonal (which is true for our setting) then the shrunk estimate is not parallel to the OLS estimate, it is a linear transformation of it.
This can lead to over-estimating the mutation effect for some transcripts at the cost of others.
Shrinkage methods that shrink the entire estimate towards 0 thus may be more appropriate in this setting.
One example of a shrunk estimator with this property is the James-Stein estimator.
