# Canonical differential expression analysis

Consider the statistical model employed in the `sleuth` R package for differential gene expression analysis.
It is described as follows:

Consider a set of transcripts, $T$, measured in $n$ samples with an experimental design matrix $X \in \mathbb{R}^{n \times p}$.
Let $Y_{ti}$ be the natural log of the abundance of transcript $t$ in sample $i$.
Given the design matrix, $X = [x_1^T; x_2^T; ... x_n^T], x_i \in \mathbb{R}^p$, we model the abundance of transcripts as

$$
Y_{ti} = x_i^T \beta_t + \epsilon_{ti}
$$

where $\epsilon_{ti} \sim N(0, \sigma_t^2)$ is the biological noise of transcript $t$ in sample $i$ and $\beta_t \in \mathbb{R}^p$ is the fixed effect of the covariates on the expression of transcript $t$.
Due to inferential noise from sequencing, the $Y_{ti}$'s are not observed directly, but indirectly through the observed perturbations $D_{ti}$.
This can be modelled as

$$
D_{ti} | Y_{ti} = Y_{ti} + \zeta_{ti}
$$

where $\zeta_{ti} \sim N(0, \tau_t^2)$ is the inferential noise of transcript $t$ in sample $i$.
Both biological and inferential noise are each IID and independent of each other, namely

$$
Cov[\epsilon_{ti}, \epsilon_{sj}] = \sigma_t^2\delta_{i,j}\delta_{t,s} \\
Cov[\zeta_{ti}, \zeta_{sj}] = \tau_t^2\delta_{i,j}\delta_{t,s} \\
Cov[\epsilon_{ti}, \zeta_{sj}] = 0 \forall s,t \forall i,j
$$

Thus, our transcript abundances are modelled as

$$
D_{t} | Y_{t} \sim N_p(X\beta_t, (\sigma_t^2 + \tau_t^2)I_n)
$$

where $I_n \in \mathbb{R}^{n \times n}$ is the identity matrix.
The goal of the differential analysis is to estimate the $|T| \cdot p$ coefficients in $\beta_t \forall t \in T$, and to determine which coefficients differ significantly from 0.
This is achieved through a Wald test or likelihood ratio test after estimating the inferential variance, $\tau_t^2$, through bootstrapping and the biological variance, $\sigma_t^2$, through dispersion estimation and shrinkage.

The estimator for the differential effect, $\hat{\beta_t}$, is obtained through the ordinary least squares estimate

$$
\hat{\beta_t} = (X^TX)^{-1}X^Td_t
$$

where $d_t$ is the observed abundances given by

$$
d_{ti} = \ln \left(\frac{k_{ti}}{\hat{s_i}} + 0.5 \right) \\
\hat{s_i} = \underset{t \in T^*}{\mathrm{median}} \frac{k_{ti}}{\sqrt[n]{\underset{j = 1}{\overset{n}{\prod }}k_{tj}}}
$$

where $k_{ti}$ is the estimated read count from `kallisto` for transcript $t$ in sample $i$ and $\hat{s_i}$ is the scaling factor for sample $i$, calculated from the set of all transcripts that pass initial filtering, $T^*$.

## Bias and variance of the OLS estimator

As shown in Supplementary Note 2 of [the `sleuth` paper](https://doi.org/10.1038/nmeth.4324), the OLS estimator is unbiased $\left( \mathbb{E} \left[ \hat{\beta_t}^{(OLS)} \right] = \beta_t \right)$.
It can also be shown that $\text{Var} \left[ \hat{\beta_t}^{(OLS)} \right] = (X^TX)^{-1} \Sigma X (X^TX)^{-1}$ for a general covariance matrix $\Sigma$.
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
