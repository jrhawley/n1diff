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
