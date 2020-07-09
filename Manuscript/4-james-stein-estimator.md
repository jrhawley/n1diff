# James-Stein Estimator

Consider the following scenario. For an $m$-variate normal distribution $Z \sim N_m(\mu, \Sigma)$ where $\mu$ is unknown and $\Sigma$ is known, if we observe a single realization of this distribution, $z$, what is a good estimator, $\hat{\mu}$, for the unknown mean $\mu$?
Stein \Cref{https://projecteuclid.org/euclid.bsmsp/1200501656} showed that the spherically symmetric estimator $\hat{\mu}_{JS} = \left( 1 - \frac{m - 2}{\Vert z \Vert^2} \right)z$ dominates the naive estimator $\hat{\mu}_{naive} = z$ for any mean $\mu$ in the case that $m \ge 3$ and $\Sigma = I_n$.
Theorem 2 of Bock \Cref{https://projecteuclid.org/download/pdf_1/euclid.aos/1176343009} extended this result to the general case to show that when $m \ge 3$, $\text{Tr}(\Sigma) \ge 2 \lambda_L$, and $0 \le c \le 2 \left( \frac{\text{Tr}(\Sigma)}{\lambda_L} - 2 \right)$

$$
\hat{\mu}_{JS} = \left( 1 - \frac{c}{z^T \Sigma^{-1} z} \right) z
$$

is the minimax estimator for $\mu$ of the mean square error $\mathbb{E} \left[ \Vert \hat{\mu} - \mu \Vert ^2\right]$ where $\lambda_L$ is the largest eigenvalue of $\Sigma$.
In the case that either $m < 3$ or $\text{Tr}(\Sigma) < 2 \lambda_L$, the naive estimator is then the minimax estimator.
The most shrinkage occurs when $c = 2 \left( \frac{\text{Tr}(\Sigma)}{\lambda_L} - 2\right)$, and the naive estimator is equivalent to the James-Stein estimator with $c = 0$.

## Applying the James-Stein estimator to differential gene expression

### A simple experimental design

Let's consider the `sleuth` model with our simple experimental design (`d ~ 1 + mutation`):

$$
D_t | Y_t \sim N_p \left( \beta_{t,0} + \mathbb{I}_{mut}\beta_{t,1}, (\sigma_t^2 + \tau_t^2)I_n \right)
$$

For our $n_{nonmut}$ non-mutated samples this is equivalent to

$$
D_t | Y_t \sim N_{n_{nonmut}} \left( \beta_{t,0}, (\sigma_t^2 + \tau_t^2)I_n \right)
$$

which can be fit with the same model process that `sleuth` typically employs.
For the single mutated sample, the model is

$$
D_t | Y_t \sim N \left( \beta_{t, 0} + \beta_{t, 1}, \sigma_t^2 + \tau_t^2 \right)
$$

The covariance matrix is the same as the mutated samples, but the mean $\beta_{t, 0} + \beta_{t, 1}$ is unknown and we have a single observation of this distribution.
If we are interested in a subset of all transcripts, $S \subset T$ (e.g. we are only concerned with transcripts near the mutated site), then we can consider the $|S|$-dimensional normal distribution and use the mean and variance estimates from the non-mutated samples.
Our model for the single mutated sample is then

$$
\Delta \sim N_{|S|}(\Beta_0 + \Beta_1, \Sigma) \\
$$

where

$$
\Beta_{i,s} = \beta_{s,i} \forall s \in S
$$

$$
\Sigma = \text{diag} \left( \max\{ \hat{\sigma}_1^2, \tilde{\sigma}_1^2 \} + \hat{\tau}_1^2, ..., \max\{ \hat{\sigma}_{|S|}^2, \tilde{\sigma}_{|S|}^2 \} + \hat{\tau}_{|S|}^2 \right) \\
$$

We switch from using coefficients $\beta_{t,i}$ to $\Beta_{i,s}$ to avoid confusion, since $\beta_{t,i} \in \mathbb{R}^p$ (a $p$-dimensional vector for each covariate in the design) whereas $\Beta_{i,s} \in \mathbb{R}^{|S|}$ (an $|S|$-dimensional vector for only a single coefficient over all transcripts in $S$).

With a single observation of this distribution, $\delta$, we can design a James-Stein estimator for the unknown effect coefficient, $\Beta_1$, that is shrunk towards 0.

$$
\hat{\Beta}_1^{(JS)} = \left( 1 - \frac{c}{(\delta - \hat{\Beta}_0)^T \Sigma^{-1} (\delta - \hat{\Beta}_0)} \right)(\delta - \hat{\Beta}_0)
$$

where $\hat{\Beta}_0$ is the estimate obtained from the non-mutated samples for all transcripts $s \in S$.

It is simple to see that $\text{Tr}(\Sigma) = \sum_{s \in S} \max\{ \hat{\sigma}_s^2, \tilde{\sigma}_s^2 \} + \hat{\tau}_s^2$ and that $\lambda_L = \max_{s \in S} \left\{ \max\{ \hat{\sigma}_s^2, \tilde{\sigma}_s^2 \} + \hat{\tau}_s^2 \right\}$.
In the case of $\text{Tr}(\Sigma) < 2 \lambda_L$ or $|S| \le 2$, we resort to the OLS estimator, which is found through the standard `sleuth` process.

### Comparison of the simple design with the OLS estimator

For this simple experimental design, the OLS estimator is given by

$$
\begin{bmatrix}
\hat{\beta}_{s,0}^{(OLS)} \\
\hat{\beta}_{s,1}^{(OLS)}
\end{bmatrix}
 = \hat{\beta}_s^{(OLS)}
 = (X^TX)^{-1}X^T d_s
 = \begin{bmatrix}
\bar{d}_s^{(nonmut)} \\
d_s^{(mut)} - \bar{d}_s^{(nonmut)}
\end{bmatrix}
$$

For our mutation coefficient, $\beta_{s,1}$, we get an OLS estimator of

$$
\hat{\beta}_{s,1}^{(OLS)} = d_s^{(mut)} - \hat{\beta}_{s,0}^{(OLS)} = \delta_s - \hat{\Beta}_{0,s}
$$

Thus, our James-Stein estimator for $\Beta_1$ is then

$$
\hat{\Beta}_1^{(JS)} = \left( 1 - \frac{c}{\left( \hat{\Beta}_1^{(OLS)} \right)^T \Sigma^{-1} \hat{\Beta}_1^{(OLS)}} \right) \hat{\Beta}_1^{(OLS)}
$$

The James-Stein estimate is parallel with the OLS estimate (which is already calculated by `sleuth`), but shrunk towards 0.

### Extending to complex experimental designs

For a more general experimental design, we can extend the above.
Given an experimental design matrix, $X \in \mathbb{R}^{n \times p}$, where $n > p$, $\text{rank}(X) = p$ and $\text{rank}(X^*) = p - 1$ where $X^* \in \mathbb{R}^{(n - 1) \times p}$ is the same design matrix but with one sample removed, a James-Stein estimator for the linear coefficient uniquely specified by the one sample is given by

$$
\hat{\Beta}_i^{(JS)} = \left( 1 - \frac{c}{\left( \hat{\Beta}_i^{(OLS)} \right)^T \Sigma^{-1} \hat{\Beta}_i^{(OLS)}} \right) \hat{\Beta}_i^{(OLS)}
$$

## Properties of the James-Stein estimator for differential expression

### Variance of the James-Stein estimator

This estimator is the minimax estimator under the mean square error, $\mathbb{E} \left[ \Vert \hat{\beta}_1^{(S)} - \beta_1^{(S)} \Vert ^2\right] = \sum_{s \in S} \mathbb{E}\left[ \left( \hat{\beta}_{1,s}^{(S)} - \beta_{1,s}^{(S)} \right)^2 \right] = \sum_{s \in S} \text{Var}\left[ \hat{\beta}_{1,s}^{(S)} \right]$.
While $\mathbb{E} \left[ \Vert \hat{\beta}_1^{(S)} - \beta_1^{(S)} \Vert ^2\right] \le \mathbb{E} \left[ \Vert \hat{\beta}_1^{(OLS)} - \beta_1^{(S)} \Vert ^2\right]$, this does not imply that $\text{Var}\left[ \hat{\beta}_{1,s}^{(S)} \right] \le \text{Var}\left[ \hat{\beta}_{1,s}^{(OLS)} \right] \forall s$.
Some transcripts may have larger variances than the naive estimator, but all transcripts in aggregate will have a smaller mean square error.
This is still desirable if the goal is to find if there is an effect on any transcripts in the set $S$, instead of a particular one within the set.
