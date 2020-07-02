# James-Stein Estimator

Consider the following scenario. For an $m$-variate normal distribution $Z \sim N_m(\mu, \Sigma)$ where $\mu$ is unknown and $\Sigma$ is known, if we observe a single realization of this distribution, $z$, what is a good estimator, $\hat{\mu}$, for the unknown mean $\mu$?
Stein \Cref{https://projecteuclid.org/euclid.bsmsp/1200501656} showed that the spherically symmetric estimator $\hat{\mu}_{JS} = \left( 1 - \frac{m - 2}{\Vert z \Vert^2} \right)z$ dominates the naive estimator $\hat{\mu}_{naive} = z$ for any mean $\mu$ in the case that $m \ge 3$ and $\Sigma = I_n$.
Theorem 2 of Bock \Cref{https://projecteuclid.org/download/pdf_1/euclid.aos/1176343009} extended this result to the general case to show that when $m \ge 3$, $\text{Tr}(\Sigma) \ge 2 \lambda_L$, and $0 \le c \le 2 \left( \frac{\text{Tr}(\Sigma)}{\lambda_L} - 2\right)$

$$
\hat{\mu}_{JS} = \left( 1 - \frac{c}{z \Sigma^{-1} z^T} \right) z
$$

is the minimax estimator for $\mu$ of the mean square error $\mathbb{E} \left[ \Vert \hat{\mu} - \mu \Vert ^2\right]$ where $\lambda_L$ is the largest eigenvalue of $\Sigma$.
In the case that either $m < 3$ or $\text{Tr}(\Sigma) < 2 \lambda_L$, the naive estimator is then the minimax estimator.
The most shrinkage occurs when $c = 2 \left( \frac{\text{Tr}(\Sigma)}{\lambda_L} - 2\right)$, and the naive estimator is equivalent to the James-Stein estimator with $c = 0$.

## Applying the James-Stein estimator to different gene expression

Let's consider the `sleuth` model with our simple experimental design:

$$
D_t | Y_t \sim N_p \left( X\beta_t, (\sigma_t^2 + \tau_t^2)I_n \right)
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
\Delta \sim N_{|S|}(\beta_0^{(S)} + \beta_1^{(S)}, \Sigma) \\
$$

where

$$
\Sigma = \text{diag} \left( \max\{ \hat{\sigma}_1^2, \tilde{\sigma}_1^2 \} + \hat{\tau}_1^2, ..., \max\{ \hat{\sigma}_{|S|}^2, \tilde{\sigma}_{|S|}^2 \} + \hat{\tau}_{|S|}^2 \right) \\
$$

With the single observation of this distribution, $\delta$, the James-Stein estimate for the unknown effect coefficient, $\beta_1^{(S)}$, is then

$$
\hat{\beta}_t^{(S)} = \left( 1 - \frac{c}{\delta \Sigma^{-1} \delta^T} \right) \delta - \hat{\beta}_0^{(S)}
$$

where $\hat{\beta}_0^{(S)}$ is the estimate obtained from the non-mutated samples for all transcripts $s \in S$.

It is simple to see that $\text{Tr}(\Sigma) = \sum_{s \in S} \max\{ \hat{\sigma}_s^2, \tilde{\sigma}_s^2 \} + \hat{\tau}_s^2$ and that $\lambda_L = \max_{s \in S} \left\{ \max\{ \hat{\sigma}_s^2, \tilde{\sigma}_s^2 \} + \hat{\tau}_s^2 \right\}$.
In the case of $\text{Tr}(\Sigma) < 2 \lambda_L$ or $|S| \le 2$, we resort to the OLS estimator, which is found through the standard `sleuth` process.

## Properties of the James-Stein estimator for differential expression

### The James-Stein estimator is biased towards 0

$$
\mathbb{E} \left[ \hat{\beta}_1^{(S)} \right] = \mathbb{E}[\delta] - c \mathbb{E} \left[ \frac{\delta}{\delta \Sigma^{-1} \delta^T} \right] - \mathbb{E} \left[ \hat{\beta}_0^{(S)} \right] \\
= \beta_0^{(S)} + \beta_1^{(S)} - c \mathbb{E} \left[ \frac{\delta}{\delta \Sigma^{-1} \delta^T} \right] - \beta_0^{(S)} \\
= \beta_1^{(S)} - c \mathbb{E} \left[ \frac{\delta}{\delta \Sigma^{-1} \delta^T} \right] \\
\overset{?}{=} \left( 1 - \frac{c \mu \mu^T}{\mu \Sigma^{-1} \mu^T} \right) \beta_1^{(S)} \\
$$

with the last step occurring by symmetry of the normal distribution about $\mu$ (and is something I still need to prove).
Note that this is true for all transcripts.

By Bock \Cref{https://projecteuclid.org/download/pdf_1/euclid.aos/1176343009}, this estimator is the minimax estimator under the mean square error.
The mean square error, $\mathbb{E} \left[ \Vert \hat{\mu} - \mu \Vert ^2\right] = \sum_{i=1}^{|S|} \mathbb{E}\left[ (\hat{\mu}_i - \mu_i)^2 \right] = \sum_{i=1}^{|S|} \text{Var}\left[ \hat{\mu}_i \right] = \text{Var}\left[ \bar{1}^T\hat{\mu} \right]$.