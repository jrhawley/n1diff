# Moments of the James-Stein estimator for differential expression

## Expected value

Due to the non-linear nature of the James-Stein estimator, we use a Taylor expansion around $\Beta_1$ to approximate the expectation.
We have

$$
\hat{\Beta}_1^{(JS)} = \left( 1 - \frac{c}{\left( \hat{\Beta}_1^{(OLS)} \right)^T \Sigma^{-1} \hat{\Beta}_1^{(OLS)}} \right) \hat{\Beta}_1^{(OLS)}
$$

where

$$
\hat{\Beta}_1^{(OLS)} \sim N_{|S|}(\Beta_1, \Sigma) \\

\Sigma_{s,s} = \left( \frac{n_{nonmut} + 1}{n_{nonmut}} \right) (\sigma_t^2 + \tau_t^2) \\
\Sigma_{s,t} = 0 \forall t \ne s
$$

Let

$$
u = \Sigma^{-1/2}\hat{\Beta}_1^{(OLS)}
$$

Then

$$
\mathbb{E}\left[ \hat{\Beta}_1^{(JS)} \right] = \mathbb{E} \left[ \hat{\Beta}_1^{(OLS)} \right] - c\Sigma^{1/2}\mathbb{E} \left[ \frac{u}{\Vert u \Vert ^2} \right] \\
 = \Beta_1 - c\Sigma^{1/2}\mathbb{E} \left[ \frac{u}{\Vert u \Vert ^2} \right] \\
$$

Expanding $\frac{u}{\Vert u \Vert ^2}$ around $a = \Sigma^{-1/2}\Beta_1$, we get

$$
\mathbb{E}\left[ \hat{\Beta}_1^{(JS)} \right]  = \Beta_1 - c\Sigma^{1/2}\mathbb{E} \left[ \frac{a}{\Vert a \Vert ^2} + \left( \frac{1}{\Vert a \Vert ^2}  - \frac{2}{\Vert a \Vert ^4} aa^T\right)(u - a) + \mathcal{O}(\Vert u - a \Vert ^2)\right] \\
 = \left(1 - \frac{c}{\Beta_1^T \Sigma^{-1} \Beta_1} \right) \Beta_1 + \mathcal{O}(\Vert u - a \Vert ^2) \\
$$

As long as

1. the number of transcripts being considered is not large, and the curse of dimensionality does not come into play, and
2. the true coefficient of variation is not large, so $\Vert u - a \Vert ^2$ is negligible compared to $\Beta_1$

the Taylor approximation should be close to $\left(1 - \frac{c}{\Beta_1^T \Sigma^{-1} \Beta_1} \right) \Beta_1$, which is an estimate of $\Beta_1$ biased towards 0.

## Variance

This estimator is the minimax estimator under the mean square error, $\mathbb{E} \left[ \Vert \hat{\Beta}_1^{(JS)} - \Beta_1 \Vert ^2\right] = \sum_{s \in S} \mathbb{E}\left[ \left( \hat{\Beta}_{1,s}^{(JS)} - \Beta_{1,s} \right)^2 \right] = \sum_{s \in S} \text{Var}\left[ \hat{\Beta}_{1,s}^{(JS)} \right]$.
While $\mathbb{E} \left[ \Vert \hat{\Beta}_1^{(JS)} - \Beta_1 \Vert ^2\right] \le \mathbb{E} \left[ \Vert \hat{\Beta}_1^{(OLS)} - \Beta_1 \Vert ^2\right]$, this does not imply that $\text{Var}\left[ \hat{\Beta}_{1,s}^{(JS)} \right] \le \text{Var}\left[ \hat{\Beta}_{1,s}^{(OLS)} \right] \forall s \in S$.
Some transcripts may have larger variances than the OLS estimator, but all transcripts in aggregate will have a smaller mean square error.
This is still desirable if the goal is to find if there is an effect on any transcripts in the set $S$, instead of a particular one within the set.

To calculate the variance for each individual transcript, we can take the same approach as above.

$$
\text{Var} \left [ \hat{\Beta}_1^{(JS)} \right] \approx \mathbb{E}\left[ \hat{\Beta}_1^{(JS)} \left( \hat{\Beta}_1^{(JS)} \right)^T \right] - \left(1 - \frac{c}{\Beta_1^T \Sigma^{-1} \Beta_1} \right)^2 \Beta_1 \Beta_1^T \\
= \Sigma^{1/2}\mathbb{E}\left[ u u^T - \frac{2c}{u^Tu}uu^T + \left( \frac{c}{u^Tu} \right)^2 uu^T \right]\Sigma^{1/2} - \left(1 - \frac{c}{\Beta_1^T \Sigma^{-1} \Beta_1} \right)^2 \Beta_1 \Beta_1^T \\
$$

where again $u = \Sigma^{-1/2} \hat{\Beta}_1^{(OLS)}$.
Expanding about $a = \Sigma^{-1/2} \Beta_1$,

$$
\text{Var} \left [ \hat{\Beta}_1^{(JS)} \right] = \left(1 - \frac{c}{\Beta_1^T \Sigma^{-1} \Beta_1} \right)^2 \Sigma + \left[ \left(1 - \frac{c}{\Beta_1^T \Sigma^{-1} \Beta_1} \right)^2 + \frac{2c}{\left( \Beta_1^T \Sigma ^{-1} \Beta_1 \right)^2 } - 1 \right] \Beta_1 \Beta_1^T + \mathcal{O}(\Vert u - a \Vert ^4) \\
$$

We see the shrinkage factor, squared, multiplying $\Sigma$, providing the possibility of a smaller estimator variance.
Unfortunately, $\Beta_1$ is unknown, so this estimator variance is a function of both the mean and variance of the transcripts under consideration.
This is in contrast to the OLS estimator, which is solely a function of the variance and experimental design.
