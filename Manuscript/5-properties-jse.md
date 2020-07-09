# Moments of the James-Stein estimator for differential expression

## Expected value

?

## Variance

This estimator is the minimax estimator under the mean square error, $\mathbb{E} \left[ \Vert \hat{\Beta}_1^{(JS)} - \Beta_1 \Vert ^2\right] = \sum_{s \in S} \mathbb{E}\left[ \left( \hat{\Beta}_{1,s}^{(JS)} - \Beta_{1,s} \right)^2 \right] = \sum_{s \in S} \text{Var}\left[ \hat{\Beta}_{1,s}^{(JS)} \right]$.
While $\mathbb{E} \left[ \Vert \hat{\Beta}_1^{(JS)} - \Beta_1 \Vert ^2\right] \le \mathbb{E} \left[ \Vert \hat{\Beta}_1^{(OLS)} - \Beta_1 \Vert ^2\right]$, this does not imply that $\text{Var}\left[ \hat{\Beta}_{1,s}^{(JS)} \right] \le \text{Var}\left[ \hat{\Beta}_{1,s}^{(OLS)} \right] \forall s \in S$.
Some transcripts may have larger variances than the OLS estimator, but all transcripts in aggregate will have a smaller mean square error.
This is still desirable if the goal is to find if there is an effect on any transcripts in the set $S$, instead of a particular one within the set.
