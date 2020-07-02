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
