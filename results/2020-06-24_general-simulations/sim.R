library("mnormt")
library("psych")
library("data.table")
library("ggplot2")


first_prod <- function(x, s_inv) {
    t_x <- t(x)
    const <- 1 / (t_x %*% s_inv %*% x)[1, 1]
    return(const * (x %*% t_x))
}

second_prod <- function(x, s_inv) {
    t_x <- t(x)
    const <- 1 / (t_x %*% s_inv %*% x)[1, 1]
    return(const^2 * (x %*% t_x))
}

jse <- function(x, s, s_inv, b0) {
    const = tr(s) / max(diag(s)) - 2
    prod_const <- (t(x) %*% s_inv %*% x)[1, 1]
    return((1 - const / prod_const) * x - b0)
}

set.seed(42)

beta_0 <- rnorm(100, runif(1, 0, 100), 1)^2
beta_1 <- rnorm(100, runif(1, 0, 100), 1)
folds <- rnorm(100, 0, 2)
sigma <- diag(folds^2)
sigma_inv <- diag(folds^(-2))

# e_1l <- lapply(1:nrow(delta), function(i) first_prod(delta[i, ], sigma_inv))
# e_1 <- Reduce("+", e_1l) / length(e_1l)
# e_2l <- lapply(1:nrow(delta), function(i) second_prod(delta[i, ], sigma_inv))
# e_2 <- Reduce("+", e_2l) / length(e_2l)
# max(abs(e_1))
# max(abs(e_2))

delta <- rmnorm(10000, beta_0 + beta_1, sigma)
b1_hat <- lapply(1:nrow(delta), function(i) {jse(delta[i, ], sigma, sigma_inv, beta_0)})

var_b1_hat_prod <- lapply(
    b1_hat,
    function(b1h) {
        v <- b1h - beta_1
        return(v %*% t(v))
    }
)

var_b1_hat <- Reduce("+", var_b1_hat_prod) / length(var_b1_hat_prod)

n_nonmut <- 11
ols_var <- (n_nonmut + 1)/n_nonmut * diag(sigma)

dt <- data.table(
    Index = 1:100,
    JSE_Var = diag(var_b1_hat),
    OLS_Var = ols_var
)

gg <- (
    ggplot(data = dt)
    + geom_point(aes(x = Index, y = OLS_Var - JSE_Var), colour = "#ba839d")
    #+ geom_point(aes(x = Index, y = JSE_Var), colour = "#924023")
)
ggsave("sim-jse.png", width = 10, height = 6, units = "cm")
