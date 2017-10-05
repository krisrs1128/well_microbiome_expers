
################################################################################
# Simulate data according to Daniela Witten's PMD paper
################################################################################

###############################################################################
## Libraries and setup
###############################################################################

library("reshape2")
library("tidyverse")
library("PMA")
set.seed(04032016)

opts <- list(
  "n" = 504, # want divisible by 8
  "p" = 20,
  "k" = 2,
  "l" = 2,
  "sigma0" = 5,
  "sigma" = 1
)

###############################################################################
## Simulate data
###############################################################################
S <- qr.Q(qr(matnorm(opts$p, opts$k, opts$sigma0))) # same source between the two matrices
W <- replicate(opts$l, matrix(0, opts$n, opts$k), simplify = F)
W[[1]][, 1] <- c(
  rep(10, opts$n / 8),
  rep(-10, opts$n / 8),
  rep(10, opts$n / 8),
  rep(-10, opts$n / 8),
  rep(0, opts$n / 2)
)
W[[1]][, 2] <- c(
  rep(10, opts$n / 4),
  rep(-10, opts$n / 4),
  rep(0, opts$n / 2)
)
W[[2]][, 1] <- c(
  rep(0, opts$n / 4),
  rep(10, opts$n / 2),
  rep(-10, opts$n / 4)
)
W[[2]][, 2] <- c(
  rep(-10, opts$n / 4),
  rep(10, opts$n / 4),
  rep(-10, opts$n / 4),
  rep(10, opts$n / 4)
)

###############################################################################
## Plot the simulated data
###############################################################################
mW <- melt(W)
colnames(mW) <- c("i", "k", "w", "table")
mW$k <- paste0("Latent Dimension ", mW$k)
mW$table[mW$table == 1] <- "X"
mW$table[mW$table == 2] <- "Y"

ggplot(mW) +
  geom_point(aes(x = i, y = w)) +
  xlab("sample index") + 
  facet_grid(table ~ k) +
  ggtitle(paste0("True latent weights W"))

## another simulation
X <- common_source_model(W, S, opts)

colnames(X[[1]]) <- paste0("X", seq_len(ncol(X[[1]])))
colnames(X[[2]]) <- paste0("Y", seq_len(ncol(X[[2]])))
x_ix <- sample(seq_len(ncol(X[[1]])), 4)
y_ix <- sample(seq_len(ncol(X[[1]])), 4)
pairs(X[[1]][, x_ix], asp = 1, main = "Four columns of X")
pairs(X[[2]][, y_ix], asp = 1, main = "Four columns of Y")

pairs(cbind(X[[1]][, x_ix[1:2]], X[[2]][, y_ix[1:2]]), asp = 1,
      main = "Two columns of X vs. Two columns of Y")

X <- lapply(X, scale)
pmd_res <- MultiCCA(lapply(X, function(x) t(x)), penalty = 10, ncomponents = 3)

mW_hat <- melt(pmd_res$ws)
colnames(mW_hat) <- c("i", "k", "w", "table")
mW_hat <- mW_hat %>% filter(k < 4)
mW_hat$k <- paste0("Recovered Dimension ", mW_hat$k)
mW_hat$table[mW_hat$table == 1] <- "X"
mW_hat$table[mW_hat$table == 2] <- "Y"

ggplot(mW_hat) +
  geom_point(aes(x = i, y = w), alpha = 0.6, size = 1) +
  facet_grid(table ~ k) +
  ggtitle(expression(paste("Recovered Weights ", hat(W), " [PMD]")))

pmd_res <- MultiCCA(lapply(X, function(x) t(x)), penalty = 1, type = "ordered")
