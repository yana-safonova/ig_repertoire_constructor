require("gtools")
require("parallel")
require("lattice")

nucl = 4

gen.samp <- function(N, real.lambda, lower.dataset.size=100, upper.datasets.size=1000) {
  stopifnot(length(real.lambda) == nucl)
  ns <- sample(lower.dataset.size:upper.datasets.size, size = N, replace = TRUE) # sizes of the datasets
  multinomial.param <- rdirichlet(N, real.lambda) # parameters for multinomial
  
  multinomial.samp <- t(mapply(function(n, prob) rmultinom(1, size=n, prob=prob), ns, split(multinomial.param, row(multinomial.param))))
  multinomial.samp
}

minuslog <- function(lambda, samp) {
  lambda.matrix <- matrix(lambda, byrow = TRUE, nrow = nrow(samp), ncol = ncol(samp))
  first <- samp + lambda.matrix
  lnP1 <- rowSums(lgamma(first))
  lnG1 <- lgamma(rowSums(first))
  
  lnP2 <- rowSums(lgamma(lambda.matrix))
  lnG2 <- lgamma(rowSums(lambda.matrix))
  
  -sum(lnP1 - lnG1 - lnP2 + lnG2)
}
parnames(minuslog) <- as.character(1:nucl)

minuslog.grad <- function(lambda, samp) {
  lambda.matrix <- matrix(lambda, byrow = TRUE, nrow = nrow(samp), ncol = ncol(samp))
  - colSums(
    digamma(samp + lambda.matrix) - digamma(rowSums(samp + lambda.matrix)) -
      digamma(lambda.matrix) + digamma(rowSums(lambda.matrix))
  )
}

### Check grad
grad.times.check = 1e1
replicate(grad.times.check, {
  x = runif(n = 4, min = 0, max = 1000)
  stopifnot(all.equal(grad(minuslog, x = x, samp = samp), 
                      minuslog.grad(lambda = x, samp = samp), tolerance = 1e-5))
})

samp <- gen.samp(N=10, real.lambda=1:nucl + 100)
start.value = list(lambda=1:nucl)
# mle2(minuslog, start=start.value, lower=list(lambda=rep(0.1, nucl)), method="L-BFGS-B", data = list(samp=samp), vecpar = TRUE, parnames = 1:10)
# mle2(minuslog, start=start.value, data=list(samp=samp), vecpar = TRUE, parnames=as.character(1:nucl))
optim(par=rep(10, nucl), gr = minuslog.grad, fn=minuslog, samp=samp, lower = rep(0.1, nucl), method = "L-BFGS-B")

mse <- function(samp, a) {
  sqrt(mean((samp - a)^2)) 
}

M <- 1e2

mse.calc <- function(M, N, real.lambda, start.value) {
  results <- mclapply(1:M, 
  #results <- replicate(M,
                       function(ind, ...) {
                       # {
                         samp <- gen.samp(N = N, real.lambda = real.lambda)
                         optim.res <- optim(par = start.value, gr = minuslog.grad, fn=minuslog, samp=samp, lower = rep(0.1, nucl), method = "L-BFGS-B")
                         optim.res$par
                       }, mc.cores=4)
  results <- simplify2array(results)
  mse(results, real.lambda)
}

Ns <- 10 * (1:20)
Ns <- 2^(7:12)
mses <- sapply(Ns, mse.calc, M = 1e2, real.lambda = 1:nucl, start.value = rep(10, nucl))
mses
plot(log(mses), log(Ns))

xyplot(log(mses) ~ log(Ns), xlab = "Sample size", ylab = "MSE",
       panel = function(x, y, ...) {
          panel.xyplot(x, y, col="black", cex=2, pch = 21)
          panel.lmline(x, y, col = "salmon")
       }
)
lm(log(mses) ~ log(Ns))
