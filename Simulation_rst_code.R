rm(list=ls())

library(devtools)
# load_all()
devtools::install_github("ruihu51/CausalSim")



require(CausalSim)

# unit test for htem.estimator()
null.sims <- FALSE

n <- 1000
set.seed(39338064)
W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
if(null.sims) {
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
  psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
  theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
} else {
  Y <- rbinom(n, size = 1, prob = mu0(A, W))
  psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
  theta0 <- var((mu0(1,W) - mu0(0,W)))
}

ret1 <- htem.estimator(A, W, Y, control = list()) # default control setting-point estimate
ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE)) # estimate with wald-type CI
ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE, conf.int.type = 'boot')) # estimate with Bootstrap CI

# default control setting-point estimate with only hybrid outputs
# $ret
# type        est        se
# 1   psi.est 0.13016493 0.7394138
# 2 theta.est 0.01530881 0.2236576
ret1 <- htem.estimator(A, W, Y, control = list())

# default control setting-point estimate with additional optional outputs
# $ret
# type        est        se
# 1   psi.est 0.12898826 0.7441527
# 2 theta.est 0.01345503 0.2204151
#
# $optional.ret
# type       est        se
# 1  psi.plug.in.est 0.1317073 0.7441527
# 2 psi.one.step.est 0.1289883 0.7441527
ret1 <- htem.estimator(A, W, Y, control = list(est.type = list('psi.est'='all', 'theta.est'='hybrid')))

# estimate with wald-type CI
# $ret
# type       est       se          ll        ul
# 1   psi.est 0.1310583 0.737582 0.085343285 0.1767732
# 2 theta.est 0.0152913 0.221020 0.001592587 0.0289900
#
# $optional.ret
# type       est       se
# 1    psi.plug.in.est 0.1338398 0.737582
# 2   psi.one.step.est 0.1310583 0.737582
# 3  theta.plug.in.est 0.0197228 0.221020
# 4 theta.one.step.est 0.0152913 0.221020
ret1 <- htem.estimator(A, W, Y, control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
                                               conf.int.type = 'Wald',
                                               conf.int = TRUE))

# estimate with Bootstrap CI
# $ret
# type        est        se          ll        ul
# 1   psi.est 0.13071134 0.7381247 0.082714492 0.1739515
# 2 theta.est 0.01493656 0.2201738 0.001910335 0.0288994
#
# $optional.ret
# type        est        se
# 1    psi.plug.in.est 0.13358361 0.7381247
# 2   psi.one.step.est 0.13071134 0.7381247
# 3  theta.plug.in.est 0.01962606 0.2201738
# 4 theta.one.step.est 0.01493656 0.2201738
ret1 <- htem.estimator(A, W, Y, control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
                                               conf.int = TRUE,
                                               conf.int.type = 'boot',
                                               n.boot = 500))

control = list(est.type = list('psi.est'='all', 'theta.est'='hybrid'))
control <- hte.measure.NullTest.control(control)

# simulation function
control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
ests.sim.1.1.sl <-  ests.sim(200, 1:500, control, out.glm=TRUE)
control = list(conf.int = TRUE,
               conf.int.type = 'boot')
ests.sim.1.2.sl <-  ests.sim(200, 1:500, control, out.glm=TRUE)

n<-1000
do.call(generate.data, list(n=n))
generate.data(n)




# simulation function
ests.sim <-  function(n_range, j_range, control, null.sims=FALSE, out.glm=TRUE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      # print(j)
      if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      W <- matrix(runif(n*3, 0, 1), ncol=3)
      A <- rbinom(n, size = 1, prob = pi0(W))
      if(null.sims) {
        Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
        psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
        theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
      } else {
        Y <- rbinom(n, size = 1, prob = mu0(A, W))
        psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
        theta0 <- var((mu0(1,W) - mu0(0,W)))
      }

      prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
      pi.hat <- prop.reg$fitted.values
      AW <- cbind(A, data.frame(W))
      mu.reg <- glm(Y ~ ., data=AW, family='binomial')
      mu.hat <- mu.reg$fitted.values
      mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                                 data.frame(W)), type = 'response')
      mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                                 data.frame(W)), type = 'response')
      mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
      if (out.glm){
        control$pi.hat = pi.hat
        control$mu.hats = mu.hats
      }
      #print(control)
      ret <- htem.estimator(A, W, Y, control = control)

      ret$pi0 <- rep(mean(pi0(W)), 2)
      ret$pi.hat <- rep(mean(pi.hat), 2)
      ret$mu0 <- rep(mean(mu0(A, W)), 2)
      ret$mu.hat <- rep(mean(A*mu.hats$mu1 + (1-A)*mu.hats$mu0), 2)
      ret$n <- rep(n, 2)
      ret$j <- rep(j, 2)
      ret$seed <- rep(seed, 2)

      ret$psi0 <- rep(psi0, 2)
      ret$theta0 <- rep(theta0, 2)

      return(ret)
    })
  })

  return(ests)
}

# numerical analysis
## case 1
## Wald-type
ests.sim.1 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE), null.sims=FALSE, out.glm=TRUE)
ests.sim.1.new <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE), null.sims=FALSE, out.glm=TRUE)
ests.sim.1.1 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE), null.sims=FALSE, out.glm=TRUE)
ests.sim.1 <- ests.sim.1.new
save(ests.sim.1, file = "ests.sim.1.RData")
save(ests.sim.1.new, file = "ests.sim.1.new.RData")
save(ests.sim.1.1, file = "ests.sim.1.1.RData")
## Bootstrap
ests.sim.2 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2, file = "ests.sim.2.RData")

ests.sim.2.new.h1 <- ests.sim(200*c(1:5), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2.new.h1, file = "ests.sim.2.new.h1.RData")
ests.sim.2.new.h2 <- ests.sim(200*c(6:7), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2.new.h2, file = "ests.sim.2.new.h2.RData")
ests.sim.2.new.h3 <- ests.sim(200*c(8:10), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2.new.h3, file = "ests.sim.2.new.h3.RData")
ests.sim.2.new <- rbind(rbind(ests.sim.2.new.h1, ests.sim.2.new.h2), ests.sim.2.new.h3)
save(ests.sim.2.new, file = "ests.sim.2.new.RData")

ests.sim.2.1.h1 <- ests.sim(200*c(1:5), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2.1.h1, file = "ests.sim.2.1.h1.RData")
ests.sim.2.1.h2 <- ests.sim(200*c(6:8), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2.1.h2, file = "ests.sim.2.1.h2.RData")
ests.sim.2.1.h3 <- ests.sim(200*c(9:10), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2.1.h3, file = "ests.sim.2.1.h3.RData")
ests.sim.2.1 <- rbind(rbind(ests.sim.2.1.h1, ests.sim.2.1.h2), ests.sim.2.1.h3)
save(ests.sim.2.1, file = "ests.sim.2.1.RData")

## plot
ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes((pi0), color=type)) +
  facet_wrap(~n, scales='free')
ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes((mu0), color=type)) +
  facet_wrap(~n, scales='free')

ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes(sqrt(n) * (est - psi0), color=type)) +
  facet_wrap(~n, scales='free')
ggplot(subset(ests.sim.1, (type %in% 'theta.est'))) +
  geom_density(aes(sqrt(n) * (est - theta0), color=type)) +
  facet_wrap(~n, scales='free')

ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes((pi0 - pi.hat), color=type)) +
  facet_wrap(~n, scales='free')
ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes((mu0 - mu.hat), color=type)) +
  facet_wrap(~n, scales='free')

## summary and plot
psi.summaries.1 <- ddply(subset(ests.sim.1, (type %in% 'psi.est')), .(n, type), summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         bias = mean(est - psi0, na.rm=TRUE),
                         var = var(est, na.rm=TRUE),
                         mse = mean((est - psi0)^2, na.rm=TRUE))
library(ggplot2)
p1 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, coverage, color=type)) +
  geom_hline(yintercept=.95)

p2 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, sqrt(n) * bias, color=type))

p3 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, n * mse, color=type))

p4 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, n * var, color=type)) +
  coord_cartesian(ylim=c(0,1))

library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)

ggplot(psi.summaries.1) +
  geom_line(aes(n, bias, color=type))

ests.sim.1$psi.coverage <- (ests.sim.1$ll <= ests.sim.1$psi0) & (ests.sim.1$psi0 <= ests.sim.1$ul)
ests.sim.1$theta.coverage <- (ests.sim.1$ll <= ests.sim.1$theta0) & (ests.sim.1$theta0 <= ests.sim.1$ul)
ggplot(subset(ests.sim.1, (type %in% 'psi.est') & (n %in% c(200, 600, 1000, 2000)) &
                (psi.coverage ==FALSE))) +
  geom_linerange(mapping = aes(x=j, ymin=ll, ymax=ul)) +
  geom_point(mapping=aes(x=j, y=psi0), size=1, col="red") +
  facet_wrap(~n, scales='free')


# case 2
null.sims <- TRUE

n <- 1000
set.seed(39338064)
W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
if(null.sims) {
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
  psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
  theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
} else {
  Y <- rbinom(n, size = 1, prob = mu0(A, W))
  psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
  theta0 <- var((mu0(1,W) - mu0(0,W)))
}

ret1 <- htem.estimator(A, W, Y, control = list())
ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE)) # estimate with wald-type CI
ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE, conf.int.type = 'boot')) # estimate with Bootstrap CI

ests.sim.3.new <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE), null.sims=TRUE, out.glm=TRUE)
save(ests.sim.3.new, file = "ests.sim.3.new.RData")

ggplot(subset(ests.sim.3, (type %in% 'psi.est'))) +
  geom_density(aes(sqrt(n) * (est - psi0), color=type)) +
  facet_wrap(~n, scales='free')
ggplot(subset(ests.sim.3, (type %in% 'theta.est'))) +
  geom_density(aes(sqrt(n) * (est - theta0), color=type)) +
  facet_wrap(~n, scales='free')

psi.summaries.3 <- ddply(subset(ests.sim.3, (type %in% 'psi.est')), .(n, type), summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         bias = mean(est - psi0, na.rm=TRUE),
                         var = var(est, na.rm=TRUE),
                         mse = mean((est - psi0)^2, na.rm=TRUE))

library(ggplot2)
p1 <- ggplot(psi.summaries.3) +
  geom_line(aes(n, coverage, color=type)) +
  geom_hline(yintercept=.95)

p2 <- ggplot(psi.summaries.3) +
  geom_line(aes(n, sqrt(n) * bias, color=type))

p3 <- ggplot(psi.summaries.3) +
  geom_line(aes(n, n * mse, color=type))

p4 <- ggplot(psi.summaries.3) +
  geom_line(aes(n, n * var, color=type)) +
  coord_cartesian(ylim=c(0,1))

library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)


# testing
cov.var <- FALSE

null.sims <- FALSE

n <- 500
set.seed(3933)
W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
if(null.sims) {
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
  psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
  theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
} else {
  Y <- rbinom(n, size = 1, prob = mu0(A, W))
  psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
  theta0 <- var((mu0(1,W) - mu0(0,W)))
}

control = list(verbose=TRUE)
hteNullTest(Y, A, W, control = list(verbose=TRUE), out.glm=FALSE, cov.var=FALSE)


test.sim.1.sl <- testing.sim(500, 1:10, control = list(), out.glm = TRUE)
test.sim.1.sl

control = list()
control <- hte.measure.NullTest.control(control)
n = length(A)

out.glm =TRUE
if (out.glm){
  prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
  pi.hat <- prop.reg$fitted.values
  AW <- cbind(A, data.frame(W))
  mu.reg <- glm(Y ~ ., data=AW, family='binomial')
  mu.hat <- mu.reg$fitted.values
  mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                             data.frame(W)), type = 'response')
  mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                             data.frame(W)), type = 'response')
  mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
  control$pi.hat = pi.hat
  control$mu.hats = mu.hats
}

cat("Estimating...")
tm0 <- proc.time()
# estimated P(A = 1 | W = w)
if(is.null(control$pi.hat)){
  prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                           newX = data.frame(W),
                           SL.library = control$pi.SL.library,
                           family = binomial(),
                           obsWeights=rep(1,n),
                           id=1:n)
  control$pi.hat <- prop.reg$SL.predict
}

# estimated mu
if(is.null(control$mu.hats)){
  AW <- cbind(A, data.frame(W))
  if(length(setdiff(Y, c(0,1))) == 0) {
    family = 'binomial'
  } else {
    family = 'gaussian'
  }
  mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                         newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                         SL.library = control$mu.SL.library,
                         family = family,
                         obsWeights=rep(1,n),
                         id=1:n)
  control$mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
}

control$mu.hat <- A * control$mu.hats$mu1 + (1-A) * control$mu.hats$mu0

# estimated tau
tau.hat <- control$mu.hats$mu1 - control$mu.hats$mu0
Z.hat <- (2*A - 1) / (A * control$pi.hat + (1-A) * (1-control$pi.hat))

# estimated theta
gamma.hat <- mean(tau.hat)

w.ecdf <- function(w) {
  vec.leq <- function(x,y) prod(x <= y)
  return(mean(apply(W, 1, vec.leq, y = w)))
}
u.vals <- apply(W, 1, w.ecdf)

# primitive function
cat("Computing Gamma and Omega...\n")
tm1 <- proc.time()
# n.new * 1 vector
w.vals <- W
Gamma.w.vals <- apply(w.vals, 1, function(w0)
  mean(apply(W, 1, function(x) prod(x <= w0))*tau.hat))
Omega.w.vals <- Gamma.w.vals - gamma.hat * u.vals

# nonparametric EIF
# n * n.new matrix
vec.eq <- function(x,y) prod(x == y)
eif.Gamma <- apply(w.vals, 1, function(w0) {
  (apply(W, 1, function(x) prod(x <= w0))) * (Z.hat * (Y - control$mu.hat) + tau.hat) -
    Gamma.w.vals[which(as.logical(apply(W, 1, vec.eq, y = w0)))]
})
eif.Omega <- apply(w.vals, 1, function(w0) {
  (apply(W, 1, function(x) prod(x <= w0)) - w.ecdf(w0)) * (Z.hat * (Y - control$mu.hat) + tau.hat - gamma.hat) -
    Omega.w.vals[which(as.logical(apply(W, 1, vec.eq, y = w0)))]
})

# one-step estimators
# n.new * 1 vector
Gamma.os.est <- colMeans(eif.Gamma) + Gamma.w.vals
Omega.os.est <- colMeans(eif.Omega) + Omega.w.vals

# testing procedure
cat("Computing statistics...\n")
tm2 <- proc.time()
# test statistics
Gamma.stat <- n^(1/2)*max(abs(Gamma.os.est))
Omega.stat <- n^(1/2)*max(abs(Omega.os.est))

# covariance matrices
n.new <- dim(w.vals)[1]
if (cov.var){
  Gamma.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
    mean(eif.Gamma[,s] * eif.Gamma[,t])
  }))
  Omega.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
    mean(eif.Omega[,s] * eif.Omega[,t])
  }))

  # quantiles
  tm3 <- proc.time()
  Gamma.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Gamma.cov.var)
  Gamma.epsilon.stats <- apply(Gamma.epsilon, 1, function(x) {max(abs(x))})
  Omega.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Omega.cov.var)
  Omega.epsilon.stats <- apply(Omega.epsilon, 1, function(x) {max(abs(x))})
}else{
  tm3 <- proc.time()
  Gamma.epsilon.stats <- replicate(control$n.boot, max(abs(t(eif.Gamma)%*%rnorm(n.new, 0, 1)/sqrt(n))))
  Omega.epsilon.stats <- replicate(control$n.boot, max(abs(t(eif.Omega)%*%rnorm(n.new, 0, 1)/sqrt(n))))
}

Gamma.pvalue <- mean(Gamma.epsilon.stats > Gamma.stat)
Gamma.quantile <- unname(quantile(Gamma.epsilon.stats, control$conf.level))
Omega.pvalue <- mean(Omega.epsilon.stats > Omega.stat)
Omega.quantile <- unname(quantile(Omega.epsilon.stats, control$conf.level))

ret <- data.frame(type = 'Gamma.stat', stat = Gamma.stat, pvalue = Gamma.pvalue,
                  quantile = Gamma.quantile)
ret <- rbind(ret,
             data.frame(type = 'Omega.stat', stat = Omega.stat, pvalue = Omega.pvalue,
                        quantile = Omega.quantile))

ret
tm4 <- proc.time()

tm4-tm0
tm4-tm3
tm3-tm2
tm2-tm1
tm1-tm0

## compare
control$n.boot <- 1000
Gamma.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Gamma.cov.var)
Gamma.epsilon.stats.0 <- apply(Gamma.epsilon, 1, function(x) {max(abs(x))})
Gamma.epsilon.stats
Gamma.epsilon.stats <- replicate(control$n.boot, max(abs(t(eif.Gamma)%*%rnorm(n.new, 0, 1)/sqrt(n))))
lines(density(Gamma.epsilon.stats.0),col='red')
plot(density(Gamma.epsilon.stats))
##
## covariance computing test
eif.Gamma <- apply(w.vals, 1, function(w0) {
  (apply(W, 1, function(x) prod(x <= w0))) * (Z.hat * (Y - control$mu.hat) + tau.hat) -
    Gamma.w.vals[which(as.logical(apply(W, 1, vec.eq, y = w0)))]
})

tm0 <- proc.time()
Gamma.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
  mean(eif.Gamma[,s] * eif.Gamma[,t])
}))
Gamma.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Gamma.cov.var)
Gamma.epsilon.stats <- apply(Gamma.epsilon, 1, function(x) {max(abs(x))})
Gamma.pvalue <- mean(Gamma.epsilon.stats > Gamma.stat)
tm1 <- proc.time()

Gamma.epsilon.stats.1 <- replicate(control$n.boot, max(abs(t(eif.Gamma)%*%rnorm(n.new, 0, 1)/sqrt(n))))
Gamma.pvalue.1 <- mean(Gamma.epsilon.stats.1 > Gamma.stat)
tm2 <- proc.time()

tm1-tm0
tm2-tm1

## testing simulation

null.sims <- FALSE

n <- 400
set.seed(3933)
W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
if(null.sims) {
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
  psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
  theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
} else {
  Y <- rbinom(n, size = 1, prob = mu0(A, W))
  psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
  theta0 <- var((mu0(1,W) - mu0(0,W)))
}
test.rst1 <- hteNullTest(Y, A, W, control = list(), out.glm=TRUE)

null.sims <- TRUE

n <- 800
set.seed(3933)
W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
if(null.sims) {
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
  psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
  theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
} else {
  Y <- rbinom(n, size = 1, prob = mu0(A, W))
  psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
  theta0 <- var((mu0(1,W) - mu0(0,W)))
}
test.rst2 <- hteNullTest(Y, A, W)


testing.sim.1 <- testing.sim(200*c(1:1), 1:500, control = list(), null.sims=FALSE, out.glm=TRUE)
testing.rst.1.h3 <- testing.sim(200*c(3:4), 1:500, control = list(), null.sims=FALSE, out.glm=TRUE)
save(testing.rst.1.h3, file = "testing.rst.1.h3.RData")
testing.rst.1.h4 <- testing.sim(200*c(5:5), 1:500, control = list(), null.sims=FALSE, out.glm=TRUE)
save(testing.rst.1.h4, file = "testing.rst.1.h4.RData")

testing.rst.1.h2 <- ret
save(testing.rst.1.h2, file = "testing.rst.1.h2.RData")

load('../HTE_Simulation_Result/testing.rst.1.h5.in.RData')
testing.rst.1.h5.in$j

# separate loop safe
control = list()
null.sims=FALSE
out.glm=TRUE
ret <- data.frame()
for (n in 2000){
  for (j in 12:100){
    # if(j %% 100 == 0)
    cat(n, j, '\n')
    seed <- sample(1e3:1e8, 1)
    while (seed %in% testing.rst.1.h5.in$seed){
      cat('resample', seed, '\n')
      seed <- sample(1e1:1e3, 1)
      cat('new seed', seed, '\n')
    }
    set.seed(seed)

    W <- matrix(runif(n*3, 0, 1), ncol=3)
    A <- rbinom(n, size = 1, prob = pi0(W))
    if(null.sims) {
      Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
      psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
      theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
    } else {
      Y <- rbinom(n, size = 1, prob = mu0(A, W))
      psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
      theta0 <- var((mu0(1,W) - mu0(0,W)))
    }

    if (length(ret)>0){
      ret.tmp <- hteNullTest(Y, A, W, control = control, out.glm = out.glm, cov.var=FALSE)

      ret.tmp$n <- rep(n, 2)
      ret.tmp$j <- rep(j, 2)
      ret.tmp$seed <- rep(seed, 2)

      ret.tmp$psi0 <- rep(psi0, 2)
      ret.tmp$theta0 <- rep(theta0, 2)

      ret <- rbind(ret, ret.tmp)
    }
    else{
      ret <- hteNullTest(Y, A, W, control = control, out.glm = out.glm, cov.var=FALSE)

      ret$n <- rep(n, 2)
      ret$j <- rep(j, 2)
      ret$seed <- rep(seed, 2)

      ret$psi0 <- rep(psi0, 2)
      ret$theta0 <- rep(theta0, 2)
    }
  }
}
testing.rst.1.h5.in.2 <- ret
save(testing.rst.1.h5.in.2, file = "testing.rst.1.h5.in.2.RData")

testing.rst.1.h4.in <- ret
save(testing.rst.1.h4.in, file = "testing.rst.1.h4.in.RData")

testing.rst.1.h4.in2 <- ret
save(testing.rst.1.h4.in2, file = "testing.rst.1.h4.in2.RData")

testing.rst.2.h1.in <- ret
save(testing.rst.2.h1.in, file = "testing.rst.2.h1.in.RData")


# simulation for luedtke's code
mmd.test = function(R,S,D.R,D.S,sig.meth='eig',num.reps=1e4,return.cutoff=FALSE){
  n = length(R)

  D.R.mat1 = matrix(rep(D.R,n),nrow=n)
  D.R.mat2 = matrix(rep(D.R,each=n),nrow=n)

  D.S.mat1 = matrix(rep(D.S,n),nrow=n)
  D.S.mat2 = matrix(rep(D.S,each=n),nrow=n)

  R.mat1 = matrix(rep(R,n),nrow=n)
  R.mat2 = matrix(rep(R,each=n),nrow=n)

  S.mat1 = matrix(rep(S,n),nrow=n)
  S.mat2 = matrix(rep(S,each=n),nrow=n)

  EE = ((2*(R.mat1-R.mat2)*(D.R.mat2-D.R.mat1) + 1 - (4*(R.mat1-R.mat2)^2-2)*D.R.mat1*D.R.mat2)*exp(-(R.mat1-R.mat2)^2)
        - ((2*(S.mat1-R.mat2)*(D.R.mat2-D.S.mat1) + 1 - (4*(S.mat1-R.mat2)^2-2)*D.S.mat1*D.R.mat2)*exp(-(S.mat1-R.mat2)^2))
        - ((2*(R.mat1-S.mat2)*(D.S.mat2-D.R.mat1) + 1 - (4*(R.mat1-S.mat2)^2-2)*D.R.mat1*D.S.mat2)*exp(-(R.mat1-S.mat2)^2))
        + (2*(S.mat1-S.mat2)*(D.S.mat2-D.S.mat1) + 1 - (4*(S.mat1-S.mat2)^2-2)*D.S.mat1*D.S.mat2)*exp(-(S.mat1-S.mat2)^2))

  # EE = exp(-(R.mat1-R.mat2)^2) - 2*exp(-(S.mat1-R.mat2)^2) + exp(-(S.mat1-S.mat2)^2)

  if(sig.meth=='eig'){
    line.means = rowMeans(EE)
    EE.ctrd = EE - matrix(rep(line.means,n),nrow=n) - matrix(rep(line.means,each=n),nrow=n) + matrix(rep(mean(line.means),n^2),nrow=n)
    num.eigs = min(200,n)
    # tmp = eigen(EE.ctrd)$values/n
    require('rARPACK')
    tmp = eigs_sym(EE.ctrd,num.eigs,which='LA')$values/n
    num.pos.eigs = num.eigs # sum(tmp>0)
    draws=c(matrix(rnorm(num.reps*num.pos.eigs)^2-1,nrow=num.reps,ncol=num.pos.eigs)%*%cbind(tmp[1:num.pos.eigs]))
  }

  # U-statistic
  diag(EE) = 0
  est = (rbind(rep(1/(n-1),n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]
  # V-statistic
  # est = (rbind(rep(1/n,n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]

  if(sig.meth=='eig'){
    pval = mean(draws>n*est)
  } else if(sig.meth=='var'){
    pval = pchisq(est/(2*var(D.R)/n)+1,df=1,lower.tail=FALSE)
  }

  return(if(!return.cutoff){
    c(est,pval)
  }else{
    c(est,pval,
      if(sig.meth=='eig'){
        quantile(draws,0.95)
      }else{
        2*var(D.R)*(qchisq(0.95,df=1)-1)})})
}

# Tests if A is used in the regression function, i.e. if
# E[Y|A,W] = E[Y|W] almost surely, or equivalently (under positivity) if
# E[Y|A=1,W] = E[Y|A=0,W] almost surely
# Returns an estimate of Psi, and a p-value
# If est.g is false, then just uses predefined function g0
# NOTE: Theory all uses bounded Y, so if Y is not bounded then at least try to make sure most of its
# mass falls in the interval [-1,1] --- can do this by scaling by a constant
est.psi.prob = function(W,A,Y,W.train=NULL,A.train=NULL,Y.train=NULL,sig.meth='var',est.g=TRUE){
  require('SuperLearner')
  SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet')
  n=length(A)

  if(is.null(W.train) | is.null(A.train) | is.null(Y.train)){
    W.train = W
    A.train = A
    Y.train = Y
  }

  # Estimate outcome regressions
  Qbar.est = SuperLearner(Y=Y.train,X=data.frame(W=W.train,A=A.train),newX=data.frame(W=rbind(W,W),A=rep(c(0,1),each=n)),SL.library=SL.library)
  Qbar.est.0 = Qbar.est$SL.predict[,1][1:n]
  Qbar.est.1 = Qbar.est$SL.predict[,1][(n+1):(2*n)]

  if(est.g){
    # gg = SuperLearner(Y=A,
    #                   X=W,
    #                   SL.library=SL.library,
    #                   family='binomial')$SL.predict[,1]
    prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
    gg <- prop.reg$fitted.values
  } else {
    gg = g0(W)
  }

  # Plug-in estimate of blip
  R = Qbar.est.1 - Qbar.est.0
  S = rep(0,n)

  D.R = A/gg * (Y-Qbar.est.1) - (1-A)/(1-gg) * (Y-Qbar.est.0)
  D.S = rep(0,n)

  return(mmd.test(R,S,D.R,D.S,sig.meth=sig.meth))
}

load('../HTE_Simulation_Result/testing.rst.1.Rdata')
load('../HTE_Simulation_Result/testing.rst.2.Rdata')

null.sims <- TRUE
sim2019rst.2 <- ldply(c(200, 400, 800, 1000), function(n) {
  ldply(1:500, function(j) {
    # print(j)
    if(j %% 100 == 0) cat(n, j, '\n')
    seed <- testing.rst.2[(testing.rst.2$n==n) & (testing.rst.2$j==j), "seed"][1]
    set.seed(seed)

    W <- matrix(runif(n*3, 0, 1), ncol=3)
    A <- rbinom(n, size = 1, prob = pi0(W))
    if(null.sims) {
      Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
      psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
      theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
    } else {
      Y <- rbinom(n, size = 1, prob = mu0(A, W))
      psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
      theta0 <- var((mu0(1,W) - mu0(0,W)))
    }

    #print(control)
    output <- est.psi.prob(W, A, Y, est.g=TRUE)
    ret <- data.frame(type = 'Gamma.stat.2019', stat = output[1],
                      pvalue = output[2], quantile = -1)

    ret$n <- n
    ret$j <- j
    ret$seed <- seed

    ret$psi0 <- psi0
    ret$theta0 <- theta0

    return(ret)
  })
})
testing.rst.1.luedtke.h1 <- sim2019rst
save(testing.rst.1.luedtke.h1, file='testing.rst.1.luedtke.h1.RData')

testing.rst.1.luedtke.h2 <- sim2019rst
save(testing.rst.1.luedtke.h2, file='testing.rst.1.luedtke.h2.RData')

testing.rst.1.luedtke.h3 <- sim2019rst
save(testing.rst.1.luedtke.h3, file='testing.rst.1.luedtke.h3.RData')

testing.rst.1.luedtke.h4 <- sim2019rst
save(testing.rst.1.luedtke.h4, file='testing.rst.1.luedtke.h4.RData')

######################################
testing.rst.1.luedtke <- rbind(rbind(rbind(testing.rst.1.luedtke.h1, testing.rst.1.luedtke.h2),
            testing.rst.1.luedtke.h3),
      testing.rst.1.luedtke.h4)
gamma.summaries.1.luedtke <- ddply(subset(testing.rst.1.luedtke, (type %in% 'Gamma.stat.2019')),
                           .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           # quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))
testing.rst.2.luedtke <- sim2019rst.2
save(testing.rst.2.luedtke, file="testing.rst.2.luedtke.RData")
gamma.summaries.2.luedtke <- ddply(subset(testing.rst.2.luedtke, (type %in% 'Gamma.stat.2019')),
                                   .(n, type), summarize,
                                   na = sum(is.na(stat)),
                                   cnt = length(stat),
                                   # quantile.reject.rate = mean(stat>quantile),
                                   pvalue.reject.rate = mean(pvalue<0.05))
load("../HTE_Simulation_Result/gamma.summaries.2.RData")
testing.rst.2.plot <- rbind(gamma.summaries.2[c(1,2,4,5),c("type", "n", "pvalue.reject.rate")],
                   gamma.summaries.2.luedtke[,c("type", "n", "pvalue.reject.rate")])

ggplot(testing.rst.2.plot, aes(x=n, y=pvalue.reject.rate, group=type, col=type)) +
  geom_line(aes(linetype=type))+
  geom_point(aes(shape=type))
# our testing seem to have lower type I error

######################################
# simulation scenario 2 in Luedtke
est.psi.prob = function(W,A,Y,W.train=NULL,A.train=NULL,Y.train=NULL,sig.meth='var',est.g=TRUE){
  require('SuperLearner')
  SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet')
  n=length(A)

  if(is.null(W.train) | is.null(A.train) | is.null(Y.train)){
    W.train = W
    A.train = A
    Y.train = Y
  }

  # Estimate outcome regressions
  Qbar.est = SuperLearner(Y=Y.train,X=data.frame(W=W.train,A=A.train),newX=data.frame(W=rbind(W,W),A=rep(c(0,1),each=n)),SL.library=SL.library)
  Qbar.est.0 = Qbar.est$SL.predict[,1][1:n]
  Qbar.est.1 = Qbar.est$SL.predict[,1][(n+1):(2*n)]

  if(est.g){
    gg = SuperLearner(Y=A,
                      X=W,
                      SL.library=SL.library,
                      family='binomial')$SL.predict[,1]
    # prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
    # gg <- prop.reg$fitted.values
  } else {
    gg = g0(W)
  }

  # Plug-in estimate of blip
  R = Qbar.est.1 - Qbar.est.0
  S = rep(0,n)

  D.R = A/gg * (Y-Qbar.est.1) - (1-A)/(1-gg) * (Y-Qbar.est.0)
  D.S = rep(0,n)

  return(mmd.test(R,S,D.R,D.S,sig.meth=sig.meth))
}

mu0.sim2 <- function(a, w, beta) 1 + beta*a*(1+w[,2]^2) + w[,1] + w[,2]

simrst2.luedtke.h2 <- ldply(seq(-0.5, 0.4, 0.1), function(beta) {
  ldply(1:500, function(j) {
    # print(j)

    n <- 100
    if(j %% 100 == 0) cat(beta, n, j, '\n')
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)

    # beta <- 0.5
    W <- data.frame(W1=rnorm(n), W2=rnorm(n))
    A <- rbinom(n, size=1, prob=0.5)
    Y <- rnorm(n, mean=mu0.sim2(A, W, beta), sd=1)
    psi0 <- mean((mu0.sim2(1, W, beta) - mu0.sim2(0, W, beta))^2)

    #print(control)
    output <- est.psi.prob(W, A, Y, est.g=TRUE)
    ret <- data.frame(type = 'Gamma.2.luedtke', stat = output[1],
                      pvalue = output[2], quantile = -1)

    ret$n <- n
    ret$j <- j
    ret$seed <- seed
    ret$beta <- beta

    ret$psi0 <- psi0

    return(ret)
  })
})

# simrst2.luedtke.h1 <- simrst2.luedtke
# simrst2.luedtke <- rbind(simrst2.luedtke.h1, simrst2.luedtke.h2)
gamma.summaries.sim2.luedtke <- ddply(subset(simrst2.luedtke, (type %in% 'Gamma.2.luedtke')),
                                   .(n, type, beta), summarize,
                                   na = sum(is.na(stat)),
                                   cnt = length(stat),
                                   # quantile.reject.rate = mean(stat>quantile),
                                   pvalue.reject.rate = mean(pvalue<0.05))
save(simrst2.luedtke, file='simrst2.luedtke.RData')
# simrst2.luedtke[(simrst2.luedtke$beta=='0.3'),]
# simrst2.luedtke$beta <- as.numeric(simrst2.luedtke$beta)
# simrst2.luedtke[(simrst2.luedtke$beta==as.character(beta)), "seed"]

# beta <- 0.3
#
control <- list(pi.SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet'),
                mu.SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet'))
out.glm <- FALSE
ret <- data.frame()
for (beta in 0.2){
  for (j in 391:500){
    # if(j %% 100 == 0)
    n <- 100
    cat(beta, n, j, '\n')
    seed <- simrst2.luedtke[(simrst2.luedtke$beta==as.character(beta)) & (simrst2.luedtke$j==j), "seed"][1]
    set.seed(seed)

    # beta <- 0.5
    W <- data.frame(W1=rnorm(n), W2=rnorm(n))
    A <- rbinom(n, size=1, prob=0.5)
    Y <- rnorm(n, mean=mu0.sim2(A, W, beta), sd=1)
    psi0 <- mean((mu0.sim2(1, W, beta) - mu0.sim2(0, W, beta))^2)
    theta0 <- var((mu0.sim2(1, W, beta) - mu0.sim2(0, W, beta)))

    if (length(ret)>0){
      ret.tmp <- hteNullTest(Y, A, W, control = control, out.glm = out.glm, cov.var=FALSE)

      ret.tmp$n <- rep(n, 2)
      ret.tmp$j <- rep(j, 2)
      ret.tmp$seed <- rep(seed, 2)
      ret.tmp$beta <- rep(beta, 2)

      ret.tmp$psi0 <- rep(psi0, 2)
      ret.tmp$theta0 <- rep(theta0, 2)

      ret <- rbind(ret, ret.tmp)
    }
    else{
      ret <- hteNullTest(Y, A, W, control = control, out.glm = out.glm, cov.var=FALSE)

      ret$n <- rep(n, 2)
      ret$j <- rep(j, 2)
      ret$seed <- rep(seed, 2)
      ret$beta <- rep(beta, 2)

      ret$psi0 <- rep(psi0, 2)
      ret$theta0 <- rep(theta0, 2)
    }

  }
}
# simrst2 <- ret
# simrst2.h3 <- ret
# simrst2.h2 <- ret
# simrst2 <- rbind(simrst2.h2, simrst2.h3)
simrst2$beta <- as.character(simrst2$beta)
save(simrst2, file='simrst2.RData')
gamma.summaries.sim2 <- ddply(subset(simrst2, (type %in% 'Gamma.stat')),
                                      .(beta, type), summarize,
                                      na = sum(is.na(stat)),
                                      cnt = length(stat),
                                      # quantile.reject.rate = mean(stat>quantile),
                                      pvalue.reject.rate = mean(pvalue<0.05))
gamma.summaries.sim2$beta <- as.numeric(gamma.summaries.sim2$beta)
# plot
sim2.plot <- rbind(gamma.summaries.sim2[,c("type", "beta", "pvalue.reject.rate")],
      gamma.summaries.sim2.luedtke[,c("type", "beta", "pvalue.reject.rate")])

ggplot(sim2.plot, aes(x=beta, y=pvalue.reject.rate, group=type, col=type)) +
  geom_line(aes(linetype=type))+
  geom_point(aes(shape=type))

# our method has slightly higher power than luedtke's method when beta is far away from zero

# DEBUG
generate.data(null.sims)
testing.sim(100, 1:500, control = list(), null.sims=FALSE, out.glm=TRUE, generate.data)

n <- 100
beta <- 0.2
j <- 390
cat(beta, n, j, '\n')

seed <- simrst2.luedtke[(simrst2.luedtke$beta==as.character(beta)) & (simrst2.luedtke$j==j), "seed"][1]
set.seed(seed)

# beta <- 0.5
W <- data.frame(W1=rnorm(n), W2=rnorm(n))
A <- rbinom(n, size=1, prob=0.5)
Y <- rnorm(n, mean=mu0.sim2(A, W, beta), sd=1)
psi0 <- mean((mu0.sim2(1, W, beta) - mu0.sim2(0, W, beta))^2)
theta0 <- var((mu0.sim2(1, W, beta) - mu0.sim2(0, W, beta)))

control <- hte.measure.NullTest.control(control <- list(pi.SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet'),
                                                        mu.SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet')))
ret <- hteNullTest(Y, A, W, control = control, out.glm = FALSE, cov.var=FALSE)

# simulation for paper
generate.data <- function(n){
  beta <- c(0, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}

n <- 250
generate.data(n)

# analysis for simulated data
pi0.array <- c()
mu0.array <- c()
psi0 <- c()
theta0 <- c()
for (i in 1:500){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  pi0.array[i] <- mean(pi0(W))
  mu0.array[i] <- mean(mu0(A, W, beta))
  psi0[i]<- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0[i] <- var((mu0(1, W, beta) - mu0(0, W, beta)))
}

plot(density(pi0.array))
plot(density(mu0.array))
plot(density(psi0))
plot(density(theta0))

# how to name data sets?
# ests.sim.1.0.25.w.sl
# ests - for estimator or CI
# 1 - for simulated data function
# 0.25 - for beta parameters in generate.data()
# w - for wald type CI
# sl - for superlearner, i.e., out.glm=FALSE

data.dict <- list(ests.sim.1.0.25.w.sl=list(beta=c(0, 0.25)))
# 1) estimator simulation
control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
tm0 <- proc.time()
ests.sim.1.0.25.w.sl <-  ests.sim(c(100, 250, 500, 750, 1000), 1:500, control, out.glm=FALSE)
tm1 <- proc.time()
cat("tm1-tm0 = ", tm1-tm0, "\n")
save(ests.sim.1.0.25.w.sl, file="ests.sim.1.0.25.w.sl.RData")

ests.sim.1.0.25.w.sl

load("../HTE_Simulation_Result/ests.sim.1.0.25.w.sl.RData")
ret <- subset(ests.sim.1.0.25.w.sl,
              (type %in% c('psi.est', 'theta.est')))[,c("type", "est", "ll", "ul", "n", "j" , "seed")]
optional.ret <- subset(ests.sim.1.0.25.w.sl,
                       (type %in% c('psi.est', 'theta.est', 'psi.plug.in.est', 'psi.one.step.est',
                       'theta.plug.in.est', 'theta.one.step.est')))[,c("type", "est", "n", "j" , "seed")]
pars.ret <- subset(ests.sim.1.0.25.w.sl,
                       (type %in% 'pars'))[,c("n", "j" , "seed", "psi0", "theta0")]
ests.sim.1 <- join(ret, pars.ret, by=c("n", "j", "seed"))
ests.sim.2 <- join(optional.ret, pars.ret, by=c("n", "j", "seed"))

psi.summaries.1 <- ddply(subset(ests.sim.1, (type %in% 'psi.est')), .(n, type), summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         cnt = length(type),
                         bias = mean(est - psi0, na.rm=TRUE),
                         var = var(est, na.rm=TRUE),
                         mse = mean((est - psi0)^2, na.rm=TRUE))

psi.summaries.2 <- ddply(subset(ests.sim.2, (type %in% c('psi.est', 'psi.plug.in.est', 'psi.one.step.est'))),
                                .(n, type), summarize, na = sum(is.na(est)),
                         # coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         bias = mean(est - psi0, na.rm=TRUE),
                         var = var(est, na.rm=TRUE),
                         mse = mean((est - psi0)^2, na.rm=TRUE))

library(ggplot2)
p1 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, coverage, color=type)) +
  geom_hline(yintercept=.95)

p2 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, sqrt(n) * bias, color=type))

p3 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, n * mse, color=type))

p4 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, n * var, color=type)) +
  coord_cartesian(ylim=c(0,1))

library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)

control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
generate.data <- function(n){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
ests.sim <-  function(n_range, j_range, control, generate_func, out.glm=FALSE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      cat(n, j, '\n')
      # if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      # generate simulated data
      simulated.data <- do.call(generate_func, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W

      # using glm to estimate
      if (out.glm){
        prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
        pi.hat <- prop.reg$fitted.values
        AW <- cbind(A, data.frame(W))
        mu.reg <- glm(Y ~ ., data=AW, family='binomial')
        mu.hat <- mu.reg$fitted.values
        mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                                   data.frame(W)), type = 'response')
        mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                                   data.frame(W)), type = 'response')
        mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
        control$pi.hat = pi.hat
        control$mu.hats = mu.hats
      }

      # using SuperLearner to estimate
      est.ret <- htem.estimator(A, W, Y, control = control)

      ret <- label.result(est.ret$ret, n, j, seed)
      optional.ret <- NULL
      if(!is.null(est.ret$optional.ret)){
        optional.ret <- label.result(est.ret$optional.ret, n, j, seed)
      }

      parameters.ret <- data.frame(type='pars', n=n, j=j, seed=seed,
                                   psi0=simulated.data$psi0,
                                   theta0=simulated.data$theta0)

      est.sim.ret <- bind_rows(ret, optional.ret, parameters.ret)
      return(est.sim.ret)
    })
  })

  return(ests)
}
tm0 <- proc.time()
ests.sim.1.25.25.w.sl <-  ests.sim(c(100, 250, 500, 750, 1000), 1:500, control, generate_func=generate.data, out.glm=FALSE)
tm1 <- proc.time()
cat("tm1-tm0 = ", tm1-tm0, "\n")
save(ests.sim.1.25.25.w.sl, file="ests.sim.1.25.25.w.sl.RData")

control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
generate.data <- function(n){
  beta <- c(0.25, 0)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
ests.sim <-  function(n_range, j_range, control, generate_func, out.glm=FALSE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      cat(n, j, '\n')
      # if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      # generate simulated data
      simulated.data <- do.call(generate_func, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W

      # using glm to estimate
      if (out.glm){
        prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
        pi.hat <- prop.reg$fitted.values
        AW <- cbind(A, data.frame(W))
        mu.reg <- glm(Y ~ ., data=AW, family='binomial')
        mu.hat <- mu.reg$fitted.values
        mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                                   data.frame(W)), type = 'response')
        mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                                   data.frame(W)), type = 'response')
        mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
        control$pi.hat = pi.hat
        control$mu.hats = mu.hats
      }

      # using SuperLearner to estimate
      est.ret <- htem.estimator(A, W, Y, control = control)

      ret <- label.result(est.ret$ret, n, j, seed)
      optional.ret <- NULL
      if(!is.null(est.ret$optional.ret)){
        optional.ret <- label.result(est.ret$optional.ret, n, j, seed)
      }

      parameters.ret <- data.frame(type='pars', n=n, j=j, seed=seed,
                                   psi0=simulated.data$psi0,
                                   theta0=simulated.data$theta0)

      est.sim.ret <- bind_rows(ret, optional.ret, parameters.ret)
      return(est.sim.ret)
    })
  })

  return(ests)
}
tm0 <- proc.time()
ests.sim.1.25.0.w.sl <-  ests.sim(c(100, 250, 500, 750, 1000), 1:500, control, generate_func=generate.data, out.glm=FALSE)
tm1 <- proc.time()
cat("tm1-tm0 = ", tm1-tm0, "\n")
save(ests.sim.1.25.0.w.sl, file="ests.sim.1.25.0.w.sl.RData")

control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
generate.data <- function(n){
  beta <- c(0, 0)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
ests.sim <-  function(n_range, j_range, control, generate_func, out.glm=FALSE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      cat(n, j, '\n')
      # if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      # generate simulated data
      simulated.data <- do.call(generate_func, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W

      # using glm to estimate
      if (out.glm){
        prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
        pi.hat <- prop.reg$fitted.values
        AW <- cbind(A, data.frame(W))
        mu.reg <- glm(Y ~ ., data=AW, family='binomial')
        mu.hat <- mu.reg$fitted.values
        mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                                   data.frame(W)), type = 'response')
        mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                                   data.frame(W)), type = 'response')
        mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
        control$pi.hat = pi.hat
        control$mu.hats = mu.hats
      }

      # using SuperLearner to estimate
      est.ret <- htem.estimator(A, W, Y, control = control)

      ret <- label.result(est.ret$ret, n, j, seed)
      optional.ret <- NULL
      if(!is.null(est.ret$optional.ret)){
        optional.ret <- label.result(est.ret$optional.ret, n, j, seed)
      }

      parameters.ret <- data.frame(type='pars', n=n, j=j, seed=seed,
                                   psi0=simulated.data$psi0,
                                   theta0=simulated.data$theta0)

      est.sim.ret <- bind_rows(ret, optional.ret, parameters.ret)
      return(est.sim.ret)
    })
  })

  return(ests)
}
tm0 <- proc.time()
ests.sim.1.0.0.w.sl <-  ests.sim(c(100, 250, 500, 750, 1000), 1:500, control, generate_func=generate.data, out.glm=FALSE)
tm1 <- proc.time()
cat("tm1-tm0 = ", tm1-tm0, "\n")
save(ests.sim.1.0.0.w.sl, file="ests.sim.1.0.0.w.sl.RData")

######################################################################
control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
load("HTEM_Simulation/Testing/testing.sim.1.25.75.sl.RData")
generate.data <- function(n){
  beta <- c(0.25, 0.75)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
seed.data <- create.seed.dict(testing.sim.1.25.75.sl)
cl <- makeCluster(getOption("cl.cores", 5))
clusterEvalQ(cl, {
  library(SuperLearner)
  library(plyr)
  library(dplyr)
})
clusterExport(cl, c("label.result", "htem.estimator", "hteNullTest", "hte.measure.NullTest.control"),
              envir=environment())
ests.sim.1.25.75.w.sl.c <- data.frame()
for (i in c(100, 250, 500, 750, 1000)){
  cat(i, '\n')
  system.time(d1 <- do.call(rbind.data.frame, parLapply(cl, 1:500, ests.sim.testing.correction, 
                                                        n_range=i, 
                                                        control=control, 
                                                        generate_func=generate.data, 
                                                        seed.data=seed.data, out.glm=FALSE)))
  ests.sim.1.25.75.w.sl.c <- rbind(ests.sim.1.25.75.w.sl.c, 
                                   d1)
}
stopCluster(cl)
save(ests.sim.1.25.75.w.sl.c, file="ests.sim.1.25.75.w.sl.c.RData")

###############################


# 2) Bootstrap simulation
control = list(est.type = list('psi.est'='hybrid', 'theta.est'='hybrid'),
               conf.int = TRUE,
               conf.int.type = 'boot',
               n.boot = 500)
generate.data <- function(n){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
ests.sim <-  function(n_range, j_range, control, generate_func, out.glm=FALSE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      cat(n, j, '\n')
      # if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      # generate simulated data
      simulated.data <- do.call(generate_func, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W

      # using glm to estimate
      if (out.glm){
        prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
        pi.hat <- prop.reg$fitted.values
        AW <- cbind(A, data.frame(W))
        mu.reg <- glm(Y ~ ., data=AW, family='binomial')
        mu.hat <- mu.reg$fitted.values
        mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                                   data.frame(W)), type = 'response')
        mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                                   data.frame(W)), type = 'response')
        mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
        control$pi.hat = pi.hat
        control$mu.hats = mu.hats
      }

      # using SuperLearner to estimate
      est.ret <- htem.estimator(A, W, Y, control = control)

      ret <- label.result(est.ret$ret, n, j, seed)
      optional.ret <- NULL
      if(!is.null(est.ret$optional.ret)){
        optional.ret <- label.result(est.ret$optional.ret, n, j, seed)
      }

      parameters.ret <- data.frame(type='pars', n=n, j=j, seed=seed,
                                   psi0=simulated.data$psi0,
                                   theta0=simulated.data$theta0)

      est.sim.ret <- bind_rows(ret, optional.ret, parameters.ret)
      return(est.sim.ret)
    })
  })

  return(ests)
}
tm0 <- proc.time()
ests.sim.1.25.25.b.sl.1000 <-  ests.sim(c(1000), 1:500, control,
                                       generate_func=generate.data, out.glm=FALSE)
tm1 <- proc.time()
cat("tm1-tm0 = ", tm1-tm0, "\n")
save(ests.sim.1.25.25.b.sl.1000, file="ests.sim.1.25.25.b.sl.1000.RData")


## for other beta settings
# beta(0,0)
load("HTEM_Simulation/Testing/testing.sim.1.0.0.sl.RData")
generate.data <- function(n){
  beta <- c(0, 0)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
seed.data <- create.seed.dict(testing.sim.1.0.0.sl)
cl <- makeCluster(getOption("cl.cores", 5))
clusterEvalQ(cl, {
  library(SuperLearner)
  library(plyr)
  library(dplyr)
})
clusterExport(cl, c("label.result", "htem.estimator", "hteNullTest", "hte.measure.NullTest.control"),
              envir=environment())
ests.sim.1.0.0.b.sl.c <- data.frame()
for (i in c(100, 250, 500, 750, 1000)){
  system.time(d1 <- do.call(rbind.data.frame, parLapply(cl, 1:500, ests.sim.testing.correction, 
                                            n_range=i, 
                                            control=control, 
                                            generate_func=generate.data, 
                                            seed.data=seed.data, out.glm=FALSE)))
  ests.sim.1.0.0.b.sl.c <- rbind(ests.sim.1.0.0.b.sl.c, 
                                 d1)
}
stopCluster(cl)
save(ests.sim.1.0.0.b.sl.c, file="ests.sim.1.0.0.b.sl.c.RData")

load("HTEM_Simulation/Testing/testing.sim.1.0.25.sl.RData")
generate.data <- function(n){
  beta <- c(0, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
seed.data <- create.seed.dict(testing.sim.1.0.25.sl)
cl <- makeCluster(getOption("cl.cores", 5))
clusterEvalQ(cl, {
  library(SuperLearner)
  library(plyr)
  library(dplyr)
})
clusterExport(cl, c("label.result", "htem.estimator", "hteNullTest", "hte.measure.NullTest.control"),
              envir=environment())
ests.sim.1.0.25.b.sl.c <- data.frame()
for (i in c(100, 250, 500, 750, 1000)){
  cat(i, '\n')
  system.time(d1 <- do.call(rbind.data.frame, parLapply(cl, 1:500, ests.sim.testing.correction, 
                                                        n_range=i, 
                                                        control=control, 
                                                        generate_func=generate.data, 
                                                        seed.data=seed.data, out.glm=FALSE)))
  ests.sim.1.0.25.b.sl.c <- rbind(ests.sim.1.0.25.b.sl.c, 
                                 d1)
}
stopCluster(cl)
save(ests.sim.1.0.25.b.sl.c, file="ests.sim.1.0.25.b.sl.c.RData")

load("HTEM_Simulation/Testing/testing.sim.1.25.75.sl.RData")
generate.data <- function(n){
  beta <- c(0.25, 0.75)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
seed.data <- create.seed.dict(testing.sim.1.25.75.sl)
cl <- makeCluster(getOption("cl.cores", 5))
clusterEvalQ(cl, {
  library(SuperLearner)
  library(plyr)
  library(dplyr)
})
clusterExport(cl, c("label.result", "htem.estimator", "hteNullTest", "hte.measure.NullTest.control"),
              envir=environment())
ests.sim.1.25.75.b.sl.c <- data.frame()
for (i in c(100, 250, 500, 750, 1000)){
  cat(i, '\n')
  system.time(d1 <- do.call(rbind.data.frame, parLapply(cl, 1:500, ests.sim.testing.correction, 
                                                        n_range=i, 
                                                        control=control, 
                                                        generate_func=generate.data, 
                                                        seed.data=seed.data, out.glm=FALSE)))
  ests.sim.1.25.75.b.sl.c <- rbind(ests.sim.1.25.75.b.sl.c, 
                                  d1)
}
stopCluster(cl)
save(ests.sim.1.25.75.b.sl.c, file="ests.sim.1.25.75.b.sl.c.RData")

ret <- subset(ests.sim.1.25.75.b.sl.c,
              (type %in% c('psi.est', 'theta.est')))[,c("type", "est", "ll", "ul", "n", "j" , "seed")]
pars.ret <- subset(ests.sim.1.25.75.b.sl.c,
                   (type %in% 'pars'))[,c("n", "j" , "seed", "psi0", "theta0")]
ests.sim.1.25.75.b.sl.ret <- join(ret, pars.ret, by=c("n", "j", "seed"))

ci.correction.b <- join(subset(ests.sim.1.25.75.b.sl.ret, (type %in% 'psi.est')),
                        subset(testing.sim.1.25.75.sl, (type %in% 'Gamma.stat'))[,c("pvalue", "n", "j" , "seed")], by=c("n", "j", "seed"))

psi.b.summaries <- ddply(ci.correction.b %>% mutate(ll.c = ifelse(pvalue>0.05, 0, ll)), .(n, type), 
                         summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         coverage.c = mean(ll.c <= psi0 & psi0 <= ul, na.rm=TRUE),
                         cnt = length(type))
library(knitr)
kable(psi.b.summaries, caption="Bootstrap CI coverage for psi0 with correction")

ci.correction.b <- join(subset(ests.sim.1.25.75.b.sl.ret, (type %in% 'theta.est')),
                        subset(testing.sim.1.25.75.sl, (type %in% 'Omega.stat'))[,c("pvalue", "n", "j" , "seed")], by=c("n", "j", "seed"))

theta.b.summaries <- ddply(ci.correction.b %>% mutate(ll.c = ifelse(pvalue>0.05, 0, ll)), .(n, type), 
                           summarize, na = sum(is.na(est)),
                           coverage = mean(ll <= theta0 & theta0 <= ul, na.rm=TRUE),
                           coverage.c = mean(ll.c <= theta0 & theta0 <= ul, na.rm=TRUE),
                           cnt = length(type))
library(knitr)
kable(theta.b.summaries, caption="Bootstrap CI coverage for theta0 with correction")



# 3) testing simulation
tm0 <- proc.time()
testing.sim.1.0.25.sl <- testing.sim(1000, 1:10, control = list(), out.glm = FALSE)
tm1 <- proc.time()
cat("tm1-tm0 = ", tm1-tm0, "\n")
save(ests.sim.1.0.25.w.sl, file="ests.sim.1.0.25.w.sl.RData")
tmp.500 <- testing.sim.1.0.25.sl
testing.sim.1.0.25.sl
gamma.summaries.1 <- ddply(subset(testing.sim.1.0.25.sl, (type %in% 'Gamma.stat')),
                                   .(n, type), summarize,
                                   na = sum(is.na(stat)),
                                   cnt = length(stat),
                                   # quantile.reject.rate = mean(stat>quantile),
                                   pvalue.reject.rate = mean(pvalue<0.05))
omega.summaries.1 <- ddply(subset(testing.sim.1.0.25.sl, (type %in% 'Omega.stat')),
                           .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           # quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))

load("../HTE_Simulation_Result/testing.sim.1.25.0.sl.100.RData")
testing.sim.1.25.0.sl.100
gamma.summaries.1 <- ddply(subset(testing.sim.1.25.0.sl.100, (type %in% 'Gamma.stat')),
                           .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           # quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))
omega.summaries.1 <- ddply(subset(testing.sim.1.25.0.sl.100, (type %in% 'Omega.stat')),
                           .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           # quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))

# try parallel computing
generate.data <- function(n){
  beta <- c(0, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}

# testing.sim.1.0.25.sl <- testing.sim(1000, 1:10, control = list(), generate_func=generate.data, out.glm = FALSE)
test <- mclapply(1:10, testing.sim, n_range=100, control = list(), generate_func=generate.data, out.glm = FALSE)

library(parallel)
numCores <- detectCores()
tm0 <- proc.time()
testing.sim.1.0.25.sl.100 <- mclapply(1:500, testing.sim, n_range=100, control = list(), generate_func=generate.data, out.glm = FALSE, mc.cores = 4)
tm1 <- proc.time()
cat("tm1-tm0 = ", tm1-tm0, "\n")

tm0 <- proc.time()
tm <- system.time(
  testing.sim.1.0.25.sl.500 <- mclapply(1:500, testing.sim, n_range=500, control = list(), generate_func=generate.data, out.glm = FALSE, mc.cores = 4)
)
cat("tm1-tm0 = ", tm1-tm0, "\n")
cat("tm1-tm0 = ", tm, "\n")


tm <- system.time(
  testing.sim.1.0.25.sl.1000 <- mclapply(1:500, testing.sim, n_range=1000, control = list(), generate_func=generate.data, out.glm = FALSE, mc.cores = 4)
)
cat("tm1-tm0 = ", tm, "\n")

save(testing.sim.1.0.25.sl.100, file="testing.sim.1.0.25.sl.100.RData")
save(testing.sim.1.0.25.sl.500, file="testing.sim.1.0.25.sl.500.RData")
save(testing.sim.1.0.25.sl.1000, file="testing.sim.1.0.25.sl.1000.RData")


tm <- system.time(
  testing.sim.1.0.25.sl.250 <- mclapply(1:500, testing.sim, n_range=250, control = list(), generate_func=generate.data, out.glm = FALSE, mc.cores = 4)
)
cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.0.25.sl.250, file="testing.sim.1.0.25.sl.250.RData")
tm <- system.time(
  testing.sim.1.0.25.sl.750 <- mclapply(1:500, testing.sim, n_range=750, control = list(), generate_func=generate.data, out.glm = FALSE, mc.cores = 4)
)
cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.0.25.sl.750, file="testing.sim.1.0.25.sl.750.RData")


#######################################
generate.data <- function(n){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
tm <- system.time(
  testing.sim.1.25.25.sl.750 <- mclapply(1:500, testing.sim, n_range=750, control = list(),
                                         generate_func=generate.data, out.glm = FALSE, mc.cores = 6)
)

cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.25.25.sl.750, file="testing.sim.1.25.25.sl.750.RData")

tm <- system.time(
  testing.sim.1.25.25.sl.250 <- mclapply(1:500, testing.sim, n_range=250, control = list(),
                                         generate_func=generate.data, out.glm = FALSE, mc.cores = 6)
)

cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.25.25.sl.250, file="testing.sim.1.25.25.sl.250.RData")

tm <- system.time(
  testing.sim.1.25.25.sl.100 <- mclapply(1:500, testing.sim, n_range=100, control = list(),
                                         generate_func=generate.data, out.glm = FALSE, mc.cores = 6)
)
cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.25.25.sl.100, file="testing.sim.1.25.25.sl.100.RData")

generate.data <- function(n){
  beta <- c(0, 0)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}

tm <- system.time(
  testing.sim.1.0.0.sl.250 <- mclapply(1:500, testing.sim, n_range=250, control = list(),
                                         generate_func=generate.data, out.glm = FALSE, mc.cores = 6)
)

cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.0.0.sl.250, file="testing.sim.1.0.0.sl.250.RData")

tm <- system.time(
  testing.sim.1.0.0.sl.100 <- mclapply(1:500, testing.sim, n_range=100, control = list(),
                                         generate_func=generate.data, out.glm = FALSE, mc.cores = 6)
)
cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.0.0.sl.100, file="testing.sim.1.0.0.sl.100.RData")

tm <- system.time(
  testing.sim.1.0.0.sl.750 <- mclapply(1:500, testing.sim, n_range=750, control = list(),
                                         generate_func=generate.data, out.glm = FALSE, mc.cores = 6)
)

cat("tm1-tm0 = ", tm, "\n")
save(testing.sim.1.0.0.sl.750, file="testing.sim.1.0.0.sl.750.RData")
# organizing simulation result
d1 <- do.call(rbind.data.frame, testing.sim.1.0.25.sl.100)
d2 <- do.call(rbind.data.frame, testing.sim.1.0.25.sl.250)
d3 <- do.call(rbind.data.frame, testing.sim.1.0.25.sl.500)
d4 <- do.call(rbind.data.frame, testing.sim.1.0.25.sl.750)
d5 <- do.call(rbind.data.frame, testing.sim.1.0.25.sl.1000)
testing.sim.1.0.25.sl <- bind_rows(d1, d2, d3, d4, d5)
save(testing.sim.1.0.25.sl, file="testing.sim.1.0.25.sl.RData")

for (index in c(100, 250, 500, 750, 1000)){
  load(paste("../HTE_Simulation_Result/testing.sim.1.25.0.sl.", as.character(index), ".RData", sep=""))
}
testing.sim.1.25.0.sl <- data.frame()
for (index in c(100, 250, 500, 750, 1000)){
  testing.sim.1.25.0.sl <- rbind(testing.sim.1.25.0.sl,
        get(paste("testing.sim.1.25.0.sl.", as.character(index), sep="")))
}
save(testing.sim.1.25.0.sl, file="testing.sim.1.25.0.sl.RData")

for (index in c(100, 250, 500, 750, 1000)){
  load(paste("../HTE_Simulation_Result/testing.sim.1.25.25.sl.", as.character(index), ".RData", sep=""))
}


d1 <- do.call(rbind.data.frame, testing.sim.1.25.25.sl.100)
d2 <- do.call(rbind.data.frame, testing.sim.1.25.25.sl.250)
# d3 <- do.call(rbind.data.frame, testing.sim.1.25.25.sl.500)
d4 <- do.call(rbind.data.frame, testing.sim.1.25.25.sl.750)
# d5 <- do.call(rbind.data.frame, testing.sim.1.25.25.sl.1000)
testing.sim.1.25.25.sl <- bind_rows(d1, d2, testing.sim.1.25.25.sl.500, d4, testing.sim.1.25.25.sl.1000)

save(testing.sim.1.25.25.sl, file="testing.sim.1.25.25.sl.RData")


for (index in c(100, 250, 500, 750, 1000)){
  load(paste("../HTE_Simulation_Result/testing.sim.1.0.0.sl.", as.character(index), ".RData", sep=""))
}


d1 <- do.call(rbind.data.frame, testing.sim.1.0.0.sl.100)
d2 <- do.call(rbind.data.frame, testing.sim.1.0.0.sl.250)
# d3 <- do.call(rbind.data.frame, testing.sim.1.25.25.sl.500)
d4 <- do.call(rbind.data.frame, testing.sim.1.0.0.sl.750)
# d5 <- do.call(rbind.data.frame, testing.sim.1.25.25.sl.1000)
testing.sim.1.0.0.sl <- bind_rows(d1, d2, testing.sim.1.0.0.sl.500, d4, testing.sim.1.0.0.sl.1000)

save(testing.sim.1.0.0.sl, file="testing.sim.1.0.0.sl.RData")


for (index in c(100, 250, 500, 750, 1000)){
  load(paste("ests.sim.1.25.25.b.sl.", as.character(index), ".RData", sep=""))
}
ests.sim.1.25.25.b.sl <- data.frame()
for (index in c(100, 250, 500, 750, 1000)){
  ests.sim.1.25.25.b.sl <- rbind(ests.sim.1.25.25.b.sl,
                                 get(paste("ests.sim.1.25.25.b.sl.", as.character(index), sep="")))
}
save(ests.sim.1.25.25.b.sl, file="ests.sim.1.25.25.b.sl.RData")

# parallel computing on windows system
generate.data <- function(n){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
library(devtools)
# load_all()
devtools::install_github("ruihu51/CausalSim")



require(CausalSim)
generate.data <- function(n){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
hte.measure.NullTest.control <- function(control){
  control.default = list(pi.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                         mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                         est.type = list('psi.est'='hybrid', 'theta.est'='hybrid'),
                         conf.int = FALSE,
                         conf.int.type = 'Wald',
                         conf.level = 0.95,
                         n.boot = 500,
                         verbose = FALSE)
  control.names <- names(control)
  if(!is.null(control.names)) {
    for (name in control.names) {
      control.default[[name]] <- control[[name]]
    }
  }
  return(control.default)
}

# Function to label result
label.result <- function(ret, n, j, seed){
  nrows <- dim(ret)[1]
  ret$n <- rep(n, nrows)
  ret$j <- rep(j, nrows)
  ret$seed <- rep(seed, nrows)
  return(ret)
}

hteNullTest <- function(Y, A, W, control = list(), out.glm=FALSE, cov.var=FALSE) {
  
  # update control parameters
  control <- hte.measure.NullTest.control(control)
  n = length(A)
  
  if (out.glm){
    prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
    pi.hat <- prop.reg$fitted.values
    AW <- cbind(A, data.frame(W))
    mu.reg <- glm(Y ~ ., data=AW, family='binomial')
    mu.hat <- mu.reg$fitted.values
    mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                               data.frame(W)), type = 'response')
    mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                               data.frame(W)), type = 'response')
    mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
    control$pi.hat = pi.hat
    control$mu.hats = mu.hats
  }
  
  if(control$verbose) cat("Estimating...\n")
  tm0 <- proc.time()
  # estimated propensity
  if(is.null(control$pi.hat)){
    prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                             newX = data.frame(W),
                             SL.library = control$pi.SL.library,
                             family = binomial(),
                             obsWeights=rep(1,n),
                             id=1:n)
    control$pi.hat <- prop.reg$SL.predict
  }
  
  # estimated outcome regression
  if(is.null(control$mu.hats)){
    AW <- cbind(A, data.frame(W))
    if(length(setdiff(Y, c(0,1))) == 0) {
      family = 'binomial'
    } else {
      family = 'gaussian'
    }
    mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                           newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                           SL.library = control$mu.SL.library,
                           family = family,
                           obsWeights=rep(1,n),
                           id=1:n)
    control$mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
  }
  
  control$mu.hat <- A * control$mu.hats$mu1 + (1-A) * control$mu.hats$mu0
  
  # estimated tau
  tau.hat <- control$mu.hats$mu1 - control$mu.hats$mu0
  Z.hat <- (2*A - 1) / (A * control$pi.hat + (1-A) * (1-control$pi.hat))
  
  # estimated theta
  gamma.hat <- mean(tau.hat)
  
  w.ecdf <- function(w) {
    vec.leq <- function(x,y) prod(x <= y)
    return(mean(apply(W, 1, vec.leq, y = w)))
  }
  u.vals <- apply(W, 1, w.ecdf)
  # primitive function
  if(control$verbose) cat("Computing Gamma and Omega...\n")
  tm1 <- proc.time()
  # n.new * 1 vector
  w.vals <- W
  Gamma.w.vals <- apply(w.vals, 1, function(w0)
    mean(apply(W, 1, function(x) prod(x <= w0))*tau.hat))
  Omega.w.vals <- Gamma.w.vals - gamma.hat * u.vals
  
  # nonparametric EIF
  # n * n.new matrix
  vec.eq <- function(x,y) prod(x == y)
  eif.Gamma <- apply(w.vals, 1, function(w0) {
    (apply(W, 1, function(x) prod(x <= w0))) * (Z.hat * (Y - control$mu.hat) + tau.hat) -
      Gamma.w.vals[which(as.logical(apply(W, 1, vec.eq, y = w0)))]
  })
  eif.Omega <- apply(w.vals, 1, function(w0) {
    (apply(W, 1, function(x) prod(x <= w0)) - w.ecdf(w0)) * (Z.hat * (Y - control$mu.hat) + tau.hat - gamma.hat) -
      Omega.w.vals[which(as.logical(apply(W, 1, vec.eq, y = w0)))]
  })
  
  # one-step estimators
  # n.new * 1 vector
  Gamma.os.est <- colMeans(eif.Gamma) + Gamma.w.vals
  Omega.os.est <- colMeans(eif.Omega) + Omega.w.vals
  
  # testing procedure
  if(control$verbose) cat("Computing statistics...\n")
  tm2 <- proc.time()
  # test statistics
  Gamma.stat <- n^(1/2)*max(abs(Gamma.os.est))
  Omega.stat <- n^(1/2)*max(abs(Omega.os.est))
  
  # covariance matrices
  n.new <- dim(w.vals)[1]
  if (cov.var){
    Gamma.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
      mean(eif.Gamma[,s] * eif.Gamma[,t])
    }))
    Omega.cov.var <- sapply(1:n.new, function(s) sapply(1:n.new, function(t) {
      mean(eif.Omega[,s] * eif.Omega[,t])
    }))
    
    # quantiles
    tm3 <- proc.time()
    Gamma.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Gamma.cov.var)
    Gamma.epsilon.stats <- apply(Gamma.epsilon, 1, function(x) {max(abs(x))})
    Omega.epsilon <- rmvnorm(n=control$n.boot, mean=rep(0, n.new), sigma = Omega.cov.var)
    Omega.epsilon.stats <- apply(Omega.epsilon, 1, function(x) {max(abs(x))})
  }else{
    tm3 <- proc.time()
    Gamma.epsilon.stats <- replicate(control$n.boot, max(abs(t(eif.Gamma)%*%rnorm(n.new, 0, 1)/sqrt(n))))
    Omega.epsilon.stats <- replicate(control$n.boot, max(abs(t(eif.Omega)%*%rnorm(n.new, 0, 1)/sqrt(n))))
  }
  
  Gamma.pvalue <- mean(Gamma.epsilon.stats > Gamma.stat)
  Gamma.quantile <- unname(quantile(Gamma.epsilon.stats, control$conf.level))
  Omega.pvalue <- mean(Omega.epsilon.stats > Omega.stat)
  Omega.quantile <- unname(quantile(Omega.epsilon.stats, control$conf.level))
  
  ret <- data.frame(type = 'Gamma.stat', stat = Gamma.stat, pvalue = Gamma.pvalue,
                    quantile = Gamma.quantile)
  ret <- rbind(ret,
               data.frame(type = 'Omega.stat', stat = Omega.stat, pvalue = Omega.pvalue,
                          quantile = Omega.quantile))
  tm4 <- proc.time()
  
  if (control$verbose){
    cat("tm1-tm0 = ", tm1-tm0, "\n")
    cat("tm2-tm1 = ", tm2-tm1, "\n")
    cat("tm3-tm2 = ", tm3-tm2, "\n")
    cat("tm4-tm3 = ", tm4-tm3, "\n")
  }
  
  ret
  
}
testing.sim <-  function(n_range, j_range, control, generate_func, out.glm=FALSE){
  ret <- data.frame()
  for (n in n_range){
    for (j in j_range){
      # if(j %% 100 == 0)
      cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)
      
      # generate simulated data
      simulated.data <- do.call(generate_func, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W
      
      if (length(ret)>0){
        ret.tmp <- hteNullTest(Y, A, W, control = control, out.glm = out.glm)
        
        ret.tmp$n <- rep(n, 2)
        ret.tmp$j <- rep(j, 2)
        ret.tmp$seed <- rep(seed, 2)
        
        ret.tmp$psi0 <- rep(simulated.data$psi0, 2)
        ret.tmp$theta0 <- rep(simulated.data$theta0, 2)
        
        ret <- rbind(ret, ret.tmp)
      }
      else{
        ret <- hteNullTest(Y, A, W, control = control, out.glm = out.glm)
        
        ret$n <- rep(n, 2)
        ret$j <- rep(j, 2)
        ret$seed <- rep(seed, 2)
        
        ret$psi0 <- rep(simulated.data$psi0, 2)
        ret$theta0 <- rep(simulated.data$theta0, 2)
      }
    }
  }
  
  return(ret)
}
library(parallel)
install.packages("doParallel")
library(doParallel)
cl <- makeCluster(getOption("cl.cores", 4))
clusterEvalQ(cl, library(SuperLearner))
clusterExport(cl, c("hteNullTest", "hte.measure.NullTest.control"), 
              envir=environment())
system.time({
  res <- parLapply(cl, 1:100, testing.sim, n_range=500, control=list(), generate_func=generate.data,
                   out.glm=FALSE)
})
res
stopCluster(cl)


###########################################################
# try to make theta larger
n <- 500
pi0.array <- c()
mu0.array <- c()
psi0 <- c()
theta0 <- c()

beta <- c(0.25, 0.75)
for (i in 1:500){
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  pi0.array[i] <- mean(pi0(W))
  mu0.array[i] <- mean(mu0(A, W, beta))
  psi0[i]<- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0[i] <- var((mu0(1, W, beta) - mu0(0, W, beta)))
}

data.frame(beta=paste("(",paste(beta, collapse=", "),")", sep=""),
                                   pi0=mean(pi0.array), mu0=mean(mu0.array),
                                   psi0=mean(psi0), theta0=mean(theta0))

# beta <- c(0.25, 0.75)
generate.data <- function(n){
  beta <- c(0.25, 0.75)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
library(parallel)
library(doParallel)
cl <- makeCluster(getOption("cl.cores", 7))
clusterEvalQ(cl, library(SuperLearner))
clusterExport(cl, c("hteNullTest", "hte.measure.NullTest.control"), 
              envir=environment())
system.time({
  testing.sim.1.25.75.sl.500 <- parLapply(cl, 1:500, testing.sim, n_range=500, 
                                         control=list(), generate_func=generate.data, out.glm=FALSE)
})
system.time({
  testing.sim.1.25.75.sl.1000 <- parLapply(cl, 1:500, testing.sim, n_range=1000, 
                                          control=list(), generate_func=generate.data, out.glm=FALSE)
})
system.time({
  testing.sim.1.25.75.sl.750 <- parLapply(cl, 1:500, testing.sim, n_range=750, 
                                          control=list(), generate_func=generate.data, out.glm=FALSE)
})
system.time({
  testing.sim.1.25.75.sl.100 <- parLapply(cl, 1:500, testing.sim, n_range=100, 
                                          control=list(), generate_func=generate.data, out.glm=FALSE)
})
system.time({
  testing.sim.1.25.75.sl.250 <- parLapply(cl, 1:500, testing.sim, n_range=250, 
                                          control=list(), generate_func=generate.data, out.glm=FALSE)
})

d1 <- do.call(rbind.data.frame, testing.sim.1.25.75.sl.100)
d2 <- do.call(rbind.data.frame, testing.sim.1.25.75.sl.250)
d3 <- do.call(rbind.data.frame, testing.sim.1.25.75.sl.500)
d4 <- do.call(rbind.data.frame, testing.sim.1.25.75.sl.750)
d5 <- do.call(rbind.data.frame, testing.sim.1.25.75.sl.1000)
testing.sim.1.25.75.sl <- bind_rows(d1, d2, d3, d4, d5)
save(testing.sim.1.25.75.sl, file="testing.sim.1.25.75.sl.RData")
getwd()

omega.summaries <- ddply(subset(testing.sim.1.25.75.sl, (type %in% 'Omega.stat')),
                         .(n, type), summarize,
                         na = sum(is.na(stat)),
                         cnt = length(stat),
                         # quantile.reject.rate = mean(stat>quantile),
                         reject.rate = mean(pvalue<0.05))

########################################
# Wald-type confidence interval with correction
control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
generate.data <- function(n){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2
  
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)
  
  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  
  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
ests.sim.testing.correction <-  function(n_range, j_range, control, 
                                         generate_func, seed.data, out.glm=FALSE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      # cat(n, j, '\n')
      # if(j %% 100 == 0) cat(n, j, '\n')
      seed <- subset(seed.data, (sample_n==n))[j,"seed"]
      cat(n, j, seed, '\n')
      set.seed(seed)
      
      # generate simulated data
      simulated.data <- do.call(generate_func, list(n=n))
      Y <- simulated.data$Y
      A <- simulated.data$A
      W <- simulated.data$W
      
      # using glm to estimate
      if (out.glm){
        prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
        pi.hat <- prop.reg$fitted.values
        AW <- cbind(A, data.frame(W))
        mu.reg <- glm(Y ~ ., data=AW, family='binomial')
        mu.hat <- mu.reg$fitted.values
        mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                                   data.frame(W)), type = 'response')
        mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                                   data.frame(W)), type = 'response')
        mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
        control$pi.hat = pi.hat
        control$mu.hats = mu.hats
      }
      
      # using SuperLearner to estimate
      est.ret <- htem.estimator(A, W, Y, control = control)
      
      ret <- label.result(est.ret$ret, n, j, seed)
      optional.ret <- NULL
      if(!is.null(est.ret$optional.ret)){
        optional.ret <- label.result(est.ret$optional.ret, n, j, seed)
      }
      
      parameters.ret <- data.frame(type='pars', n=n, j=j, seed=seed,
                                   psi0=simulated.data$psi0,
                                   theta0=simulated.data$theta0)
      
      est.sim.ret <- bind_rows(ret, optional.ret, parameters.ret)
      return(est.sim.ret)
    })
  })
  
  return(ests)
}

create.seed.dict <- function(x){
  x %>%
    subset(type %in% "Gamma.stat") %>%
    select(c("n", "j", "seed")) %>%
    rename(sample_n=n, j=j, seed=seed)
}
seed.data <- create.seed.dict(testing.sim.1.25.25.sl)


load("HTEM_Simulation/Testing/testing.sim.1.25.25.sl.RData")
system.time(
  ests.sim.1.25.25.w.sl.c <-  ests.sim.testing.correction(c(100, 250, 500, 750, 1000), 1:500, 
                                       control, 
                                       generate_func=generate.data, 
                                       seed.data=testing.sim.1.25.25.sl, 
                                       out.glm=FALSE)
)

############################### 
library(parallel)
library(doParallel)
cl <- makeCluster(getOption("cl.cores", 7))
clusterEvalQ(cl, {
  library(SuperLearner)
  library(plyr)
  library(dplyr)
  })
clusterExport(cl, c("label.result", "htem.estimator", "hteNullTest", "hte.measure.NullTest.control"),
              envir=environment())
system.time({
  ests.sim.1.25.25.w.sl.c.1000 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                            n_range=1000, 
                                          control=control, 
                                          generate_func=generate.data, 
                                          seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.w.sl.c.750 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                            n_range=750, 
                                            control=control, 
                                            generate_func=generate.data, 
                                            seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.w.sl.c.500 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                           n_range=500, 
                                           control=control, 
                                           generate_func=generate.data, 
                                           seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.w.sl.c.250 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                           n_range=250, 
                                           control=control, 
                                           generate_func=generate.data, 
                                           seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.w.sl.c.100 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                           n_range=100, 
                                           control=control, 
                                           generate_func=generate.data, 
                                           seed.data=seed.data, out.glm=FALSE)
})
stopCluster(cl)
d1 <- do.call(rbind.data.frame, ests.sim.1.25.25.w.sl.c.100)
d2 <- do.call(rbind.data.frame, ests.sim.1.25.25.w.sl.c.250)
d3 <- do.call(rbind.data.frame, ests.sim.1.25.25.w.sl.c.500)
d4 <- do.call(rbind.data.frame, ests.sim.1.25.25.w.sl.c.750)
d5 <- do.call(rbind.data.frame, ests.sim.1.25.25.w.sl.c.1000)
ests.sim.1.25.25.w.sl.c <- bind_rows(d1, d2, d3, d4, d5)
save(ests.sim.1.25.25.w.sl.c, file="ests.sim.1.25.25.w.sl.c.RData")

# reconstruct
ret <- subset(ests.sim.1.25.25.w.sl.c,
              (type %in% c('psi.est', 'theta.est')))[,c("type", "est", "ll", "ul", "n", "j" , "seed")]
optional.ret <- subset(ests.sim.1.25.25.w.sl.c,
                       (type %in% c('psi.est', 'theta.est', 'psi.plug.in.est', 'psi.one.step.est',
                                    'theta.plug.in.est', 'theta.one.step.est')))[,c("type", "est", "n", "j" , "seed")]
pars.ret <- subset(ests.sim.1.25.25.w.sl.c,
                   (type %in% 'pars'))[,c("n", "j" , "seed", "psi0", "theta0")]
ests.sim.1.25.25.w.sl.ret <- join(ret, pars.ret, by=c("n", "j", "seed"))
ests.sim.1.25.25.w.sl.optional.ret <- join(optional.ret, pars.ret, by=c("n", "j", "seed"))

# combining testing result
# psi-wald
ci.correction.w <- join(subset(ests.sim.1.25.25.w.sl.ret, (type %in% 'psi.est')),
     subset(testing.sim.1.25.25.sl, (type %in% 'Gamma.stat'))[,c("pvalue", "n", "j" , "seed")], by=c("n", "j", "seed"))

psi.w.summaries <- ddply(ci.correction.w %>% mutate(ll.c = ifelse(pvalue>0.05, 0, ll)), .(n, type), 
                         summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         coverage.c = mean(ll.c <= psi0 & psi0 <= ul, na.rm=TRUE),
                         cnt = length(type))
library(knitr)
kable(psi.w.summaries, caption="Wald-type CI coverage for psi0 with correction")

ci.correction.w %>% 
  mutate(ll.c = ifelse(pvalue>0.05, 0, ll),
         psi.coverage = (ll <= psi0) & (psi0 <= ul),
         psi.coverage.c = (ll.c <= psi0) & (psi0 <= ul)) %>%
  filter((!psi.coverage)|(!psi.coverage.c)) %>%
  filter(n==1000) %>%
  select(c(type, est, ll, ul, ll.c, psi0, pvalue))

# psi-bootstrap
control = list(est.type = list('psi.est'='hybrid', 'theta.est'='hybrid'),
               conf.int = TRUE,
               conf.int.type = 'boot',
               n.boot = 500)
cl <- makeCluster(getOption("cl.cores", 7))
clusterEvalQ(cl, {
  library(SuperLearner)
  library(plyr)
  library(dplyr)
})
clusterExport(cl, c("label.result", "htem.estimator", "hteNullTest", "hte.measure.NullTest.control"),
              envir=environment())
system.time({
  ests.sim.1.25.25.b.sl.c.1000 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                            n_range=1000, 
                                            control=control, 
                                            generate_func=generate.data, 
                                            seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.b.sl.c.750 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                           n_range=750, 
                                           control=control, 
                                           generate_func=generate.data, 
                                           seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.b.sl.c.500 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                           n_range=500, 
                                           control=control, 
                                           generate_func=generate.data, 
                                           seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.b.sl.c.250 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                           n_range=250, 
                                           control=control, 
                                           generate_func=generate.data, 
                                           seed.data=seed.data, out.glm=FALSE)
})
system.time({
  ests.sim.1.25.25.b.sl.c.100 <- parLapply(cl, 1:500, ests.sim.testing.correction, 
                                           n_range=100, 
                                           control=control, 
                                           generate_func=generate.data, 
                                           seed.data=seed.data, out.glm=FALSE)
})
stopCluster(cl)
d1 <- do.call(rbind.data.frame, ests.sim.1.25.25.b.sl.c.100)
d2 <- do.call(rbind.data.frame, ests.sim.1.25.25.b.sl.c.250)
d3 <- do.call(rbind.data.frame, ests.sim.1.25.25.b.sl.c.500)
d4 <- do.call(rbind.data.frame, ests.sim.1.25.25.b.sl.c.750)
d5 <- do.call(rbind.data.frame, ests.sim.1.25.25.b.sl.c.1000)
ests.sim.1.25.25.b.sl.c <- bind_rows(d1, d2, d3, d4, d5)
save(ests.sim.1.25.25.b.sl.c, file="ests.sim.1.25.25.b.sl.c.RData")

ret <- subset(ests.sim.1.25.25.b.sl.c,
              (type %in% c('psi.est', 'theta.est')))[,c("type", "est", "ll", "ul", "n", "j" , "seed")]
pars.ret <- subset(ests.sim.1.25.25.b.sl.c,
                   (type %in% 'pars'))[,c("n", "j" , "seed", "psi0", "theta0")]
ests.sim.1.25.25.b.sl.ret <- join(ret, pars.ret, by=c("n", "j", "seed"))

ci.correction.b <- join(subset(ests.sim.1.25.25.b.sl.ret, (type %in% 'psi.est')),
                        subset(testing.sim.1.25.25.sl, (type %in% 'Gamma.stat'))[,c("pvalue", "n", "j" , "seed")], by=c("n", "j", "seed"))

psi.b.summaries <- ddply(ci.correction.b %>% mutate(ll.c = ifelse(pvalue>0.05, 0, ll)), .(n, type), 
                         summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         coverage.c = mean(ll.c <= psi0 & psi0 <= ul, na.rm=TRUE),
                         cnt = length(type))
library(knitr)
kable(psi.b.summaries, caption="Bootstrap CI coverage for psi0 with correction")

# theta-wald
ci.correction.w <- join(subset(ests.sim.1.25.25.w.sl.ret, (type %in% 'theta.est')),
                        subset(testing.sim.1.25.25.sl, (type %in% 'Omega.stat'))[,c("pvalue", "n", "j" , "seed")], by=c("n", "j", "seed"))

theta.w.summaries <- ddply(ci.correction.w %>% mutate(ll.c = ifelse(pvalue>0.05, 0, ll)), .(n, type), 
                         summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= theta0 & theta0 <= ul, na.rm=TRUE),
                         coverage.c = mean(ll.c <= theta0 & theta0 <= ul, na.rm=TRUE),
                         cnt = length(type))
library(knitr)
kable(theta.w.summaries, caption="Wald-type CI coverage for theta0 with correction")

ci.correction.w %>% 
  mutate(ll.c = ifelse(pvalue>0.05, 0, ll),
         theta.coverage = (ll <= theta0) & (theta0 <= ul),
         theta.coverage.c = (ll.c <= theta0) & (theta0 <= ul)) %>%
  filter((!theta.coverage)&(theta.coverage.c)) %>%
  filter(n==1000) %>%
  select(c(type, est, ll, ul, ll.c, theta0, pvalue, theta.coverage, theta.coverage.c))

# theta-bootstrap
ci.correction.b <- join(subset(ests.sim.1.25.25.b.sl.ret, (type %in% 'theta.est')),
                        subset(testing.sim.1.25.25.sl, (type %in% 'Omega.stat'))[,c("pvalue", "n", "j" , "seed")], by=c("n", "j", "seed"))

theta.b.summaries <- ddply(ci.correction.b %>% mutate(ll.c = ifelse(pvalue>0.05, 0, ll)), .(n, type), 
                         summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= theta0 & theta0 <= ul, na.rm=TRUE),
                         coverage.c = mean(ll.c <= theta0 & theta0 <= ul, na.rm=TRUE),
                         cnt = length(type))
library(knitr)
kable(theta.b.summaries, caption="Bootstrap CI coverage for theta0 with correction")

ests.sim.1.25.25.b.sl.ret %>%
  subset((type %in% 'theta.est') & (n %in% c(250, 500, 750, 1000))) %>% 
  mutate(theta.coverage = (ll <= theta0) & (theta0 <= ul)) %>%
  filter(!theta.coverage) %>%
  select(c(type, est, ll, ul, theta0, theta.coverage, n)) %>%
  group_by(n) %>%
  summarise(n())

ests.sim.1.25.25.b.sl.ret %>%
  subset((type %in% 'theta.est') & (n %in% c(250, 500, 750, 1000))) %>% 
  mutate(theta.coverage = ((ll <= theta0) & (theta0 <= ul))) %>%
  filter(!theta.coverage) %>%
  ggplot()+
  geom_linerange(mapping = aes(x=j, ymin=ll, ymax=ul)) +
  # geom_point(mapping=aes(x=j, y=theta0), size=1, col="red") +
  facet_wrap(~n, scales='free') + 
  ggtitle(label="Bootstrap CIs failed to cover theta0")                                                             

ests.sim.1.25.25.b.sl.ret %>%
  subset((type %in% 'theta.est') & (n %in% c(250))) %>% 
  mutate(theta.coverage = ((ll <= theta0) & (theta0 <= ul))) %>%
  arrange(theta.coverage) %>%
  mutate(ll=ifelse(theta.coverage, -1, ll),
         ul=ifelse(theta.coverage, -1, ul),
         ranking = rank(theta.coverage, ties.method = 'first')) %>%
  # filter(!theta.coverage) %>%
  ggplot()+
  geom_linerange(mapping = aes(x=ranking, ymin=ll, ymax=ul)) +
  ylim(c(0,0.6))

ests.sim.1.25.25.b.sl.ret %>%
  subset((type %in% 'theta.est') & (n %in% c(1000))) %>% 
  mutate(theta.coverage = ((ll <= theta0) & (theta0 <= ul))) %>%
  arrange(theta.coverage) %>%
  mutate(ll=ifelse(theta.coverage, -1, ll),
         ul=ifelse(theta.coverage, -1, ul),
         ranking = rank(theta.coverage, ties.method = 'first')) %>%
  # filter(!theta.coverage) %>%
  ggplot()+
  geom_linerange(mapping = aes(x=ranking, ymin=ll, ymax=ul)) +
  ylim(c(0,0.2))
