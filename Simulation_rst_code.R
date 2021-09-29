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
save(ests.sim.1, file = "ests.sim.1.RData")

## Bootstrap
ests.sim.2 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)
save(ests.sim.2, file = "ests.sim.2.RData")

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
