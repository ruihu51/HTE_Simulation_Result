rm(list=ls())

library(devtools)
load_all()
devtools::install_github("ruihu51/CausalSim")

# single estimation
null.sims <- FALSE

n <- 2000
set.seed(3287481)
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

ret <- est.psi(A, W, Y, func_1 = "SL.glm", func_2 = "SL.glm")
ret
est.psi(A, W, Y, func_1 = "SL.gam", func_2 = "SL.gam")
est.psi(A, W, Y, func_1 = "SL.earth", func_2 = "SL.earth")

func_1 = "SL.gam"
n = length(A)
prop.reg1 <- do.call(func_1, list(Y=A, X = data.frame(W),
                                  newX = data.frame(W),
                                  family = binomial(),
                                  obsWeights=rep(1,n),
                                  id=1:n))
# prop.reg1 <- SL.earth(Y=A, X = data.frame(W), newX = data.frame(W), family = binomial(), obsWeights=rep(1,n), id=1:n)
pi.hat <- prop.reg1$pred
# pi.hat <- prop.reg$fitted.v

#  pi.hat mean square error
mean((pi0(W) - pi.hat)^2)

# estimated mu
func_2 = 'SL.gam'
AW <- cbind(A, data.frame(W))
if(out.glm) {
  mu.reg  <- glm(Y ~ ., data=AW, family='binomial')
  mu.hat <- mu.reg$fitted.values
  mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
  mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
} else {
  # mu.reg <- SL.earth(Y=Y, X = data.frame(cbind(A, W)), newX = rbind(data.frame(cbind(A=1, W)),data.frame(cbind(A=0, W))), family = binomial(), obsWeights=rep(1,n), id=1:n)
  mu.reg <- do.call(func_2, list(Y=Y, X = data.frame(cbind(A, W)),
                                 newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                                 family = binomial(),
                                 obsWeights=rep(1,n),
                                 id=1:n))
  mu1.hat <- mu.reg$pred[1:n]
  mu0.hat <- mu.reg$pred[-(1:n)]
  mu.hat <- A * mu1.hat + (1-A) * mu0.hat
}

#  mu.hat mean square error
mean((mu0(A, W) - mu.hat)^2)

psi0 <- mean((mu0(1,W) - mu0(0,W))^2)

# estimated tau
tau.hat <- mu1.hat - mu0.hat
Z.hat <- (2*A - 1) / (A * pi.hat + (1-A) * (1-pi.hat))

mean( 2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2)

plot(tau.hat^2 - (mu0(1,W) - mu0(0,W))^2)

# estimated theta
gamma.hat <- mean(tau.hat)

# plug-in
plug.in.theta <- mean((tau.hat-gamma.hat)^2)
se.theta <- sd( 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - plug.in.theta)
ret <- rbind(ret,
             data.frame(type = 'Plug-in (Theta)', est = plug.in.theta, ll=plug.in.theta - qnorm(.975) * se.theta / sqrt(n), ul=plug.in.theta + qnorm(.975) * se.theta / sqrt(n)))

# one-step
one.step.est.theta <- mean( 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)
ret <- rbind(ret,
             data.frame(type = 'One-step (Theta)', est = one.step.est.theta, ll=one.step.est.theta - qnorm(.975) * se.theta / sqrt(n), ul=one.step.est.theta + qnorm(.975) * se.theta / sqrt(n)))

plot(2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)



sum(subset(ests.sim.1, (type %in% c("One-step") & n==1000))$est>0)
library(dplyr)
subset(ests.sim.1, (type %in% c("One-step"))) %>%
  group_by(n) %>%
  summarise_all(~sum(est > 0))

subset(ests.sim.1, (type %in% c("One-step") & n==2000 & !(ll <= psi0 & psi0 <= ul)))
subset(ests.sim.1, (type %in% c("One-step")))
subset(summaries.psi, (type %in% c("One-step")))

ggplot(subset(ests, (type %in% c("Plug-in", "One-step", "TMLE-new")))) +
  geom_density(aes(sqrt(n) * (est - psi0), color=type)) +
  facet_wrap(~n, scales='free')

library(ggplot2)
ggplot(summaries.psi) +
  geom_line(aes(n, coverage, color=type)) +
  geom_hline(yintercept=.95)

ggplot(summaries.psi) +
  geom_line(aes(n, sqrt(n) * bias, color=type))

ggplot(summaries.psi) +
  geom_line(aes(n, n * mse, color=type))

ggplot(summaries.psi) +
  geom_line(aes(n, n * var, color=type)) +
  coord_cartesian(ylim=c(0,1))

summaries <- ddply(subset(ests, (type %in% c("Plug-in", "One-step", "TMLE-new"))), .(n, type), summarize,
                   na = sum(is.na(est)),
                   coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                   bias = mean(est - psi0, na.rm=TRUE),
                   var = var(est, na.rm=TRUE),
                   mse = mean((est - psi0)^2, na.rm=TRUE))


# testing simulation result
library(plyr)
load('testing.rst.1.h2.RData')
load('testing.rst.1.h3.RData')
load('testing.rst.1.h4.in.RData')
load('testing.rst.1.h5.in.RData')
load('testing.rst.1.h4.in2.RData')
load('testing.rst.1.h5.in.2.RData')
load('testing.rst.1.h5.in3.RData')
load('testing.rst.1.h5.in4.RData')
load('testing.rst.1.h5.in5.RData')
load('testing.rst.1.h5.in6.RData')
load('testing.rst.1.h5.in7.RData')

testing.rst.1.h4 <- rbind(testing.rst.1.h4.in, testing.rst.1.h4.in2)
save(testing.rst.1.h4, file = 'testing.rst.1.h4.Rdata')

testing.rst.1.h5 <- rbind(rbind(rbind(rbind(rbind(rbind(testing.rst.1.h5.in, testing.rst.1.h5.in.2),
                                testing.rst.1.h5.in3),
                          testing.rst.1.h5.in4),
                          testing.rst.1.h5.in5),testing.rst.1.h5.in6),testing.rst.1.h5.in7)
save(testing.rst.1.h5, file = 'testing.rst.1.h5.Rdata')

length(unique(testing.rst.1.h5$seed))
testing.rst.1.h5$j

testing.rst.1.h5 <- testing.rst.1.h5.in3

testing.rst.1 <- rbind(rbind(rbind(testing.rst.1.h2, testing.rst.1.h3),testing.rst.1.h4), testing.rst.1.h5)
testing.rst.1
save(testing.rst.1, file = 'testing.rst.1.Rdata')

gamma.summaries.1 <- ddply(subset(testing.rst.1, (type %in% 'Gamma.stat')), .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))

omega.summaries.1 <- ddply(subset(testing.rst.1, (type %in% 'Omega.stat')), .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))

load('testing.rst.2.h1.in.RData')
load('testing.rst.2.h2.RData')
load('testing.rst.2.h3.RData')
load('testing.rst.2.h4.in1.RData')

testing.rst.2 <- rbind(rbind(rbind(testing.rst.2.h1.in, testing.rst.2.h2),
                       testing.rst.2.h3),testing.rst.2.h4.in1)
testing.rst.2

gamma.summaries.2 <- ddply(subset(testing.rst.2, (type %in% 'Gamma.stat')), .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))

omega.summaries.2 <- ddply(subset(testing.rst.2, (type %in% 'Omega.stat')), .(n, type), summarize,
                           na = sum(is.na(stat)),
                           cnt = length(stat),
                           quantile.reject.rate = mean(stat>quantile),
                           pvalue.reject.rate = mean(pvalue<0.05))






