rm(list=ls())

library(devtools)
devtools::install_github("ruihu51/CausalSim")
require(CausalSim)

load("Data/Estimator/ests.sim.2.25.75.w.sl3.c.RData")
ret <- subset(ests.sim.2.25.75.w.sl3.c,
              (type %in% c('psi.est', 'theta.est')))[,c("type", "est", "se", "ll", "ul", "n", "j" , "seed")]
optional.ret <- subset(ests.sim.2.25.75.w.sl3.c,
                       (type %in% c('psi.est', 'theta.est', 'psi.plug.in.est', 'psi.one.step.est',
                                    'theta.plug.in.est', 'theta.one.step.est')))[,c("type", "est", "n", "j" , "seed")]
pars.ret <- subset(ests.sim.2.25.75.w.sl3.c,
                   (type %in% 'pars'))[,c("n", "j" , "seed", "psi0", "theta0")]

ests.sim.2.25.75.w.sl.ret <- join(ret, pars.ret, by=c("n", "j", "seed"))
ests.sim.2.25.75.w.sl.optional.ret <- join(optional.ret, pars.ret, by=c("n", "j", "seed"))

ests.sim.2.25.75.w.sl.ret %>% tail()

seed.data <- create.seed.dict(testing.sim.2.25.75.sl3)
# generate.data.1 <- function(n){
  beta <- c(0.25, 0.75)
  pi0 <- function(w) expit(0.5+(w[,1]/3))
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  # psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  psi0 <- (beta[1]+beta[2]*(5/6))^2+(beta[2]^2*(1/5-1/9))+beta[2]^2*1/4
  # theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
  theta0 <- (beta[2]^2*(1/5-1/9))+beta[2]^2*1/4

  tau0 <- beta[1] + beta[2]*(w[,1]^2+w[,3])
  gamma0 <- beta[1] + (5/6)*beta[2]

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0, mu0=mu0, pi0=pi0)
  return(simulated.data)
}

ests <- ldply(c(100, 250, 500, 750, 1000, 1500, 2000), function(n) {
  ldply(1:1000, function(j) {

    seed <- subset(seed.data, (sample_n==n))[j,"seed"]
    cat(n, j, seed, '\n')
    set.seed(seed)

    # generate simulated data
    beta <- c(0.25, 0.75)
    pi0 <- function(w) expit(0.5+(w[,1]/3))
    mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

    W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
    A <- rbinom(n, size=1, prob=pi0(W))
    Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

    psi0 <- (beta[1]+beta[2]*(5/6))^2+(beta[2]^2*(1/5-1/9))+beta[2]^2*1/4
    theta0 <- (beta[2]^2*(1/5-1/9))+beta[2]^2*1/4

    tau0 <- beta[1] + beta[2]*(W[,1]^2+W[,3])
    gamma0 <- beta[1] + (5/6)*beta[2]

    Z <- (2*A - 1) / (A * pi0(W) + (1-A) * (1-pi0(W)))

    psi.eif <- 2 * tau0 * (Y - mu0(A, W, beta)) * Z + tau0^2 - psi0
    psi.se <- sd(psi.eif)

    theta.eif <- 2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z + (tau0-gamma0)^2 - theta0
    theta.se <- sd(theta.eif)

    psi.ret <- data.frame(type = 'psi.true', true.se = psi.se)
    theta.ret <- data.frame(type = 'theta.true', true.se = theta.se)
    ret <- rbind(psi.ret, theta.ret)
    ret <- label.result(ret, n, j, seed)

    return(ret)
  })
})

psi.empirical.sd.ret <- subset(ests.sim.2.25.75.w.sl.ret, type %in% c('psi.est')) %>%
  group_by(n) %>%
  summarise(empirical.sd=sd(sqrt(n)*(est-psi0)))

psi.sd.ret <- left_join(subset(ests.sim.2.25.75.w.sl.ret, type %in% c('psi.est')),
                        psi.empirical.sd.ret,
                        by=c("n"))

psi.sd.ret %>%
  ggplot()+
  geom_density((aes(se))) +
  geom_vline(aes(xintercept=empirical.sd),
             color="blue", linetype="dashed", size=1) +
  facet_wrap(~n, scales='free')

theta.empirical.sd.ret <- subset(ests.sim.2.25.75.w.sl.ret, type %in% c('theta.est')) %>%
  group_by(n) %>%
  summarise(empirical.sd=sd(sqrt(n)*(est-theta0)))

theta.sd.ret <- left_join(subset(ests.sim.2.25.75.w.sl.ret, type %in% c('theta.est')),
          theta.empirical.sd.ret,
          by=c("n"))

theta.sd.ret %>%
  ggplot()+
  geom_density((aes(se))) +
  geom_vline(aes(xintercept=empirical.sd),
             color="blue", linetype="dashed", size=1) +
  facet_wrap(~n, scales='free')

ggplot(theta.sd.ret) +
  geom_boxplot(aes(as.factor(n), se, color=type)) +
  # ggtitle(label=paste("convergence of psi", paste(beta, collapse = "_"))) +
  xlab("n") +
  theme(legend.key.size = unit(0.1,"line"))

theta.se.summary <- subset(ests.sim.2.25.75.w.sl.ret, type %in% c('theta.est')) %>%
  group_by(n) %>%
  summarise(cnt=n(),
            se.mean = mean(se),
            se.median = median(se))

join(theta.se.summary, theta.empirical.sd.ret,
     by=c("n"))


psi.se.summary <- subset(ests.sim.2.25.75.w.sl.ret, type %in% c('psi.est')) %>%
  group_by(n) %>%
  summarise(cnt=n(),
            se.mean = mean(se),
            se.median = median(se))

join(psi.se.summary, psi.empirical.sd.ret,
     by=c("n"))


# sd(sqrt(2000)*ests.sim.2.25.75.w.sl.ret %>%
#   filter(type %in% c('theta.est')) %>%
#   filter(n==2000) %>%
#   .[["est"]]-0.191)

psi.ret <- inner_join(subset(ests.sim.2.25.75.w.sl.ret, type %in% c('psi.est')),
                      subset(ests, type %in% c("psi.true")), by=c("n", "j", "seed"))

theta.ret <- inner_join(subset(ests.sim.2.25.75.w.sl.ret, type %in% c('theta.est')),
                      subset(ests, type %in% c("theta.true")), by=c("n", "j", "seed"))

psi.ret %>%
  ggplot()+
  geom_abline(slope=1)+
  geom_point(mapping=aes(x=true.se, y=se), size=1, col="red")+
  facet_wrap(~n, scales='free')

theta.ret %>%
  ggplot()+
  geom_abline(slope=1)+
  geom_point(mapping=aes(x=true.se, y=se), size=1, col="red")+
  xlim(c(0,4))+
  ylim(c(0,4))+
  facet_wrap(~n, scales='free')



psi.ret %>% head()

n <- 2000
j <- 1
seed <- subset(seed.data, (sample_n==n))[j,"seed"]
cat(n, j, seed, '\n')
set.seed(seed)

# generate simulated data
beta <- c(0.25, 0.75)
pi0 <- function(w) expit(0.5+(w[,1]/3))
mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
A <- rbinom(n, size=1, prob=pi0(W))
Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

psi0 <- (beta[1]+beta[2]*(5/6))^2+(beta[2]^2*(1/5-1/9))+beta[2]^2*1/4
theta0 <- (beta[2]^2*(1/5-1/9))+beta[2]^2*1/4

tau0 <- beta[1] + beta[2]*(W[,1]^2+W[,3])
gamma0 <- beta[1] + (5/6)*beta[2]

Z <- (2*A - 1) / (A * pi0(W) + (1-A) * (1-pi0(W)))

psi.eif <- 2 * tau0 * (Y - mu0(A, W, beta)) * Z + tau0^2 - psi0
psi.se <- sd(psi.eif)

theta.eif <- 2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z + (tau0-gamma0)^2 - theta0
theta.se <- sd(theta.eif)

psi.ret <- data.frame(type = 'psi.true', true.se = psi.se)
theta.ret <- data.frame(type = 'theta.true', true.se = theta.se)
ret <- rbind(psi.ret, theta.ret)
ret <- label.result(ret, n, j, seed)

seed.data %>%
  group_by(sample_n) %>%
  summarise(n_distinct(seed))

psi.summaries <- ddply(subset(ests.sim.2.25.75.w.sl.optional.ret,
                              (type %in% c('psi.est', 'psi.plug.in.est', 'psi.one.step.est'))),
                       .(n, type), summarize, na = sum(is.na(est)),
                       cnt = length(est),
                       bias = mean(est - psi0, na.rm=TRUE),
                       var = var(est, na.rm=TRUE),
                       mse = mean((est - psi0)^2, na.rm=TRUE))


htem.estimator.mu <- function(A, W, Y, control = list()){

  # update control parameters
  control <- hte.measure.NullTest.control(control)
  n = length(A)

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
    # mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
    #                        newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
    #                        SL.library = control$mu.SL.library,
    #                        family = family,
    #                        obsWeights=rep(1,n),
    #                        id=1:n)
    mu0.reg <- SuperLearner(Y=Y[which(A==0)], X = data.frame(W[which(A==0),]),
                            newX = data.frame(W),
                            SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                            family = family,
                            obsWeights=rep(1,length(which(A==0))),
                            id=1:length(which(A==0)))
    mu0.hat2 <- mu0.reg$SL.predict
    mu1.reg <- SuperLearner(Y=Y[which(A==1)], X = data.frame(W[which(A==1),]),
                            newX = data.frame(W),
                            SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                            family = family,
                            obsWeights=rep(1,length(which(A==1))),
                            id=1:length(which(A==1)))
    mu1.hat2 <- mu1.reg$SL.predict
    control$mu.hats <- data.frame(mu1=mu1.hat2, mu0=mu0.hat2)
  }
  control$mu.hat <- A * control$mu.hats$mu1 + (1-A) * control$mu.hats$mu0

  # estimated tau
  tau.hat <- control$mu.hats$mu1 - control$mu.hats$mu0
  Z.hat <- (2*A - 1) / (A * control$pi.hat + (1-A) * (1-control$pi.hat))

  # 1. estimated psi
  # a) plug-in
  psi.plug.in <- mean(tau.hat^2)
  psi.eif.hat <- 2 * tau.hat * (Y - control$mu.hat) * Z.hat + tau.hat^2 - psi.plug.in
  psi.se <- sd(psi.eif.hat)

  # b) one-step
  psi.one.step.est <- mean(2 * tau.hat * (Y - control$mu.hat) * Z.hat + tau.hat^2)

  # c) hybrid
  psi.est <- ifelse(psi.one.step.est>0, psi.one.step.est, psi.plug.in)

  # 2. estimated theta
  gamma.hat <- mean(tau.hat)

  # a) plug-in
  theta.plug.in <- mean((tau.hat-gamma.hat)^2)
  theta.eif.hat <- 2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - theta.plug.in
  theta.se <- sd(theta.eif.hat)

  # b) one-step
  theta.one.step.est <- mean(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)

  # c) hybrid
  theta.est <- ifelse(theta.one.step.est>0, theta.one.step.est, theta.plug.in)

  # estimator
  est.type.names <- names(control$est.type)
  psi.ret <- data.frame()
  theta.ret <- data.frame()
  optional.ret <- data.frame()

  # psi estimator
  if('psi.est' %in% est.type.names){
    psi.ret <- data.frame(type = 'psi.est', est = psi.est, se = psi.se)
    if (control$est.type[["psi.est"]] == 'all'){
      optional.ret <- rbind(optional.ret,
                            data.frame(type = 'psi.plug.in.est', est = psi.plug.in, se = psi.se))
      optional.ret <- rbind(optional.ret,
                            data.frame(type = 'psi.one.step.est', est = psi.one.step.est, se = psi.se))
    }
  }

  # theta estimator
  if('theta.est' %in% est.type.names){
    theta.ret <- data.frame(type = 'theta.est', est = theta.est, se = theta.se)
    if (control$est.type[["theta.est"]] == 'all'){
      optional.ret <- rbind(optional.ret,
                            data.frame(type = 'theta.plug.in.est', est = theta.plug.in, se = theta.se))
      optional.ret <- rbind(optional.ret,
                            data.frame(type = 'theta.one.step.est', est = theta.one.step.est, se = theta.se))
    }
  }
  ret <- rbind(psi.ret, theta.ret)

  # confidence interval
  if(control$conf.int){
    if(tolower(control$conf.int.type) == "wald"){
      # Wald-type CI
      ret$ll = ret$est - qnorm(1-(1-control$conf.level)/2) * ret$se / sqrt(n)
      ret$ul = ret$est + qnorm(1-(1-control$conf.level)/2) * ret$se / sqrt(n)
    } else {
      # Bootstrap CI
      boot.ests <- sapply(1:control$n.boot, function(x) {
        boot.inds <- sample(1:n, n, replace=TRUE)

        boot.pi.hat <- control$pi.hat[boot.inds]
        boot.mu.hats <- control$mu.hats[boot.inds,]
        boot.ret <- htem.estimator(A[boot.inds], W[boot.inds,], Y[boot.inds],
                                   control = list(pi.hat = boot.pi.hat,
                                                  mu.hats = boot.mu.hats,
                                                  conf.int = FALSE))
        boot.ret$ret$est
      })
      boot.ci <- apply(boot.ests, 1, quantile, c((1-control$conf.level)/2, 1-(1-control$conf.level)/2))
      ret$ll = c(boot.ci[,1][[1]], boot.ci[,2][[1]])
      ret$ul = c(boot.ci[,1][[2]], boot.ci[,2][[2]])
    }
  }

  if (length(optional.ret)>0){
    est.ret = list(ret=ret, optional.ret=optional.ret)
  }else{
    est.ret = list(ret=ret)
  }
  return(est.ret)
}
