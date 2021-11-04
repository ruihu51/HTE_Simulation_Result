rm(list=ls())

# library(plyr)
# install.packages('gam')
library(gam)

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

# simulation for paper
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

control <- list()
out.glm <- FALSE
ret <- data.frame()
# n=100
for (n in 100){
  for (j in 1:500){
    # if(j %% 100 == 0)
    cat(n, j, '\n')
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)

    # generate simulated data
    simulated.data <- do.call(generate.data, list(n=n))
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
testing.sim.1.25.0.sl.100 <- ret
save(testing.sim.1.25.0.sl.100, file = "testing.sim.1.25.0.sl.100.RData")

control <- list()
out.glm <- FALSE
ret <- data.frame()
# n=250
for (n in 250){
  for (j in 1:500){
    # if(j %% 100 == 0)
    cat(n, j, '\n')
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)

    # generate simulated data
    simulated.data <- do.call(generate.data, list(n=n))
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
testing.sim.1.25.0.sl.250 <- ret
save(testing.sim.1.25.0.sl.250, file = "testing.sim.1.25.0.sl.250.RData")

control <- list()
out.glm <- FALSE
ret <- data.frame()
# n=500
for (n in 500){
  for (j in 1:500){
    # if(j %% 100 == 0)
    cat(n, j, '\n')
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)

    # generate simulated data
    simulated.data <- do.call(generate.data, list(n=n))
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
testing.sim.1.25.0.sl.500 <- ret
save(testing.sim.1.25.0.sl.500, file = "testing.sim.1.25.0.sl.500.RData")

control <- list()
out.glm <- FALSE
ret <- data.frame()
# n=750
for (n in 750){
  for (j in 1:500){
    # if(j %% 100 == 0)
    cat(n, j, '\n')
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)

    # generate simulated data
    simulated.data <- do.call(generate.data, list(n=n))
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
testing.sim.1.25.0.sl.750 <- ret
save(testing.sim.1.25.0.sl.750, file = "testing.sim.1.25.0.sl.750.RData")

control <- list()
out.glm <- FALSE
ret <- data.frame()
# n=1000
for (n in 1000){
  for (j in 1:500){
    # if(j %% 100 == 0)
    cat(n, j, '\n')
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)

    # generate simulated data
    simulated.data <- do.call(generate.data, list(n=n))
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
testing.sim.1.25.0.sl.1000 <- ret
save(testing.sim.1.25.0.sl.1000, file = "testing.sim.1.25.0.sl.1000.RData")
