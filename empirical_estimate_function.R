#
#
#  estimate function.
#
#  12/7/10. Tieming
#
#




estimate <- function(simd, alpha.obs, delta.obs, s2.ej.obs, I, J, mcn){
## This is the function to analyze the simualted data by two methods.
## d: a list containing (1) data; (2) t.diff; (3) true.max.
## p: proportion of the 0 treatment difference in mixture model.
## I: number of replicates in CRD or number of blocks in RCBD.
## J: total number of genes.
## exp.design: design of the experiment, either CRD or RCBD.
##
  
  ## Parameter estimation:
  ##
  ## (1) Estimate d0 and s20:
  df <- I*3-3
  limma.para <- limma.para.f(log(s2.ej.obs), df)
  d0  <- limma.para[1]
  s20 <- limma.para[2]
  ## (1) Other hyper-parameters:
  pos.expected.s2.ej <- (df*s2.ej.obs + d0*s20)/(df+d0-2)
  pos.expected.s2.j  <- pos.expected.s2.ej*2/I
  para.est.1 <- p.solution(alpha.obs, pos.expected.s2.j/4)
  para.lik.1 <- p.likhood.solution(alpha.obs, pos.expected.s2.j/4,
  										start=c(para.est.1[c("p","mu","s2")]))
  para.est.2 <- p.solution(delta.obs, pos.expected.s2.j*3/4)
  para.lik.2 <- p.likhood.solution(delta.obs, pos.expected.s2.j*3/4,
  										start=c(para.est.2[c("p","mu","s2")]))
  para <- c(para.lik.1, para.lik.2)
  names(para) <- c("p0","alpha.mu","alpha.s2","pi0","delta.mu","delta.s2")
         
                  
  ## Approximate the posterior coefficients for the joint mixture distribution
  ## of (alpha, delta | alpha.obs, delta.obs, sample error varince) in 
  ## four parts for all genes.
  ##    
  c1 <- 0
  c2 <- 0
  c3 <- 0
  c4 <- 0
  for (mcn.cnt in 1:mcn){
    v.e.rand  <- 1/rgamma(J, shape=(df+d0)/2, rate=(df*s2.ej.obs+d0*s20)/2)
    #v.e.rand   <- rinvgamma(J, shape=(df+d0)/2, scale=(df*s2.ej.obs+d0*s20)/2)
    s2.j.rand <- v.e.rand*2/I
    
    c1 <- c1 + para["p0"]*para["pi0"]*
          dnorm(alpha.obs, 0, sqrt(s2.j.rand/4))*
          dnorm(delta.obs, 0, sqrt(s2.j.rand*3/4))
    c2 <- c2 + (1-para["p0"])*para["pi0"]*
          dnorm(alpha.obs, para["alpha.mu"], sqrt(para["alpha.s2"]+s2.j.rand/4))*
          dnorm(delta.obs, 0, sqrt(s2.j.rand*3/4))
    c3 <- c3 + para["p0"]*(1-para["pi0"])*
          dnorm(alpha.obs, 0, sqrt(s2.j.rand/4))*
          dnorm(delta.obs, para["delta.mu"], sqrt(para["delta.s2"]+s2.j.rand*3/4))
    c4 <- c4 + (1-para["p0"])*(1-para["pi0"])*
          dnorm(alpha.obs, para["alpha.mu"], sqrt(para["alpha.s2"]+s2.j.rand/4))*
          dnorm(delta.obs, para["delta.mu"], sqrt(para["delta.s2"]+s2.j.rand*3/4))
  }
  c <- c1 + c2 + c3 + c4
  c1 <- c1/c
  c2 <- c2/c
  c3 <- c3/c
  c4 <- c4/c
  
  
  ## Randomly draw samples from posterior joint mixture distribution
  ## for alpha and delta. 
  ##
  
  ## (1) Initial values for the posterior probability of heterosis. 
  hph.cnt     <- rep(0, J)
  lph.cnt     <- rep(0, J)
  mph.cnt     <- rep(0, J)
  ## (2) Initial values for the posterior mean of |alpha|, hph, lph, mph.
  abs.alpha.pos <- rep(0, J)
  delta.pos     <- rep(0, J)
  #hph.pos       <- rep(0, J)
  #lph.pos       <- rep(0, J)
 
  ## The "for" loop for random draws from the posterior distribution.
  for(mcn.cnt in 1:mcn){
  	alpha.draw <- rep(0,J)
  	delta.draw <- rep(0,J)
    ru <- runif(J)
    alpha.nonzero.index <- ((ru>c1 & ru<=(c1+c2)) | (ru>(c1+c2+c3)))
    alpha.draw[alpha.nonzero.index] <- 
           rnorm(sum(alpha.nonzero.index),
           (alpha.obs[alpha.nonzero.index]/(pos.expected.s2.j[alpha.nonzero.index]/4)+
           para["alpha.mu"]/para["alpha.s2"])/
           (1/(pos.expected.s2.j[alpha.nonzero.index]/4)+1/para["alpha.s2"]),
           sqrt(1/(1/(pos.expected.s2.j[alpha.nonzero.index]/4)+
           1/para["alpha.s2"])))
    delta.nonzero.index <- (ru>(c1+c2))
    delta.draw[delta.nonzero.index] <-
           rnorm(sum(delta.nonzero.index),
           (delta.obs[delta.nonzero.index]/(pos.expected.s2.j[delta.nonzero.index]*3/4)+
           para["delta.mu"]/para["delta.s2"])/
           (1/(pos.expected.s2.j[delta.nonzero.index]*3/4)+1/para["delta.s2"]),
           sqrt(1/(1/(pos.expected.s2.j[delta.nonzero.index]*3/4)+
           1/para["delta.s2"])))
    ## (1) Estimate posterior probability of heterosis.
    hph.cnt <- hph.cnt + as.numeric(delta.draw > abs(alpha.draw))
    lph.cnt <- lph.cnt + as.numeric(delta.draw < -abs(alpha.draw))
    mph.cnt <- mph.cnt + as.numeric(delta.draw != 0)    
    ## (2) Estimate for posterior mean of |alpha|, hph, lph, mph.
    abs.alpha.pos <- (abs.alpha.pos*(mcn.cnt-1) + abs(alpha.draw))/mcn.cnt
    delta.pos     <- (delta.pos*(mcn.cnt-1) + delta.draw)/mcn.cnt
  }
  ## (1) (Continued) Estimate posterior probability of heterosis.
  hph.p <- hph.cnt/mcn
  lph.p <- lph.cnt/mcn
  mph.p <- mph.cnt/mcn
  
  
  
  
  ## Compute FPr and TPr for ROC for this (one) simulation.
  ##
  ## (1) Sample average method. 
  ##     MSE (Mean Square Error) is proportional to 
  ##     error variance (s2.ej.obs).
  t.mph.mse <- delta.obs/sqrt(s2.ej.obs)
  t.hph.mse <- (delta.obs-abs(alpha.obs))/sqrt(s2.ej.obs)
  t.lph.mse <- (-abs(alpha.obs)-delta.obs)/sqrt(s2.ej.obs)
  ## (2) Empirical Bayes method.
  ##     Use hph.p, lph.p and mph.p.
  ## (3) Ranking of genes.
  ## (3.1) HPH:
  sig.hph     <- simd$hph.gene
  ord.hph.mse <- order(t.hph.mse, decreasing=TRUE)
  ord.hph.eb  <- order(hph.p, decreasing=TRUE)
  hph         <- cbind(sig.hph[ord.hph.mse], sig.hph[ord.hph.eb])
  hph.roc     <- gen.roc(hph, m1=sum(simd$hph.gene), m0=J-sum(simd$hph.gene))
  ## (3.2) LPH:
  sig.lph 	   <- simd$lph.gene
  ord.lph.mse <- order(t.lph.mse, decreasing=TRUE)
  ord.lph.eb  <- order(lph.p, decreasing=TRUE)
  lph         <- cbind(sig.lph[ord.lph.mse], sig.lph[ord.lph.eb])
  lph.roc     <- gen.roc(lph, m1=sum(simd$lph.gene), m0=J-sum(simd$lph.gene))
  ## (3.3) MPH:
  sig.mph     <- simd$mph.gene
  ord.mph.mse <- order(abs(t.mph.mse), decreasing=TRUE)
  ord.mph.eb  <- order(mph.p, decreasing=TRUE)
  mph         <- cbind(sig.mph[ord.mph.mse], sig.mph[ord.mph.eb])
  mph.roc     <- gen.roc(mph, m1=sum(simd$mph.gene), m0=J-sum(simd$mph.gene))
  
  return(list(para=para,
              hph.p=hph.p, lph.p=lph.p, mph.p=mph.p,
              abs.alpha.pos=abs.alpha.pos, delta.pos=delta.pos,
              hph.roc=hph.roc, lph.roc=lph.roc, mph.roc=mph.roc))
}




limma.para.f <- function(z, df){
# Given the REML estimates of residual variances, compute the
# hyper-parameters d0 and s20 by limma method.
#
  y <- var(z)-trigamma(df/2)
  x <- seq(0.01,df,0.01)
  hat.y <- trigamma(x)
  d0 <- 2*x[order(abs(y-hat.y))[1]]
  s20 <- exp(mean(z)-digamma(df/2)+digamma(d0/2)-log(d0/df))
  return(c(d0,s20))
}


p.solution <- function(y, s2.j){
# This is the function of estimate the optimal solution for
# parameters to achieve the maximum likelihood while giving
# the constrain that p is between 0 and 1.
#
  m1 <- mean(y)
  m2 <- mean(y^2)
 
  fdensity <- function(p){
  	s <- max((((m2 - mean(s2.j))/(1-p)) - (m1/(1-p))^2),0)
  	
  	sum(log(p*dnorm(y, 0, sqrt(s2.j)) + 
       (1-p)*dnorm(y, m1/(1-p), sqrt(s+s2.j))))
  }
   
  psol <- optimize(fdensity, interval=c(0,1), maximum=TRUE)$maximum
  mu.est <- min(m1/(1-psol), max(y))
  mu.est <- max(mu.est, min(y))
  s2.est <- max((m2 - mean(s2.j))/(1-psol) - mu.est^2,0.00001)
  return(c(p=psol,mu=mu.est,s2=s2.est))
}




p.likhood.solution <- function(y, s2.j, start){
# This is the function of estimate the optimal solution for
# parameters to achieve the maximum likelihood while giving
# the constrain that p is between 0 and 1.
#
  fdensity <- function(para){
    p  <- para[1]
    mu <- para[2]
    s2 <- para[3]
    
    negloglik <- -sum(log(p*dnorm(y, 0, sqrt(s2.j))+
                          (1-p)*dnorm(y, mu, sqrt(s2+s2.j))))
    return(negloglik)  
  }  
  
  para <- nlminb(start=start, objective=fdensity,
          lower=c(0,min(y),0), upper=c(1,max(y),Inf))$par 
  names(para) <- c("p","mu","s2")
  return(para)
}


gen.roc <- function(gt, m1, m0){
# This function computes the TPr and FPr given
# the vector of "gene" sorted by either the p values
# or the t values.
  tpr <- apply(gt, 2, cumsum)/m1
  fpr <- apply(abs(gt-1), 2, cumsum)/m0
  return(data.frame(avg.tpr=tpr[,1], avg.fpr=fpr[,1], 
          eb.tpr=tpr[,2], eb.fpr=fpr[,2]))
}








