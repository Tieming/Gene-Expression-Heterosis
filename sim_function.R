#
#
#  Simulation
#
#  sim function.
#
#  11/17/10. Tieming
#
#


sim <- function(para, I, J, seed){
# This is the function to generate data.
#
  set.seed(seed)

  ########
  # Simulate true values of alpha, delta and s2.ej from the
  # parent distribution.
  ########
  alpha       <- rep(0, J)
  samp        <- sample(1:J, ceiling((1-para["p0"])*J))
  alpha[samp] <- rt(length(samp), df=para["t.dg"], ncp=para["t.shift"])
  
  delta        <- rep(0, J)
  samp         <- sample(1:J, ceiling((1-para["pi0"])*J))
  random.index <- sample(1:length(samp), size=length(samp)/2)
  delta[samp][random.index]  <- rnorm(length(samp)/2, mean=para["n1.mean"], sd=sqrt(para["n1.s2"]))
  delta[samp][-random.index] <- rnorm(length(samp)/2, mean=para["n2.mean"], sd=sqrt(para["n2.s2"]))
  
  s2.ej     <- para["d0"]*para["s20"]*rinvchisq(J, df=para["d0"])
  df        <- 3*I-3
  s2.ej.obs <- s2.ej*rchisq(J, df)/df
  
  # Simulate alpha.hat and delta.hat based on alpha and delta.
  s2.j      <- s2.ej*2/I
  alpha.hat <- rnorm(J, alpha, sqrt(s2.j/4))
  delta.hat <- rnorm(J, delta, sqrt(s2.j*3/4))
  
  # True heterosis.
  hph <- delta - abs(alpha)
  hph.gene <- rep(0, J)
  hph.gene[delta > abs(alpha)] <- 1
  
  lph <- -abs(alpha) - delta
  lph.gene <- rep(0, J)
  lph.gene[delta < -abs(alpha)] <- 1
  
  mph <- delta
  mph.gene <- rep(0, J)
  mph.gene[delta != 0] <- 1
    
  return(list(alpha=alpha, delta=delta, 
              alpha.hat=alpha.hat, delta.hat=delta.hat,
              s2.ej=s2.ej,
              s2.ej.obs=s2.ej.obs,
              hph=hph, hph.gene=hph.gene,
              lph=lph, lph.gene=lph.gene,
              mph=mph, mph.gene=mph.gene))
}
  














