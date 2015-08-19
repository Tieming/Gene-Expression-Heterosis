##
## Multiple Simulations to Collect Data for Tables and Plots.
##
## 12/7/11. Tieming
##
##


setwd("/Users/Tieming/research/proj2/manuscript-10-9-11-simulation-prob-models/Approximate-Posterior")


##
## Simulation starts here.
##

install.packages(c("MCMCpack","coda","lattice","MASS","geoR"))
library("lattice")
library("MASS")
library("coda")
library("MCMCpack")
library("geoR")



start <- Sys.time()
sim.total <- 100
J0 <- 5000
sim.para <- c(0.8, 2, 0.01, 0.6, -0.05, 0.2, 0, 0.2, 2.8, 0.025)
names(sim.para) <- c("p0","t.dg","t.shift","pi0","n1.mean","n1.s2","n2.mean","n2.s2","d0","s20")
## Bias and MSE of Heterosis for SA and EB methods.
hph.sa.bias <- 0
lph.sa.bias <- 0
mph.sa.bias <- 0
hph.eb.bias <- 0
lph.eb.bias <- 0
mph.eb.bias <- 0
hph.sa.mse  <- 0
lph.sa.mse  <- 0
mph.sa.mse  <- 0
hph.eb.mse  <- 0
lph.eb.mse  <- 0
mph.eb.mse  <- 0
## Average bias across genes for MPH.
mph.sa.per.sim.bias <- NULL
mph.eb.per.sim.bias <- NULL
## FPr vs. TPr (ROC) curves.
hph.roc <- NULL
lph.roc <- NULL
mph.roc <- NULL
## Posterior Probability of Exhibiting Heterosis.
estd.hph    <- NULL
estd.mph    <- NULL
estd.lph    <- NULL
simd.gene.hph   <- NULL
simd.gene.mph   <- NULL
simd.gene.lph   <- NULL



seed <- 51147602
for (sim.time in 1:sim.total){
  cat("simulation ", sim.time, fill=TRUE)
  seed <- seed + sim.time
  
  ## simulate. 
  source("/Users/Tieming/research/proj2/manuscript-10-9-11-simulation-prob-models/Approximate-Posterior/sim_function.R")
  simd  <- sim(para=sim.para, I=3, J=J0, seed=seed)

  ## estimate.
  source("/Users/Tieming/research/proj2/manuscript-10-9-11-simulation-prob-models/Approximate-Posterior/empirical_estimate_function.R")
  estd  <- estimate(simd, alpha.obs=simd$alpha.hat, 
                      delta.obs=simd$delta.hat,
                      s2.ej.obs=simd$s2.ej.obs,
                      I=3, J=J0, mcn=1000)


   ## Table for Bias and MSE of SA and EB methods.
   hph.sa.bias <- (hph.sa.bias*(sim.time-1) + mean(simd$delta.hat-abs(simd$alpha.hat)-simd$hph))/sim.time
   lph.sa.bias <- (lph.sa.bias*(sim.time-1) + mean(-abs(simd$alpha.hat)-simd$delta.hat-simd$lph))/sim.time
   mph.sa.bias <- (mph.sa.bias*(sim.time-1) + mean(simd$delta.hat-simd$mph))/sim.time
   hph.eb.bias <- (hph.eb.bias*(sim.time-1) + mean(estd$delta.pos-estd$abs.alpha.pos-simd$hph))/sim.time
   lph.eb.bias <- (lph.eb.bias*(sim.time-1) + mean(-estd$abs.alpha.pos-estd$delta.pos-simd$lph))/sim.time
   mph.eb.bias <- (mph.eb.bias*(sim.time-1) + mean(estd$delta.pos-simd$mph))/sim.time
   hph.sa.mse  <- (hph.sa.mse*(sim.time-1) + mean((simd$delta.hat-abs(simd$alpha.hat)-simd$hph)^2))/sim.time
   lph.sa.mse  <- (lph.sa.mse*(sim.time-1) + mean((-abs(simd$alpha.hat)-simd$delta.hat-simd$lph)^2))/sim.time
   mph.sa.mse  <- (mph.sa.mse*(sim.time-1) + mean((simd$delta.hat-simd$mph)^2))/sim.time
   hph.eb.mse  <- (hph.eb.mse*(sim.time-1) + mean((estd$delta.pos-estd$abs.alpha.pos-simd$hph)^2))/sim.time
   lph.eb.mse  <- (lph.eb.mse*(sim.time-1) + mean((-estd$abs.alpha.pos-estd$delta.pos-simd$lph)^2))/sim.time
   mph.eb.mse  <- (mph.eb.mse*(sim.time-1) + mean((estd$delta.pos-simd$mph)^2))/sim.time
 
   ## Average bias across genes of MPH for each simulation.
   mph.sa.per.sim.bias <- c(mph.sa.per.sim.bias, mean(simd$delta.hat-simd$mph))
   mph.eb.per.sim.bias <- c(mph.eb.per.sim.bias, mean(estd$delta.pos-simd$mph))
       
   ## Average of Ranked Estimation Error.
   if(sim.time==1){
   	  hph.rank.err.sa <- sort(simd$delta.hat-abs(simd$alpha.hat)-simd$hph)
   	  lph.rank.err.sa <- sort(-abs(simd$alpha.hat)-simd$delta.hat-simd$lph)
   	  mph.rank.err.sa <- sort(simd$delta.hat-simd$mph)
   	  hph.rank.err.eb <- sort(estd$delta.pos-estd$abs.alpha.pos-simd$hph)
   	  lph.rank.err.eb <- sort(-estd$abs.alpha.pos-estd$delta.pos-simd$lph)
   	  mph.rank.err.eb <- sort(estd$delta.pos-simd$mph)
   	}
   hph.rank.err.sa <- (hph.rank.err.sa*(sim.time-1)+sort(simd$delta.hat-abs(simd$alpha.hat)-simd$hph))/sim.time
   lph.rank.err.sa <- (lph.rank.err.sa*(sim.time-1)+sort(-abs(simd$alpha.hat)-simd$delta.hat-simd$lph))/sim.time
   mph.rank.err.sa <- (mph.rank.err.sa*(sim.time-1)+sort(simd$delta.hat-simd$mph))/sim.time
   hph.rank.err.eb <- (hph.rank.err.eb*(sim.time-1)+sort(estd$delta.pos-estd$abs.alpha.pos-simd$hph))/sim.time
   lph.rank.err.eb <- (lph.rank.err.eb*(sim.time-1)+sort(-estd$abs.alpha.pos-estd$delta.pos-simd$lph))/sim.time
   mph.rank.err.eb <- (mph.rank.err.eb*(sim.time-1)+sort(estd$delta.pos-simd$mph))/sim.time
 
 
   ## ROC Plots:
   if(sim.time ==1){
   	  hph.roc <- estd$hph.roc
   	  lph.roc <- estd$lph.roc
   	  mph.roc <- estd$mph.roc
   	}
   hph.roc <- cbind(hph.roc, estd$hph.roc)
   lph.roc <- cbind(lph.roc, estd$lph.roc)
   mph.roc <- cbind(mph.roc, estd$mph.roc)
 
   ## (FDR) vs. (Estimated FDR) Plots:
   ## HPH:
   estd.hph <- cbind(estd.hph, estd$hph.p)
   simd.gene.hph <- cbind(simd.gene.hph, simd$hph.gene)
   ## LPH:
   estd.lph <- cbind(estd.lph, estd$lph.p)
   simd.gene.lph <- cbind(simd.gene.lph, simd$lph.gene)
   ## MPH:
   estd.mph <- cbind(estd.mph, estd$mph.p)
   simd.gene.mph <- cbind(simd.gene.mph, simd$mph.gene)
}



Sys.time() - start
 


### Bias and MSE for SA and EB methods.
###
###
hph.sa.bias
hph.eb.bias
lph.sa.bias
lph.eb.bias
mph.sa.bias
mph.eb.bias
hph.sa.mse
hph.eb.mse
lph.sa.mse
lph.eb.mse
mph.sa.mse
mph.eb.mse



### Posterior Means for |alpha|, MPH, HPH, LPH -- Box Plots:
###
### (1) HPH:
pdf(file="prob-hph-box.pdf")
par(mar=c(5,5,1,1))
boxplot(hph.rank.err.sa, hph.rank.err.eb, names=c("Sample Average","Empirical Bayes"), ylab="Average of Ranked Estimation Error", cex.axis=1.7, cex.lab=2, cex=1.2, lwd=2)
abline(h=0, lty="dashed", col="darkgray", lwd=3)
dev.off()
### (2) LPH:
pdf(file="prob-lph-box.pdf")
par(mar=c(5,5,1,1))
boxplot(lph.rank.err.sa, lph.rank.err.eb, names=c("Sample Average","Empirical Bayes"), ylab="Average of Ranked Estimation Error", cex.axis=1.7, cex.lab=2, cex=1.2, lwd=2)
abline(h=0, lty="dashed", col="darkgray", lwd=3)
dev.off()
### (3) MPH:
pdf(file="prob-mph-box.pdf")
par(mar=c(5,5,1,1))
boxplot(mph.rank.err.sa, mph.rank.err.eb, names=c("Sample Average","Empirical Bayes"), ylab="Average of Ranked Estimation Error", cex.axis=1.7, cex.lab=2, cex=1.2, lwd=2)
abline(h=0, lty="dashed", col="darkgray", lwd=3)
dev.off()



### ROC curve (FPr vs. TPr):
###
### (1) HPH:
hph.o <- roc.fcn(hph.roc, dg=3, sim.total)
pdf(file="prob-hph-roc.pdf")
par(mar=c(5,5,1,1))
plot(hph.o$fpr, hph.o$eb.tpr, ylab="TPr", xlab="FPr", cex.lab=1.7, cex.axis=1.7, lwd=3, type="l", col="red", ylim=range(c(hph.o$eb.tpr,hph.o$sa.tpr)))
lines(hph.o$fpr, hph.o$sa.tpr, lwd=3)
legend("bottomright", c("Sample Average","Empirical Bayes"), col=c("black","red"), lwd=3, bty="n", cex=1.7)
dev.off()
### (2) LPH:
lph.o <- roc.fcn(lph.roc, dg=3, sim.total)
pdf(file="prob-lph-roc.pdf")
par(mar=c(5,5,1,1))
plot(lph.o$fpr, lph.o$eb.tpr, ylab="TPr", xlab="FPr", cex.lab=1.7, cex.axis=1.7, lwd=3, type="l", col="red", ylim=range(c(lph.o$eb.tpr,lph.o$sa.tpr)))
lines(lph.o$fpr, lph.o$sa.tpr, lwd=3)
legend("bottomright", c("Sample Average","Empirical Bayes"), col=c("black","red"), lwd=3, bty="n", cex=1.7)
dev.off()
### (3) MPH:
mph.o <- roc.fcn(mph.roc, dg=3, sim.total)
pdf(file="prob-mph-roc.pdf")
par(mar=c(5,5,1,1))
plot(mph.o$fpr, mph.o$eb.tpr, ylab="TPr", xlab="FPr", cex.lab=1.7, cex.axis=1.7, lwd=3, type="l", col="red", ylim=range(c(mph.o$eb.tpr,mph.o$sa.tpr)))
lines(mph.o$fpr, mph.o$sa.tpr, lwd=3)
legend("bottomright", c("Sample Average","Empirical Bayes"), col=c("black","red"), lwd=3, bty="n", cex=1.7)
dev.off()




roc.fcn <- function(tpr.fpr, dg=3, sim.total){
## This function computes the average TPr given FPr for both 
## Sample Average Method and Empirical Bayes Method.
##
  tpr.fpr <- round(tpr.fpr, digits=dg)
  step    <- 1/(10^dg)
  sa.tpr  <- NULL
  eb.tpr  <- NULL
  sa.tpr.empty.index <- NULL
  eb.tpr.empty.index <- NULL
  sa.fpr.matrix  <- tpr.fpr[,seq(2, dim(tpr.fpr)[2], 4)]
  eb.fpr.matrix  <- tpr.fpr[,seq(4, dim(tpr.fpr)[2], 4)] 
  index <- 0
  for(step.val in seq(0, 0.05, step)){
  	index <- index + 1
  	next.sa <- NULL
  	next.eb <- NULL
   for(k in 1:sim.total){
      next.place.sa <- sum(sa.fpr.matrix[,k]<=step.val)
      next.sa <- c(next.sa, tpr.fpr[next.place.sa,1+4*(k-1)])
      next.place.eb <- sum(eb.fpr.matrix[,k]<=step.val)
      next.eb <- c(next.eb, tpr.fpr[next.place.eb,3+4*(k-1)])
   }
   next.sa   <- mean(next.sa)
   sa.tpr    <- c(sa.tpr, next.sa)
   next.eb   <- mean(next.eb)
   eb.tpr    <- c(eb.tpr, next.eb)
  }   
  fpr <- seq(0, 0.05, step)
  return(list(sa.tpr=sa.tpr, eb.tpr=eb.tpr, fpr=fpr))
}







### (FDR) vs. (Estiamted FDR) Plots:
###
### (1) HPH:
each.fdr <- fdr.one.exp(estd.hph, simd.gene.hph, sim.total, J=J0)
avg.fdr  <- fdr.avg.exp.update(each.fdr$true.fdr, each.fdr$est.fdr, sim.total, step=0.001, dg=3)
### (2) LPH:
each.fdr <- fdr.one.exp(estd.lph, simd.gene.lph, sim.total, J=J0)
avg.fdr  <- fdr.avg.exp.update(each.fdr$true.fdr, each.fdr$est.fdr, sim.total, step=0.001, dg=3)
### (3) MPH:
each.fdr <- fdr.one.exp(estd.mph, simd.gene.mph, sim.total, J=J0)
avg.fdr  <- fdr.avg.exp.update(each.fdr$true.fdr, each.fdr$est.fdr, sim.total, step=0.001, dg=3)

pdf(file="prob-mph-fdr.pdf")
par(mar=c(5,5,1,1))
plot(avg.fdr$true.fdr, avg.fdr$est.fdr, xlab="TRUE FDR", ylab="ESTIMATED FDR", type="l", lwd=2, xlim=c(0,0.25), ylim=c(0,0.25), cex.lab=1.7, cex.axis=1.7)
abline(a=0,b=1, lty="dashed", col="darkgray", lwd=3)
lines(avg.fdr$true.fdr, avg.fdr$est.fdr, lwd=3)
dev.off()


fdr.one.exp <- function(ppde, true.positive, sim.total, J){
# calculate FDR and estimated FDR based on PPDE 
# for each experiment/simulation.
#
  for(sim in 1:sim.total){
    ord <- order(ppde[,sim], decreasing=TRUE)
    ppde[,sim] <- ppde[ord,sim]
    true.positive[,sim] <- true.positive[ord,sim]
  }
  true.fdr <- apply(1-true.positive, 2, cumsum)/seq(1,J,1)
  colnames(true.fdr) <- seq(1,sim.total,1)
  est.fdr  <- apply(1-ppde, 2, cumsum)/seq(1,J,1)
  colnames(est.fdr) <- seq(1,sim.total,1)
  return(list(true.fdr=true.fdr, est.fdr=est.fdr))
}


fdr.avg.exp <- function(true.fdr, est.fdr, sim.total, step, dg){
# Based on the true FDR and estimated FDR for each experiment/simulation,
# calculate the average FDR across experiments/simulations.
#
  true.fdr <- signif(true.fdr, digits=dg)
  est.fdr  <- signif(est.fdr, digits=dg)
  avg.true.fdr <- NULL
  for(sim in 1:sim.total){
  	tmp.true.fdr <- NULL
    for(inv.step in seq(0,0.25,step)){
      tmp.place <- sum(est.fdr[,sim]<=inv.step)
      if(tmp.place==0){
        tmp.true.fdr <- c(tmp.true.fdr, 0)
      }else{
        tmp.true.fdr <- c(tmp.true.fdr, true.fdr[tmp.place,sim])
      }
    }
    avg.true.fdr <- cbind(avg.true.fdr, tmp.true.fdr)
  }
  avg.true.fdr <- apply(avg.true.fdr, 1, mean)
  avg.est.fdr  <- seq(0,0.25,step)
  return(list(true.fdr=avg.true.fdr, est.fdr=avg.est.fdr))
}


fdr.avg.exp.update <- function(true.fdr, est.fdr, sim.total, step, dg){
# Based on the true FDR and estimated FDR for each experiment/simulation,
# calculate the average FDR across experiments/simulations.
#
  true.fdr <- signif(true.fdr, digits=dg)
  est.fdr  <- signif(est.fdr, digits=dg)
  max.true.fdr <- min(max(true.fdr), 0.25)
  avg.true.fdr <- NULL
  for(sim in 1:sim.total){
  	tmp.true.fdr <- NULL
    for(inv.step in seq(0,max.true.fdr,step)){
      tmp.place <- sum(est.fdr[,sim]<=inv.step)
      if(tmp.place==0){
        tmp.true.fdr <- c(tmp.true.fdr, 0)
      }else{
        tmp.true.fdr <- c(tmp.true.fdr, true.fdr[tmp.place,sim])
      }
    }
    avg.true.fdr <- cbind(avg.true.fdr, tmp.true.fdr)
  }
  avg.true.fdr <- apply(avg.true.fdr, 1, mean)
  avg.est.fdr  <- seq(0,max.true.fdr,step)
  return(list(true.fdr=avg.true.fdr, est.fdr=avg.est.fdr))
}



###
### compute probability conclusions in the 'real data analysis' section
###

# simulate alpha according to the estimated hyper-parameters.
J <- 52973
u <- runif(n=J)
alpha.sim <- rep(0, J) 
alpha.sim[u>alfalfa.estimate$para["p0"]] <- rnorm(sum(u>alfalfa.estimate$para["p0"]),
                           mean=alfalfa.estimate$para["alpha.mu"], 
                           sd=sqrt(alfalfa.estimate$para["alpha.s2"]))

# simulate delta according to the estimated hyper-parameters.
delta.sim <- rep(0, J)
u <- runif(n=J)
delta.sim[u>alfalfa.estimate$para["pi0"]] <- rnorm(sum(u>alfalfa.estimate$para["pi0"]),
                           mean=alfalfa.estimate$para["delta.mu"],
                           sd=sqrt(alfalfa.estimate$para["delta.s2"]))

sum(delta.sim > abs(alpha.sim))/J
sum(delta.sim < -abs(alpha.sim))/J







 




