rm(list=ls())

library(MASS)
library(GIGrvg)

source("samplers_functions.R")

mydata <- read.table(file="data1.txt", header=TRUE)

m <- nrow(mydata)
X <- cbind(rep(1, m), mydata$X1, mydata$X2, mydata$X3)
Y <- mydata$Y
D <- mydata$d

nburn <- 50000
nsim <- 50000
nthin <- 10
nsize <- nsim/nthin

# Fay-Herriot
set.seed(12345)
FH_time <- system.time(FH_results <- FH_sampler(X, Y, D, nburn, nsim, nthin))
FH_theta.est <- apply(FH_results$Theta.chain, 1, mean)
FH_theta.sd <- apply(FH_results$Theta.chain, 1, sd)
FH_B.chain <- array(0, dim=c(m, nsize))
for (i in 1:nsize) FH_B.chain[,i] <- D/(FH_results$Sigma2.chain[i] + D)
FH_B.est <- apply(FH_B.chain, 1, mean)
FH_U.est <- apply(FH_results$Theta.chain - X%*%FH_results$Beta.chain,1,mean)

# Datta and Mandal
set.seed(12345)
DM_time <- system.time(DM_results <- DM_sampler(X, Y, D, nburn, nsim, nthin))
DM_theta.est <- apply(X %*% DM_results$Beta.chain + DM_results$V.chain, 1, mean)
DM_theta.sd <- apply(X %*% DM_results$Beta.chain + DM_results$V.chain, 1, sd)
DM_B.chain <- array(0, dim=c(m, nsize))
for (i in 1:nsize)
{
	 DM_B.chain[,i] <- 1 - DM_results$Sigma2.chain[i]/(D+DM_results$Sigma2.chain[i])*DM_results$p.chain[i]/(DM_results$p.chain[i] + (1-DM_results$p.chain[i])*sqrt((DM_results$Sigma2.chain[i] + D)/D)*exp(-0.5*(Y - X%*%DM_results$Beta.chain[,i])^2*DM_results$Sigma2.chain[i]/(D*(D+DM_results$Sigma2.chain[i]))))
}
DM_B.est <- apply(DM_B.chain, 1, mean)
DM_U.est <- apply(DM_results$V.chain, 1, mean)

# Horseshoe
set.seed(12345)
HS_time <- system.time(HS_results <- TPB_sampler(0.5, 0.5, X, Y, D, nburn, nsim, nthin))
HS_theta.est <- apply(X %*% HS_results$Beta.chain + HS_results$U.chain, 1, mean)
HS_theta.sd <- apply(X %*% HS_results$Beta.chain + HS_results$U.chain, 1, sd)
HS_B.est <- apply(D/(t(t(HS_results$Lambda2.chain)*HS_results$Tau2.chain) + D), 1, mean)
HS_U.est <- apply(HS_results$U.chain,1,mean)

# Strawderman-Berger
set.seed(12345)
SB_results <- TPB_sampler(1, 0.5, X, Y, D, nburn, nsim, nthin)
SB_theta.est <- apply(X %*% SB_results$Beta.chain + SB_results$U.chain, 1, mean)
SB_theta.sd <- apply(X %*% SB_results$Beta.chain + SB_results$U.chain, 1, sd)
SB_B.est <- apply(D/(t(t(SB_results$Lambda2.chain)*SB_results$Tau2.chain) + D), 1, mean)
SB_U.est <- apply(SB_results$U.chain,1,mean)

# Normal-Exponential-Gamma
set.seed(12345)
NEG_results <- TPB_sampler(1, 0.75, X, Y, D, nburn, nsim, nthin)
NEG_theta.est <- apply(X %*% NEG_results$Beta.chain + NEG_results$U.chain, 1, mean)
NEG_theta.sd <- apply(X %*% NEG_results$Beta.chain + NEG_results$U.chain, 1, sd)
NEG_B.est <- apply(D/(t(t(NEG_results$Lambda2.chain)*NEG_results$Tau2.chain) + D), 1, mean)
NEG_U.est <- apply(NEG_results$U.chain,1,mean)

# Laplace
set.seed(12345)
LA_time <- system.time(LA_results <- NG_sampler(1, X, Y, D, nburn, nsim, nthin))
LA_theta.est <- apply(X %*% LA_results$Beta.chain + LA_results$U.chain, 1, mean)
LA_theta.sd <- apply(X %*% LA_results$Beta.chain + LA_results$U.chain, 1, sd)
LA_B.est <- apply(D/(t(t(LA_results$Lambda2.chain)*LA_results$Tau2.chain) + D), 1, mean)
LA_U.est <- apply(LA_results$U.chain,1,mean)

# Normal-Gamma
set.seed(12345)
NG_results <- NG_sampler(0.5, X, Y, D, nburn, nsim, nthin)
NG_theta.est <- apply(X %*% NG_results$Beta.chain + NG_results$U.chain, 1, mean)
NG_theta.sd <- apply(X %*% NG_results$Beta.chain + NG_results$U.chain, 1, sd)
NG_B.est <- apply(D/(t(t(NG_results$Lambda2.chain)*NG_results$Tau2.chain) + D), 1, mean)
NG_U.est <- apply(NG_results$U.chain,1,mean)

# WLS
WLS_Beta.est <- solve(t(X) %*% diag(1/D) %*% X) %*%(t(X) %*% diag(1/D) %*% Y)
WLS_theta.est <- X %*% WLS_Beta.est
WLS_stdres <- (Y - WLS_theta.est)/sqrt(D)

save(m, X, Y, D, WLS_Beta.est, WLS_theta.est, HS_results, SB_results, NEG_results, LA_results, NG_results, DM_results, FH_results, file="data1_results.RData")

theta.true <- scan("data1_theta0.txt")
FH_dm <- dev_measure(FH_theta.est, theta.true)
DM_dm <- dev_measure(DM_theta.est, theta.true)
HS_dm <- dev_measure(HS_theta.est, theta.true)
SB_dm <- dev_measure(SB_theta.est, theta.true)
NEG_dm <- dev_measure(NEG_theta.est, theta.true)
LA_dm <- dev_measure(LA_theta.est, theta.true)
NG_dm <- dev_measure(NG_theta.est, theta.true)
WLS_dm <- dev_measure(WLS_theta.est, theta.true)
direct_dm <- dev_measure(Y, theta.true)

##### DIC #####
dic <- function(Y, Theta.chain, D)
{
	n <- dim(Theta.chain)[2]
	Theta.est <- apply(Theta.chain, 1, mean)
	D.chain <- rep(0, n)
	for (i in 1:n) D.chain[i] <- -2*sum(dnorm(Y, mean=Theta.chain[,i], sd=sqrt(D), log=TRUE))
	D.thetabar <- -2*sum(dnorm(Y, mean=Theta.est, sd=sqrt(D), log=TRUE))
	
	2*mean(D.chain) - D.thetabar
}
# D(\bar \theta) + 2*p_D; p_D = \bar D - D(\bar theta); D(\theta) = -2*log-likelihood
n.model <- 7
dic.all <- rep(0, n.model)

# FH
Theta.chain <- FH_results$Theta.chain
dic.all[1] <- dic(Y, Theta.chain, D)
# DM
Theta.chain <- X%*%DM_results$Beta.chain + DM_results$V.chain
dic.all[2] <- dic(Y, Theta.chain, D)
# HS
Theta.chain <- X%*%HS_results$Beta.chain + HS_results$U.chain
dic.all[3] <- dic(Y, Theta.chain, D)
# SB
Theta.chain <- X%*%SB_results$Beta.chain + SB_results$U.chain
dic.all[4] <- dic(Y, Theta.chain, D)
# NEG
Theta.chain <- X%*%NEG_results$Beta.chain + NEG_results$U.chain
dic.all[5] <- dic(Y, Theta.chain, D)
# LA
Theta.chain <- X%*%LA_results$Beta.chain + LA_results$U.chain
dic.all[6] <- dic(Y, Theta.chain, D)
# NG
Theta.chain <- X%*%NG_results$Beta.chain + NG_results$U.chain
dic.all[7] <- dic(Y, Theta.chain, D)

###### figures


# figure 5
pdf("shrinkage_factor.pdf",width=5, height=5)
plot(abs(WLS_stdres), HS_B.est, pch=16, col=1, ylim=c(0,1), xlab="|Standardized Residuals|", ylab="Shrinkage Factors")
points(abs(WLS_stdres), LA_B.est, pch=17, col=2)
points(abs(WLS_stdres), FH_B.est, pch=8, col=4)
points(abs(WLS_stdres), DM_B.est, pch=6, col=3)
legend("bottomleft", legend=c("HS", "LA", "DM", "FH"), pch=c(16, 17, 6, 8), col=c(1,2,3,4))
dev.off()

# figure 6
pdf("U_est.pdf",width=5, height=5)
plot(HS_U.est, pch=16, col="black", ylim=c(-2,6), xlab="States", ylab="")
points(LA_U.est, pch=17, col="red")
points(DM_U.est, pch=6, col="green")
points(FH_U.est, pch=8, col="blue")
abline(h=0, lty=2)
abline(v=22, lty=2)
text(22,4.6, "Massachusetts",pos=1)
legend("topright", legend=c("HS", "LA", "DM", "FH"), pch=c(16,17,6,8), col=c("black", "red", "green","blue"))
dev.off()

















