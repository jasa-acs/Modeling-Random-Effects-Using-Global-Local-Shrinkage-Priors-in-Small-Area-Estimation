rm(list=ls())

library(MASS)
library(GIGrvg)

source("samplers_functions.R")

mydata <- read.table(file="data2.txt", header=TRUE)

Y <- mydata$y
m <- nrow(mydata)
X <- cbind(rep(1, m), mydata$x)
D <- mydata$D

nburn <- 50000
nsim <- 50000
nthin <- 10
nsize <- nsim/nthin

set.seed(12345)
DM_time <- system.time(DM_results <- DM_sampler(X, Y, D, nburn, nsim, nthin))

set.seed(12345)
LA_time <- system.time(LA_results <- NG_sampler(1, X, Y, D, nburn, nsim, nthin))

set.seed(12345)
HS_time <- system.time(HS_results <- TPB_sampler(0.51, 0.5, X, Y, D, nburn, nsim, nthin))

set.seed(12345)
FH_time <- system.time(FH_results <- FH_sampler(X, Y, D, nburn, nsim, nthin))

save(m, X, Y, D, FH_time, DM_time, HS_time, LA_time, FH_results, DM_results, LA_results, HS_results, file="data2_results.RData")

load("data2_results.RData")

# dic
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

n.model <- 4
dic.all <- rep(0, n.model)

# DM
Theta.chain <- X%*%DM_results$Beta.chain + DM_results$V.chain
dic.all[1] <- dic(Y, Theta.chain, D)
# LA
Theta.chain <- X%*%LA_results$Beta.chain + LA_results$U.chain
dic.all[2] <- dic(Y, Theta.chain, D)
# HS
Theta.chain <- X%*%HS_results$Beta.chain + HS_results$U.chain
dic.all[3] <- dic(Y, Theta.chain, D)
# FH
Theta.chain <- FH_results$Theta.chain
dic.all[4] <- dic(Y, Theta.chain, D)

##### figures
library("maptools")
library("fields")
library("RColorBrewer")

ctyid1 <- mydata$FIPS

shp <- readShapePoly("cb_2013_us_county_20m/cb_2013_us_county_20m.shp")
shp_div <- readShapePoly("cb_2013_us_division_20m/cb_2013_us_division_20m.shp")


CO <- as.character(shp@data$COUNTYFP)
ST <- as.character(shp@data$STATEFP)
ctyid2 <- as.numeric(paste(ST, CO, sep=""))
shp@data$FIPS <- ctyid2

plotid <- match(ctyid1, ctyid2)

shp2 <- shp[plotid,] # shape file that only contains 3141 counties; the plot order is the same as in data2.txt

ak_id <- shp2@data$STATEFP=="02"
hi_id <- shp2@data$STATEFP=="15"
cont_id <- shp2@data$STATEFP!="02" & shp2@data$STATEFP!="15"

proj4string(shp2) <- CRS("+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")
shp_proj <- spTransform(shp2, CRS("+init=epsg:3349"))

proj4string(shp_div) <- CRS("+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")
shp_div_proj <- spTransform(shp_div, CRS("+init=epsg:3349"))

DM_theta.est <- apply(X %*% DM_results$Beta.chain + DM_results$V.chain, 1, mean)
DM_delta.est <- apply(DM_results$Delta.chain, 1, mean)

LA_theta.est <- apply(X %*% LA_results$Beta.chain + LA_results$U.chain, 1, mean)
LA_U.est <- apply(LA_results$U.chain,1,mean)


# Figure 7: map of theta
ncolors <- 50
plotclr <- heat.colors(ncolors, alpha = 1)

plotvar <- c(LA_theta.est)
breakpoints <- quantile(plotvar, probs=seq(from=0, to=1, length=ncolors+1))
colornum <- findInterval(plotvar, breakpoints, all.inside=T)
colcode <- plotclr[colornum]

pdf(file = "theta_map1_all_EXP_div.pdf", width = 8, height = 5)
par(mar=c(1,1,1,1))
layout_mat <- matrix(1, nrow=6, ncol=8)
layout_mat[5:6,1:2] <- 2
layout_mat[6,3] <- 3
layout_mat[,8] <- 4
layout(layout_mat,widths=c(rep(1, 7), 0.3))
plot(shp_proj[cont_id,], col=colcode[cont_id], lwd=0.1)
plot(shp_div_proj, add=TRUE, lwd=1.5)
plot(shp_proj[ak_id,], col=colcode[ak_id], lwd=0.1)
plot(shp_div_proj, add=TRUE, lwd=1.5)
plot(shp_proj[hi_id,], col=colcode[hi_id], lwd=0.1)
plot(shp_div_proj, add=TRUE, lwd=1.5)
image.scale(z=plotclr, zlim=range(plotvar), col=plotclr, breaks=breakpoints, horiz=F)
dev.off()

# Figure 8: map of absolute value of random effects
ncolors <- 9

plotvar <- abs(LA_U.est)
breakpoints <- quantile(plotvar, probs=seq(from=0, to=1, length=ncolors+1))
plotclr <- rev(brewer.pal(ncolors, "OrRd"))
colornum <- findInterval(plotvar, breakpoints, all.inside=T)
colcode <- plotclr[colornum]

pdf(file = "u_map1_all_EXP_div.pdf", width = 8, height = 5)
par(mar=c(1,1,1,1))
layout_mat <- matrix(1, nrow=6, ncol=7)
layout_mat[5:6,1:2] <- 2
layout_mat[6,3] <- 3
layout(layout_mat)
plot(shp_proj[cont_id,], col=colcode[cont_id], lwd=0.1)
plot(shp_div_proj, add=TRUE, lwd=1.5)
legend("bottomright",legend=leglabs(round(breakpoints,3)),fill=plotclr,bty="n", horiz=FALSE)
plot(shp_proj[ak_id,], col=colcode[ak_id], lwd=0.1)
plot(shp_div_proj, add=TRUE, lwd=1.5)
plot(shp_proj[hi_id,], col=colcode[hi_id], lwd=0.1)
plot(shp_div_proj, add=TRUE, lwd=1.5)
dev.off()

range(LA_theta.est)
LA_med <- median(LA_theta.est)
quantile(LA_theta.est, probs=c(0.25, 0.75))
shp_proj@data[which.max(LA_theta.est),]
shp_proj@data[which.min(LA_theta.est),]

state_id <- floor(shp_proj@data$FIPS/1000)
states <- unique(state_id)
nstate <- length(states)

gr_q3_rate <- rep(0, nstate)
le_q1_rate <- rep(0, nstate)

for (i in 1:nstate)
{
	counties <- state_id==states[i]
	ncounty <- sum(counties)
	gr_q3_rate[i] <- sum(LA_theta.est[counties] > 0.1886) / ncounty
	le_q1_rate[i] <- sum(LA_theta.est[counties] < 0.1111) / ncounty
}

o <- order(gr_q3_rate)
gr_q3_rate[o]
states[o]

o <- order(le_q1_rate)
le_q1_rate[o]
states[o]