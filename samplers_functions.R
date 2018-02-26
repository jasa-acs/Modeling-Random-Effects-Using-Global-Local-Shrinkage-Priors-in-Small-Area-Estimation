
DM_sampler <- function(X, Y, D, nburn, nsim, nthin)
{
	# parameters
	m <- nrow(X)
	np <- ncol(X)
	a <- 2*mean(D)
	b <- 3
	c <- 1
	d <- 4
	
	# initial values
	V <- rep(1,m)
	Sigma2 <- 1
	Delta <- rep(1,m)
	Beta <- rep(1,np)
	p <- 0.5
	
	# chain containers
	V.chain <- array(0, dim=c(m,nsim/nthin))
	Sigma2.chain <- rep(1, nsim/nthin)
	Delta.chain <- array(1,dim=c(m,nsim/nthin))
	Beta.chain <- array(0, dim=c(np,nsim/nthin))
	p.chain <- rep(0, nsim/nthin)
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		# update Sigma2
		Sigma2 <- 1/rgamma(1,shape=b+sum(Delta)/2,rate=a+sum(Delta*V*V)/2)
		
		# update p
		p <- rbeta(1, shape1=c+sum(Delta), shape2=d+m-sum(Delta))
		
		# update Delta
		p.delta <- p/(p+(1-p)*sqrt((Sigma2+D)/D)*exp(-0.5*(Y-X%*%Beta)^2*Sigma2/(D+Sigma2)/D))
		for (i in 1:m) Delta[i] <- rbinom(1,1,p.delta[i])
		
		# update Beta
		Xstd <- X/sqrt(D)
		sigma <- solve(t(Xstd)%*%Xstd)
		mean <- apply(X*(Y-Delta*V)/D,2,sum)
		mean <- sigma%*%mean
		Beta <- mvrnorm(1, mu=mean,Sigma=sigma)
		
		# update V
		mean.v <- Sigma2*(Y-X%*%Beta)/(Sigma2+D)
		var.v <- Sigma2*D/(Sigma2+D)
		V <- rnorm(m, mean=mean.v, sd=sqrt(var.v))*Delta
		
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			Sigma2.chain[(index-nburn)/nthin] <- Sigma2
			Delta.chain[,(index-nburn)/nthin] <- Delta
			p.chain[(index-nburn)/nthin] <- p
			V.chain[,(index-nburn)/nthin] <- V
			Beta.chain[,(index-nburn)/nthin] <- Beta		
		}	
	}
	
	list(Beta.chain=Beta.chain, V.chain=V.chain, Delta.chain=Delta.chain, p.chain=p.chain, Sigma2.chain=Sigma2.chain)
}

CDM_sampler <- function(alpha, X, Y, D, nburn, nsim, nthin)
{
	# parameters
	m <- nrow(X)
	np <- ncol(X)
	
	# initial values
	Theta <- rep(1,m)
	A <- c(0.1, 0.5)
	Delta <- rep(1,m)
	Beta <- rep(1,np)
	p <- 0.5
	
	# chain containers
	V.chain <- array(0, dim=c(m,nsim/nthin))
	A.chain <- array(1, dim=c(2,nsim/nthin))
	Delta.chain <- array(1,dim=c(m,nsim/nthin))
	Beta.chain <- array(0, dim=c(np,nsim/nthin))
	p.chain <- rep(0, nsim/nthin)
	Theta.chain <- array(0, dim=c(m,nsim/nthin))
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		# update Theta
		mean.theta <- ((X*D)%*%Beta + A[2-Delta]*Y)/(D + A[2-Delta])
		sd.theta <- sqrt((D*A[2-Delta])/(D+A[2-Delta]))
		Theta <- rnorm(m, mean=mean.theta, sd=sd.theta)
		
		# update Beta
		Xstd <- X/sqrt(A[2-Delta])
		sigma.beta <- solve(t(Xstd)%*%Xstd)
		mean.beta <- sigma.beta%*%t(Xstd)%*%(Theta/sqrt(A[2-Delta]))
		Beta <- mvrnorm(1, mu=mean.beta,Sigma=sigma.beta)
		
		# update p
		p <- rbeta(1, shape1=1+sum(Delta), shape2=1+m-sum(Delta))
		
		# update Delta
		ppart1 <- p*exp(-(Theta - X%*%Beta)^2/(2*A[1]))/sqrt(A[1])
		ppart2 <- (1-p)*exp(-(Theta - X%*%Beta)^2/(2*A[2]))/sqrt(A[2])
		p.delta <- ppart1/(ppart1+ppart2)
		Delta <- rbinom(m,1,p.delta)
		
		# update A1 A2
		A[1] <- 1/rgamma(1, shape=alpha[1]+sum(Delta)/2-1, rate=sum(Delta*(Theta - X%*%Beta)^2)/2)
		while(A[1] > A[2]) A[1] <- 1/rgamma(1, shape=alpha[1]+sum(Delta)/2-1, rate=sum(Delta*(Theta - X%*%Beta)^2)/2)
	
		A[2] <- 1/rgamma(1, shape=alpha[2]+sum(1-Delta)/2-1, rate=sum((1-Delta)*(Theta - X%*%Beta)^2)/2)
		while(A[1] > A[2]) A[2] <- 1/rgamma(1, shape=alpha[2]+sum(1-Delta)/2-1, rate=sum((1-Delta)*(Theta - X%*%Beta)^2)/2)
	
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			A.chain[,(index-nburn)/nthin] <- A
			Delta.chain[,(index-nburn)/nthin] <- Delta
			p.chain[(index-nburn)/nthin] <- p
			Theta.chain[,(index-nburn)/nthin] <- Theta
			Beta.chain[,(index-nburn)/nthin] <- Beta
			V.chain[,(index-nburn)/nthin] <- Theta - X%*%Beta
		}	
	}
	
	list(Beta.chain=Beta.chain, V.chain=V.chain, Delta.chain=Delta.chain, p.chain=p.chain, A.chain=A.chain)
}

FH_sampler <- function(X, Y, D, nburn, nsim, nthin)
{
	# parameters
	m <- nrow(X)
	np <- ncol(X)
	
	# initial values
	Theta <- rep(1,m)
	Beta <- rep(1, np)
	Sigma2 <- 1
	
	# chain containers
	Theta.chain <- array(0, dim=c(m,nsim/nthin))
	Beta.chain <- array(0, dim=c(np,nsim/nthin))
	Sigma2.chain <- rep(0, nsim/nthin)
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		
		# update Sigma2
		Sigma2 <- 1/rgamma(1,shape=m/2-1, rate=sum((Theta-X%*%Beta)^2)/2)
		
		# update Beta
		mean.Beta <- solve(t(X)%*%X)%*%t(X)%*%Theta
		var.Beta <- Sigma2*solve(t(X)%*%X)
		Beta <- mvrnorm(1, mu=mean.Beta, Sigma=var.Beta)
		
		# update Theta
		var.Theta <- 1/(1/Sigma2+1/D)
		mean.Theta <- var.Theta*(Y/D+X%*%Beta/Sigma2)
		Theta <- rnorm(m, mean=mean.Theta, sd=sqrt(var.Theta))
		
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			Sigma2.chain[(index-nburn)/nthin] <- Sigma2
			Theta.chain[,(index-nburn)/nthin] <- Theta
			Beta.chain[,(index-nburn)/nthin] <- Beta		
		}	
	}
	
	list(Beta.chain=Beta.chain, Theta.chain=Theta.chain, Sigma2.chain=Sigma2.chain)
}

TPB_sampler <- function(p, q, X, Y, D, nburn, nsim, nthin)
{
	# parameters
	m <- nrow(X)
	np <- ncol(X)
	
	a <- 10^-10
	b <- 10^-10
	phi <- 1
	
	# initial values
	Lambda2 = rep(1,m)
	Tau2 = 1
	Z = rep(1,m)
	Beta = rep(1,np)
	U = rep(1,m)
	
	# chain containers
	Lambda2.chain = array(0, dim=c(m,nsim/nthin))
	Tau2.chain = rep(0, nsim/nthin)
	Z.chain = array(0,dim=c(m,nsim/nthin))
	Beta.chain = array(0, dim=c(np,nsim/nthin))
	U.chain = array(0, dim=c(m,nsim/nthin))
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		# update Tau2
		Tau2 <- 1/rgamma(1,shape=(a+m)/2,rate=sum(U^2/Lambda2)/2+b/2)
		
		# update Z
		Z <- rgamma(m, shape=p+q, rate=phi+Lambda2)
		
		# update Lambda2
		for (i in 1:m) Lambda2[i] = rgig(1,p-0.5,U[i]^2/Tau2,2*Z[i])
		
		# update Beta
		Xstd <- X/sqrt(D)
		sigma <- solve(t(Xstd)%*%Xstd)
	    mean <- apply(X*(Y-U)/D,2,sum)
		mean <- sigma%*%mean
		Beta <- mvrnorm(1, mu=mean,Sigma=sigma)
		
		# update U
		sigma2 <- 1/(1/D+1/Lambda2/Tau2)
		mean <- sigma2 * (Y-X%*%Beta)/D
		U <- rnorm(m, mean=mean, sd=sqrt(sigma2))
	
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			Tau2.chain[(index-nburn)/nthin] <- Tau2
			Lambda2.chain[,(index-nburn)/nthin] <- Lambda2
			Z.chain[,(index-nburn)/nthin] <- Z
			U.chain[,(index-nburn)/nthin] <- U
			Beta.chain[,(index-nburn)/nthin] <- Beta		
		}	
	}
	
	list(Beta.chain=Beta.chain, U.chain=U.chain, Tau2.chain=Tau2.chain, Lambda2.chain=Lambda2.chain)
}

NG_sampler <- function(a, X, Y, D, nburn, nsim, nthin)
{
	m <- nrow(X)
	np <- ncol(X)
	
	b <- 1
	p <- 10^-10
	q <- 10^-10
	
	Lambda2 <- rep(1,m)
	Tau2 <- 1
	Beta <- rep(1,np)
	U <- rep(1,m)
		
	Lambda2.chain <- array(0, dim=c(m,nsim/nthin))
	Tau2.chain <- rep(0, nsim/nthin)
	Beta.chain <- array(0, dim=c(np,nsim/nthin))
	U.chain <- array(0, dim=c(m,nsim/nthin))
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		# update Tau2
		Tau2 <- 1/rgamma(1,shape=(p+m)/2,rate=sum(U^2/Lambda2)/2+q/2)
			
		# update Lambda2
		for (i in 1:m) Lambda2[i] <- rgig(1,a-0.5,U[i]^2/Tau2,2*b)
		
		# update Beta
		Xstd <- X/sqrt(D)
		sigma <- solve(t(Xstd)%*%Xstd)
	    	mean <- apply(X*(Y-U)/D,2,sum)
		mean <- sigma%*%mean
		Beta <- mvrnorm(1, mu=mean,Sigma=sigma)
		
		# update U
		sigma2 <- 1/(1/D+1/Lambda2/Tau2)
		mean <- sigma2 * (Y-X%*%Beta)/D
		U <- rnorm(m, mean=mean, sd=sqrt(sigma2))
	
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			Tau2.chain[(index-nburn)/nthin] <- Tau2
			Lambda2.chain[,(index-nburn)/nthin] <- Lambda2
			U.chain[,(index-nburn)/nthin] <- U
			Beta.chain[,(index-nburn)/nthin] <- Beta		
		}	
	}
	
	list(Beta.chain=Beta.chain, U.chain=U.chain, Tau2.chain=Tau2.chain, Lambda2.chain=Lambda2.chain)
	
}

dev_measure <- function(thetahat, theta0)
{
	AAD <- mean(abs(thetahat - theta0))
	ASD <- mean((thetahat-theta0)^2)
	ARB <- mean(abs(thetahat-theta0)/abs(theta0))
	ASRB <- mean((thetahat-theta0)^2/theta0^2)
	
	list(AAD=AAD, ASD=ASD, ARB=ARB, ASRB=ASRB)
}

image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}
