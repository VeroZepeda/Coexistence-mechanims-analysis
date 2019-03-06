#In the comments we indicate some equation numbers. The equations can be found in Zepeda and Martorell (2019) Ecology.
#This code needs:
#An object "lambdas" which is a matrix of lambdas with the years in the columns and the species in the rows. 
#An object "alfas" which is a matrix of per capita competition coefients with the associated species in the columns and the focal species in the rows. 
#An object "para" which is a vector with the values of parameter a of each focal species (Eq. 1 of the main text, Eq. S5.6).
#An object "parb" which is a vector with the values of parameter b of each focal species (Eq. 1 of the main text, Eq. S5.6).

#................................Simulations with the original model (using Eq 1. in the main text).....................................#


#Proyects population dynamics  of all species. This is the population growth model. iter= number of iterations, n0=initial number of individuals, If n0 = 999, all species are initialized to 0,1 individuals.
simorig=function(lambdas,alfas,para,parb,iter=10000,n0=999){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	if(n0==999) n0=rep(0.1,riq)
	sal=matrix(nrow=riq,ncol=iter)
	yr=floor(runif(iter)*nyr)+1
	for(i in 1:iter){
		lam=lambdas[,yr[i]]
		num=lam*n0
		den=(1+alfas%*%n0)^exp(para+parb*log(lam))
		n0=num/den
		sal[,i]=n0
		
	}
	sal
}

#Calculates long term, low-density population growth rates for a single focal species j.
calcr0=function(lambdas,alfas,para,parb,n0,j,iter=2000000){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sim=matrix(nrow=riq,ncol=iter)
	sal=1:iter
	yr=floor(runif(iter)*nyr)+1
	for(i in 1:iter){
		lam=lambdas[,yr[i]]
		num=lam*n0
		den=(1+alfas%*%n0)^exp(para+parb*log(lam))
		n0=num/den
		sim[,i]=n0
	}
	for(i in 1:iter){
		lam=lambdas[j,yr[i]]
		den=(1+alfas[j,]%*%sim[,i])^exp(para[j]+parb[j]*log(lam))
		sal[i]=log(lam/den)
	}
	sal
}


#Calculates long term, low-density population growth rates for all the species in the community.A burn-in of two thousand iterations is included. In case a species has been extincted during the transients but still can coexist with the rest in the stable state (seemingly does not occur), densities after the burn in are slightly increased to allow re-invasion 
invadiv=function(lambdas,alfas,para,parb){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sal=matrix(NA,nrow=riq,ncol=2000000)
	for(i in 1:riq){
		nn=rep(0.1,riq)
		nn[i]=0
		sim=simorig(lambdas,alfas,para,parb,iter=2000,n0=nn)
		nn=sim[,2000]+0.01
		nn[i]=0
		sal[i,]=calcr0(lambdas,alfas,para,parb,nn,i,iter=2000000)
	}
	sal
}


#Calculates long term, low-density population growth rates for all the focal species but changing the values of parameter a and parameter b for the focal species.
#para2 was calculated as follows: 
#para2 = para+parb*mean(log(lambdas))
#parb2 = rep(0,riq), riq is the number of focal species. 
invadiv2=function(lambdas,alfas,para,parb,para2,parb2){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sal=matrix(NA,nrow=riq,ncol=2000000)
	for(i in 1:riq){
		nn=rep(0.1,riq)
		nn[i]=0
		sim=simorig(lambdas,alfas,para,parb,iter=2000,n0=nn)
		nn=sim[,2000]+0.01
		nn[i]=0
		sal[i,]=calcr0(lambdas,alfas,para2,parb2,nn,i,iter=2000000)
	}
	sal
}



#.....................................................LINEAR APROXIMATION.........................................................#

#Simulates the population dynamics for a focal species by using lineal aproximation for curly E and curly C. It returns Enes, that contains population sizes of all study species, Eboni and Cboni that are the environmental and competiton standarized parameters (curly E and curly C). Time series are in the rows and species identiy in the columns. 

simulin=function(lambdas,alfas,para,parb,iter,n0){
	lambdas=log(lambdas)
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	Eprom=rowMeans(lambdas)
	gamma=parb/(parb*Eprom-1)
	Q=exp(para+parb*Eprom)
	Enes=matrix(ncol=riq,nrow=iter)
	efE=Enes
	efC=Enes
	yr=floor(runif(iter)*nyr)+1
	for(i in 1:iter){
		Ecrud=lambdas[,yr[i]]
		Ccrud=log(1+alfas%*%n0)
		Estd=Ecrud-exp(para+parb*Ecrud)*Eprom/Q
		Elin=(1-parb*Eprom)*(Ecrud-Eprom)
		Clin=Q*Ccrud-Eprom
		n0=n0*exp(Elin-Clin+gamma*Elin*Clin)
		Enes[i,]=n0
		efE[i,]=Elin
		efC[i,]=Clin
	}
	list("Enes"=Enes,"Eboni"=efE,"Cboni"=efC)
}

#Does the same as simulin but now for each study species as invader. It returns the same objects (Enes, Cboni and Eboni) but now as arrays of iterations (after a burn-in). Species (invaders) constitute the third dimension of the arrays. tira= Burn-in iterations.
simulini=function(lambdas,alfas,para,parb,iter,tira){
	riq=dim(alfas)[1]
	enes=array(NA,dim=c(iter+tira,riq,riq))
	eboni=enes
	cboni=enes
	for(i in 1:riq){
		n0=rep(0.1,riq)
		n0[i]=0
		res=simulin(lambdas,alfas,para,parb,(iter+tira),n0)
		enes[,,i]=res$Enes
		eboni[,,i]=res$Eboni
		cboni[,,i]=res$Cboni
	}
	enes=enes[-(1:tira),,]
	eboni=eboni[-(1:tira),,]
	cboni=cboni[-(1:tira),,]
	list("Enes"=enes,"Eboni"=eboni,"Cboni"=cboni)
}


#.....................................................QUADRATIC APROXIMATION.........................................................#

#Simulates the population dynamics for a focal species by using quadratic aproximation for curly E and curly C. It returns Enes, that contains population sizes of all study species, Eboni and Cboni that are the environmental and competiton standarized parameters. Time series are in the rows and species identity in the columns. 

simucuad=function(lambdas,alfas,para,parb,iter,n0){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	lambdas=log(lambdas)
	Eprom=rowMeans(lambdas)
	gamma=parb/(parb*Eprom-1)
	Q=exp(para+parb*Eprom)
	Enes=matrix(ncol=riq,nrow=iter)
	efE=Enes
	efC=Enes
	yr=floor(runif(iter)*nyr)+1
	for(i in 1:iter){
		Ecrud=lambdas[,yr[i]]
		Ccrud=log(1+alfas%*%n0)
		Estd=Ecrud-exp(para+parb*Ecrud)*Eprom/Q
		Ecuad=(1-parb*Eprom)*(Ecrud-Eprom)-0.5*parb^2*Eprom*(Ecrud-Eprom)^2
		Ccuad=Q*Ccrud-Eprom
		n0=n0*exp(Ecuad-Ccuad+gamma*Ecuad*Ccuad)
		Enes[i,]=n0
		efE[i,]=Ecuad
		efC[i,]=Ccuad
	}
	list("Enes"=Enes,"Eboni"=efE,"Cboni"=efC)
}



#The same as simulini but now with the quadratic aproximation.
simucuadi=function(lambdas,alfas,para,parb,iter,tira){
	riq=dim(alfas)[1]
	enes=array(NA,dim=c(iter+tira,riq,riq))
	eboni=enes
	cboni=enes
	for(i in 1:riq){
		n0=rep(0.1,riq)
		n0[i]=0
		res=simucuad(lambdas,alfas,para,parb,(iter+tira),n0)
		enes[,,i]=res$Enes
		eboni[,,i]=res$Eboni
		cboni[,,i]=res$Cboni
	}
	enes=enes[-(1:tira),,]
	eboni=eboni[-(1:tira),,]
	cboni=cboni[-(1:tira),,]
	list("Enes"=enes,"Eboni"=eboni,"Cboni"=cboni)
}


#This function calculates the weights qir (Eq. 5.26) for all species which are required to calculate deltaE, deltaC and deltaI. They are in matrix matq, that contains q values in each column for each invader species. It also returns Nast, that contains population sizes at the equilibrium. NOTE that Nast may contain negative values. See the published paper for the strategy followed when this happens.
calcq=function(lambdas,alfas,para,parb){
	riq=dim(alfas)[1]
	lambdas=log(lambdas)
	Eprom=rowMeans(lambdas)
	Nast=matrix(ncol=riq,nrow=(riq-1))
	matq=Nast
	Q=exp(para+parb*Eprom)
	R=exp(Eprom/Q)-1
	FI=matrix(ncol=(riq-1),nrow=(riq-1))
	for(i in 1:riq){
		submat=alfas[-i,-i]
		alfinv=alfas[i,-i]
		subR=R[-i]
		Nast[,i]=solve(submat)%*%subR
		denominvasora=alfinv%*%Nast[,i]+1
		denomresid=submat%*%Nast[,i]+1
		coefinv=Q[i]/denominvasora
		coefres=Q[-i]/denomresid
		fi=alfinv*coefinv
		for(j in 1:(riq-1)){
			FI[j,]=coefres[j]*submat[j,]
		}
		matq[,i]=fi%*%solve(FI)
	}
	list("matq"=matq, "Nast"=Nast)
}




#Decompose log-term, low density growth rate into the contributions of the different coexistence mechanisms (Eq. S5.19). datq are the wieghts (matq) obtained from the previous function. lin and cuad are the outputs of simulini and simucuadi.  The output of eq19 is a list with the vectors containing the contributions of different coexistence mechanisms for each species in the community. ri=log-term, low density growth rate.
eq19=function(lambdas,parb,lin,cuad,datq){
	riq=dim(lambdas)[1]
	lambdas=log(lambdas)
	Eprom=rowMeans(lambdas)
	gamma=parb/(parb*Eprom-1)
	EspEi=rep(0,riq)
	EspCi=EspEi
	jii=EspEi
	EspEr=rep(0,riq)
	EspCr=EspEi
	jir=EspEi
	for(i in 1:riq){
		EspEi[i]=colMeans(cuad$Eboni[,,i])[i]
		EspCi[i]=colMeans(cuad$Cboni[,,i])[i]
		jii[i]=gamma[i]*colMeans(lin$Eboni[,,i]*lin$Cboni[,,i])[i]
	}
	for(i in 1:riq){
		pagiE=cuad$Eboni[,,i]
		EspEri=colMeans(pagiE[,-i])
		pagiC=cuad$Eboni[,,i]
		EspCri=colMeans(pagiC[,-i])
		pagiji=lin$Eboni[,,i]*lin$Cboni[,,i]
		jiri=colMeans(pagiji[,-i])
		cq=datq$matq[,i]
		EspEr=EspEri%*%cq
		EspCr=EspCri%*%cq
		jir=sum(jiri*cq*gamma[-i])
	}
	
	list("deltaE"=EspEi-EspEr,"deltaC"=EspCi-EspCr,"deltaI"=jii-jir,"ri"=EspEi-EspEr-(EspCi-EspCr)+jii-jir)
}

#Calculates matrix psi (Eq. S5.27) which allow us to estimate relative non-linearity. 

calcpsi=function(lambdas,alfas,para,parb,datq){
	riq=dim(alfas)[1]
	lambdas=log(lambdas)
	Eprom=rowMeans(lambdas)
	Nast=matrix(ncol=riq,nrow=(riq-1))
	matq=Nast
	Q=exp(para+parb*Eprom)
	R=exp(Eprom/Q)-1
	FI=matrix(ncol=(riq-1),nrow=(riq-1))
	FI2=array(NA,dim=c(riq-1,riq-1,riq))
	PSI=FI2
	for(k in 1:riq){
		submat=alfas[-k,-k]
		alfinv=alfas[k,-k]
		dim(alfinv)=c(1,(riq-1))
		matalf=t(alfinv)%*%alfinv
		subR=R[-k]
		Nast[,k]=solve(submat)%*%subR
		denominvasora=as.numeric(alfinv%*%Nast[,k]+1)
		FI2[,,k]=-Q[k]/denominvasora^2*matalf
		}
	for(k in 1:riq){
		FI2k=FI2[,,k]
		FI2r=FI2[,,-k]
		cq=datq$matq[,k]
		for(i in 1:(riq-1)){
			for(j in 1:(riq-1)){
				PSI[i,j,k]=FI2k[i,j]-sum(cq*FI2r[i,j,])	
				}
			}
		}
	PSI/2
	}

#Estimates the contribution of relative non-linearity for each species in the community (Eq. S5.29). 		
calcDN=function(datpsi,lin){
	riq=dim(datpsi)[3]
	DN=rep(NA,riq)
	for(i in 1:riq){
		enes=lin$Enes[,-i,i]
		DN[i]=sum(diag(datpsi[,,i]%*%cov(enes)))
	}
	DN
}


#The following lines can be used to run the code and obtan the constributions of the coexistence mechanisms:

lin=simulini(lambdas,alfas,para,parb,1000000,10000)
cuad=simucuadi(lambdas,alfas,para,parb,1000000,10000)
coefq=calcq(lambdas,alfas,para,parb)
datpsi=calcpsi(lambdas,alfas,para,parb,coefq)

deltas=eq19(lambdas,parb,lin,cuad,coefq)
DN=calcDN(datpsi,lin)

		

