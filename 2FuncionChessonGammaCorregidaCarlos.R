###Última actualización Enero 2019
#ANÁLISIS DEL MODELO SIN APROXIMACIÓN CUADRÁTICA
#Simulaciones con el modelo original
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

#calcula tasa de crecimiento de largo plazo y baja densidad para una especie j
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


#Hace lo que calcr0 para las 19 especies
invadiv=function(lambdas,alfas,para,parb){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sal=matrix(NA,nrow=riq,ncol=2000000)
	for(i in 1:riq){
		nn=rep(0.1,riq)
		nn[i]=0
		sim=simorig(lambdas,alfas,para,parb,iter=2000,n0=nn)
		nn=sim[,2000]+0.1
		nn[i]=0
		sal[i,]=calcr0(lambdas,alfas,para,parb,nn,i,iter=2000000)
	}
	sal
}


#Hace lo que calcr0 para las 19 especies, pero usando valores de a y b diferentes para la especie invasora.
invadiv2=function(lambdas,alfas,para,parb,para2,parb2){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sal=matrix(NA,nrow=riq,ncol=2000000)
	for(i in 1:riq){
		nn=rep(0.1,riq)
		nn[i]=0
		sim=simorig(lambdas,alfas,para,parb,iter=2000,n0=nn)
		nn=sim[,2000]+0.1
		nn[i]=0
		sal[i,]=calcr0(lambdas,alfas,para2,parb2,nn,i,iter=2000000)
	}
	sal
}



#################################################
#ANÁLISIS DE LA APROXIMACIÓN CUADRÁTICA AL MODELO


#simula la dinámica de la comunidad usando la aproximación lineal. Regresa Enes, que tiene los tamaños poblacionales de todas las especies, con las especies en las columnas y el tiempo en los renglones. Eboni y Cboni contienen los parametros ambiental y de competencia estandarizados y con aproximación lineal en el mismo formato.
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


#simula la dinámica de la comunidad usando la aproximación cuadrática. Regresa Enes, que tiene los tamaños poblacionales de todas las especies, con las especies en las columnas y el tiempo en los renglones. Eboni y Cboni contienen los parametros ambiental y de competencia estandarizados y con aproximación cuadrática en el mismo formato.
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

#hace lo mismo que simulin, pero para todas las posibles especies como invasoras. Regresa los mismos objetos, pero ahora son arrays de iter iteraciones (después de un burn-in de tira ietarciones) POR especies POR cada especie como invasora.
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

#como simulini, pero con aproximación cuadrática.
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


#calcula los coeficientes q para todas las especies. Están contenidos en matq, que contiene en cada columna los valores de q para cada especie invasora. Adiocionalmente regresa Nast, que contiene los tamaños poblacionales en el equilibrio.
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

eq18=function(lambdas,parb,lin,cuad){
	riq=dim(lambdas)[1]
	lambdas=log(lambdas)
	Eprom=rowMeans(lambdas)
    gamma=parb/(parb*Eprom-1)
	ri=rep(0,riq)
	for(i in 1:riq){
		EspE=colMeans(cuad$Eboni[,,i])
		EspC=colMeans(cuad$Cboni[,,i])
		ji=colMeans(lin$Eboni[,,i]*lin$Cboni[,,i])
		ri[i]=EspE[i]-EspC[i]+gamma[i]*ji[i]
	}
	ri	
}

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
		
calcDN=function(datpsi,lin){
	riq=dim(datpsi)[3]
	DN=rep(NA,riq)
	for(i in 1:riq){
		enes=lin$Enes[,-i,i]
		DN[i]=sum(diag(datpsi[,,i]%*%cov(enes)))
	}
	DN
}





		

