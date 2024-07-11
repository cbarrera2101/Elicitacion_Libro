######################################
## Capítulo 3
######################################
## 1)
# Ejemplo sobre la edad de elicitada a mano

edad<-39:56
h.edad<-c(0,0,1.5,3,6,9,12,15,15,15,9,9,9,9,6,3,0,0)

prob.edad<-h.edad/sum(h.edad)

(media.elicitada<-sum(edad*prob.edad))
(var.elicitada<-sum((edad-media.elicitada)^2*prob.edad))

# Solución vía simulación
muestra<-sample(edad,10000,replace=T,prob=prob.edad)
mean(muestra)
var(muestra)
median(muestra)
quantile(muestra, probs=c(0.05,0.5,0.95))

plot(edad,cumsum(prob.edad),type='s',ylab='Probabilidad',xlab='Edad')
points(edad,pnorm(edad,mean=media.elicitada,sd=sqrt(var.elicitada)),
       type='l')
title(main='Función de Distribución Acumulada') 

## 2)
temp<-scan()
0        500000 
0.0625   950000  
0.125    1100000 
0.25     1250000 
0.50     1500000 
0.75     1900000 
0.875    2200000 
0.9375   2450000 
1.000    3000000

temp<-matrix(temp,ncol=2,byrow=T)
plot(temp[,2],temp[,1],type='b',ylab='Cuantil',xlab='Cantidad')
abline(h=temp[,1],lty=2)
densidad<-function(temp,marcaX='XXX',marcaY='YYY'){
  k<-nrow(temp)-1
  altura<-NULL
  for(i in 1:k){
    altura<-c(altura,(temp[i+1,1]-temp[i,1])/(temp[i+1,2]-temp[i,2]))
  }
  plot(0,0,ylim=c(0,max(altura)),xlim=c(temp[1,2],temp[k+1,2]),
       ylab=marcaY,xlab=marcaX)
  
  for(i in 1:k) polygon(c(temp[i,2],temp[i,2],temp[i+1,2],
                          temp[i+1,2],temp[i,2]),c(0,altura[i],altura[i],0,0),col='grey')
}

densidad(temp,marcaX='Cantidad',marcaY='Densidad')

######################################
## Capítulo 5
######################################
## 3)
# Función auxiliar que dada una probabilidad nos da el percentil
x.s<-function(y.s,x1,x2,p1,p2){
  res<-lm(c(p1,p2)~c(x1,x2))$coefficients
  xx<-(y.s-res[1])/res[2]
  return(xx)
}

# Distribución Apriori Simulada
dist.apriori.sim<-function(x.s,p.s,Nsim=1000,n.equ=10,funcion='mean'){
  simulados<-NULL
  y<-sort(runif(Nsim*n.equ))
  for(i in 2:length(x.s)){
    y.temp<-y[y<p.s[i]]
    n<-length(y.temp)
    y<-y[-(1:n)]
    temp<-x.s(y.temp,x.s[i-1],x.s[i],p.s[i-1],p.s[i])
    simulados<-c(simulados,temp)
  }# fin if
  simulados<-matrix(sample(simulados),ncol=n.equ)
  resu<-apply(simulados,1,FUN=funcion)
  
  plot(density(resu),main='Distribución Apriori Simulada',ylab='Densidad')
  return(resu)
}# fin dist.apriori.sim

res1<-dist.apriori.sim(estatura,probas)
summary(res1)
sd(res1)

## 4)
# Determinación de los parámetros de una Normal
# a partir de una muestra de tamaño n
# representada en intervalos
#

estimaNormal <- function(parametros,frecu=frecu,limites=limites){
  media<-parametros[1]
  vari<-parametros[2]
  dev.tip<-sqrt(vari)
  proba.esti<-frecu/sum(frecu)
  error<-0.0
  for(i in 1:length(frecu)){
    error<-error+(pnorm(limites[i+1],mean=media,sd=dev.tip)-pnorm(limites[i],mean=media,sd=dev.tip)-proba.esti[i])^2
  }
  return(error)
}

limites<-c(140,160,170,180,190)
frecus<-c(10,40,45,5)

#estimaNormal(c(170,20),frecu=frecus,limites=limites)
(resu<-optim(c(170,25),estimaNormal,method='L-BFGS-B',
             lower=c(160,0.01),upper=c(185,100),
             limites=limites,frecu=frecus))

# Nivel de seguridad en nuestro conocimiento n.seg=5
param.opt<-resu$par
calcula.media.y.precision<-function(xx) c(mean(xx),1/var(xx))

n.seg<-5
Nsim<-1000
tempo<-matrix(apply(matrix(rnorm(n.seg*Nsim,mean=param.opt[1],
                                 sd=sqrt(param.opt[2])),ncol=n.seg),1,calcula.media.y.precision),
              ncol=2,byrow=T)

library(MASS)
fitdistr(tempo[,2],'gamma')
fitdistr(tempo[,1],'normal')

histograma<-function(frecu,limi){
  probas<-frecu/sum(frecu)
  altura<-probas/(limi[-1]-limi[-length(limi)])
  cotax<-c(0.90*min(limi),1.1*max(limi))
  cotay<-c(0.0,1.2*max(altura))
  plot(cotax,cotay,type='n')
  for(i in 1:length(frecu)){
    polygon(c(limi[i],limi[i],limi[i+1],limi[i+1]),c(0,altura[i],altura[i],0),col='gray')
  }
}

histograma(frecus,limites)

points(xx<-seq(min(limites),max(limites),length=100),
       dnorm(xx,mean=param.opt[1],sd=sqrt(param.opt[2])),type='l')

## 5)
ajuste.beta<-function(teta,valores=valores){
  alfa<-teta[1]
  beta<-teta[2]
  cuantil0.05<-valores[1]
  cuantil0.95<-valores[3]
  media<-valores[2]
  cuant1.teo<-qbeta(0.05,alfa,beta)
  cuant2.teo<-qbeta(0.95,alfa,beta)
  media.teo<-alfa/(alfa+beta)
  res<-(cuantil0.05-cuant1.teo)^2
  +(cuantil0.95-cuant2.teo)^2
  +(media-media.teo)^2
  return(res)
}

valores<-c(0.005,0.01,0.02)

optim(c(1,1),ajuste.beta,method="L-BFGS-B",
      lower=c(1,1)/1000000,upper=c(10,10),
      valores=valores)

## 6)
error<-NA
acumulado<-NA
media.sin<-NA
media.corr<-NA
medias<-seq(0.5,4,length=20)

for(i in medias){
  proba<-dpois(0:20,i)
  acumu<-1-sum(proba[1:6])
  media<-sum((0:5)*proba[1:6])+6.0*(1-sum(proba[1:6]))
  media.sin<-c(media.sin,media)
  error<-c(error,i-media)
  acumulado<-c(acumulado,acumu)
}
acumulado<-acumulado[-1]
error<-error[-1]
media.sin<-media.sin[-1]

# Relación entre la media y el error
plot(medias,error,xlab=expression(lambda),ylab='Error')
title(main='Error que se comete con el truncamiento \n 
en la estimación de la media')
plot(acumulado,error)

## 7)
calcula.lambda<-function(proba){
  acumu<-1-sum(proba[1:6])
  media<-sum((0:5)*proba[1:6])+6.0*(1-sum(proba[1:6]))
  media.cor<-media-0.0007817+ 0.4164908*acumu
  + 2.3313949*acumu^2
  return(media.cor)
}

# Generación de la multinomial
temp<-scan()
170    250   300   180   60   35   5 

res.multi<-rmultinom(2000,1000,temp)/1000
lambdas<-apply(res.multi,2,calcula.lambda)
hist(lambdas,freq=F,xlab=expression(lambda),
     main='Distribución Apriori',ylab='Frecuencia')
summary(lambdas)
fitdistr(lambdas,'gamma')

## 8)
(x<-rnorm(100,mean=170,sd=7.307592))
x<-seq(145,200,by=1)
y<-dnorm(x,mean=170,sd=7.307592)
y2<-dnorm(x,mean=170,sd=7.307592/sqrt(20))
plot(x,y,type='l',ylab='Densidad',xlab='Estatura (en cm)',
     ylim=c(0,0.3))
points(x,y2,type='l',lty=2)
title(main='Distribución Apriori del Promedio')
legend(160,0.31,c('Dist. Poblacional Elicitada',
                  'Dist. Apriori de la Media'),
       lty=c(1,2))

## 9)
# Elicitación de los parámetros del modelo lineal 
# via la distr. poblacional

x.s<-function(y.s,x1,x2,p1,p2){
  res<-lm(c(p1,p2)~c(x1,x2))$coefficients
  xx<-(y.s-res[1])/res[2]
  return(xx)
}

# Distribución Apriori Simulada a un nivel de x
resp.apriori.sim<-function(y.s,p.s,Nsim=1000,n.equ=10){
  simulados<-NULL
  y<-sort(runif(Nsim*n.equ))
  for(i in 2:length(y.s)){
    y.temp<-y[y<p.s[i]]
    n<-length(y.temp)
    y<-y[-(1:n)]
    temp<-x.s(y.temp,y.s[i-1],y.s[i],p.s[i-1],p.s[i])
    simulados<-c(simulados,temp)
  }# fin if
  simulados<-matrix(sample(simulados),ncol=n.equ)
  return(simulados)
}# fin resp.apriori.sim a un nivel de x

probas<-c(0,0.10,0.25,0.5,0.75,0.90,1)
# Para estatura=160
peso<-c(40,48,55,59,64,65,78)
X160<-resp.apriori.sim(peso,probas)

# Para estatura=170
peso<-c(55,63,65,68,73,75,90)
X170<-resp.apriori.sim(peso,probas)

# Para estatura=180
peso<-c(65,72,75,77,80,85,98)
X180<-resp.apriori.sim(peso,probas)

n.equ=10
YY<-cbind(X160,X170,X180)
estaturas<-c(rep(160,n.equ),rep(170,n.equ),rep(180,n.equ))

# función del modelo lineal
coefi.modelo<-function(y,x)lm.fit(cbind(1,x),y)$coefficients

Betas<-apply(YY,1,coefi.modelo,estaturas)

dim(Betas)
plot(Betas[1,],Betas[2,],
     xlab=expression(beta[0]),
     ylab=expression(beta[1]),main='Distribución Conjunta')

plot(density(Betas[1,]),main='Distribución apriori',
     xlab=expression(beta[0]),ylab='Densidad')
plot(density(Betas[2,]),main='Distribución apriori',
     xlab=expression(beta[1]),ylab='Densidad')

apply(Betas,1,mean)
apply(Betas,1,sd)
var(t(Betas))

require(ash)
x <- t(Betas)
ab <- matrix( c(-200,0,0,2), 2, 2)    # interval [-5,5) x [-5,5)
nbin <- c( 20, 20)                    # 400 bins
bins <- bin2(x, ab, nbin)             # bin counts,ab,nskip

m <- c(5,5)
f <- ash2(bins,m)
image(f$x,f$y,f$z)
contour(f$x,f$y,f$z,add=TRUE)

## 10)
temp<-scan()
0.0625  0.36 
0.125 0.42 
0.25  0.50 
0.50  0.60 
0.75  0.68 
0.875  0.75
0.9375  0.80

temp<-matrix(temp,ncol=2,byrow=T)
temp<-rbind(c(0,0),temp,c(1,1))
Area<-temp[2:9,1]-temp[1:8,1]
base<-temp[2:9,2]-temp[1:8,2]
Altura<-Area/base

plot(1,1,ylim=c(0,max(Altura)*1.3),xlim=c(0,1),type='n',
     ylab='',xlab='')
title(ylab='Densidad',xlab=substitute(pi))
title(main='Distribución Subjetiva sobre Médicos Alcohólicos')
for(i in 1:length(Altura)){
  rect(temp[i,2],0,temp[i+1,2],Altura[i],col='grey')
}

## 11)
# Determinación de los parámetros de una Normal
# a partir de de una muestra de tamaño n
# representada en intervalos

estimaNormal <- function(parametros,frecu=frecu,limites=limites){
  media<-parametros[1]
  vari<-parametros[2]
  dev.tip<-sqrt(vari)
  proba.esti<-frecu/sum(frecu)
  error<-0.0
  for(i in 1:length(frecu)){
    error<-error+(pnorm(limites[i+1],mean=media,sd=dev.tip)-
                    pnorm(limites[i],mean=media,sd=dev.tip)-proba.esti[i])^2
  }
  return(error)
}	

# Determinación de la apriori de beta0, beta1 y sigma2
# en el modelo Peso=beta0+beta1 Estatura + error

# Consideremos los hombres que miden exactamente 170cm
# De 1000 hombres en una muestra hipotética 
limites<-c(45,55,60,70,80,90)
frecus<-c(0,10,450,440,0)

(resu.170<-optim(c(70,25),estimaNormal,method='L-BFGS-B',
                 lower=c(40,0.01),upper=c(100,100),
                 limites=limites,frecu=frecus))

param.opt170<-resu.170$par	

# Consideremos los hombres que miden exactamente 180cm
# De 1000 hombres en una muestra hipotética 
limites<-c(55,65,70,75,80,85,95,110)
frecus<-c(0,10,100,400,450,35,5)

(resu.180<-optim(c(80,25),estimaNormal,method='L-BFGS-B',
                 lower=c(50,0.01),upper=c(110,100),
                 limites=limites,frecu=frecus))

param.opt180<-resu.180$par

## 12)
calcula.parame.model<-function(ii,nivel1,nivel2,
                               media1,media2,var1,var2,
                               n.seg1,n.seg2){
  xx<-c(rep(nivel1,n.seg1),rep(nivel2,n.seg2))
  yy<-c(rnorm(n.seg1,mean=media1,sd=sqrt(var1)),
        rnorm(n.seg2,mean=media2,sd=sqrt(var2)))
  res.lm<-lm(yy~xx)
  res<-c(res.lm$coefficients,mean((res.lm$residuals)^2))
  return(res)
}
# Nivel de seguridad en nuestro conocimiento n.seg170=25
n.seg170<-15

# Nivel de seguridad en nuestro conocimiento n.seg180=20
n.seg180<-10	

## 13)
# genera muestra
muestra<-function(propor,edad,n.eq){
  resu<-NA
  for(ii in 1:length(propor))resu<-rbind(resu,
                                         cbind(sample(c(1,0),n.eq[ii],replace=T,
                                                      prob=c(propor[ii],1-propor[ii])),edad[ii]))
  return(resu)
} # fin muestra

propi<-c(0.075,0.20,0.6,0.95)
neq<-c(5,10,20,20)
edad<-c(10,12,14,16)

betas<-NULL
for(i in 1:1000){
  result<- muestra(propi,edad,neq)
  yy<-result[,1]
  xx<-result[,2]
  tempo<-glm(yy~xx,family=binomial)$coefficients
  betas<-rbind(betas,tempo)
} #fin for

library(robust)
## Loading required package: fit.models
covRob(betas)
## Call:
covRob(data = betas)
solve(matrix(c(16.470, -1.15549,-1.15549,0.08185),ncol=2))

######################################
## Capítulo 6
######################################
## 14)
Construye.Apriori.Normal.IC<-function(LI,LS,n,nivel){
  media<-(LI+LS)/2
  alfa<-1-nivel
  alfa.2<-alfa/2
  gl<-n-1
  t.teo<-qt(1-alfa.2,gl)
  S2<-(LS-media)^2*n/t.teo^2
  precision<-n/S2
  LI.var<-gl*S2/qchisq(1-alfa.2,gl)
  LS.var<-gl*S2/qchisq(alfa.2,gl)
  LI.preci<-1/LS.var
  LS.preci<-1/LI.var
  
  a.mini<-function(parame,LI=LI,LS=LS,nivel=nivel){
    a0<-parame[1]
    b0<-parame[2]
    alfa<-1-nivel
    alfa.2<-alfa/2
    res<-(LI-qgamma(alfa.2,a0,rate=b0))^2+(LS-qgamma(1-alfa.2,a0,rate=b0))^2
    return(res)
  }		
  res<-optim(c(1,1),a.mini,method ="L-BFGS-B",lower=c(0.0001,0.0001),upper=c(Inf,Inf),LI=LI,LS=LS,nivel=nivel)
  alfa0<-res$par[1]
  beta0<-res$par[2]
  list(media=media,precision=precision,alfa0=alfa0,beta0=beta0)
}	

Construye.Apriori.Normal.IC(10,20,20,0.95)

## 15)
puntos<-c(800,780,600,500,490,200,150,100)
puestos<-order(puntos,decreasing=T)
un.partido<-function(x)sample(c(x[1],x[3]),1,prob=c(x[2],x[4]))

# Arg. aaaa: argumento dummy
torneo<-function(aaaa,puntos,puestos){
  n<-length(puntos)
  ciclos<-log2(n) 
  # En cada etapa se eliminan la mitad de los jugadores
  datos<-cbind(puestos,puntos)
  for(i in ciclos:1){
    n<-2^i
    n.m<-n/2
    partidos<-matrix(cbind(datos[1:n.m,],datos[n:(n.m+1),]),nrow=n.m)
    res<-apply(partidos,1,un.partido)
    datos<-datos[datos[,1] %in% res,]
  } #Fin del for
  return(res)
} # Fin función torneo

# Torneos ganados por cada jugador en 1000 torneos
(res<-table(apply(matrix(1,nrow=1000,ncol=1),1,torneo,
                  puntos,puestos)))
(proba.est<-res/sum(res))

## 16)
pies.dat<-scan()
24.2  9.4 5.5 3.0 3.2 1
21.7  8.5 6.1 3.2 2.6 0
25.4  9.6 5.5 4.0 3.1 1
25.0 10.1 5.3 3.5 2.7 1
22.0  8.5 5.7 3.1 2.7 0
25.9  9.3 6.1 4.3 3.3 1
22.2  8.6 5.2 3.9 2.9 0
21.7  8.4 5.0 3.2 2.3 0
25.5  9.2 6.1 3.3 3.2 1
24.4  9.4 4.7 3.6 2.8 1

pies.dat<-matrix(pies.dat,ncol=6,byrow=T)[,-6]
n<-nrow(pies.dat)
S<-cov(pies.dat)
x.b<-colMeans(pies.dat)
# Media apriori
mu.0<-c(1,1,1,1,1)
tau.0<-diag(rep(1,5))
nu.0<-1
alfa.0<-1

# Para. aposteriori
mu.1<-(nu.0*mu.0+n*x.b)/(nu.0+n)
tau.1<-tau.0+S+nu.0*n/(nu.0+n)*(mu.0-x.b)%*%t(((mu.0-x.b)))

library(MCMCpack)
library(MASS)
library(hdrcde)
result<-NULL
varianzas<-array(NA,dim=c(nrow(tau.1),ncol(tau.1),1000))

for(i in 1:1000){
  # Genera R de la Wishart
  R<-rwish(alfa.0+n,tau.1)
  Prec.1<-(nu.0+n)*R
  res<-mvrnorm(1,mu.1,temp<-solve(Prec.1))
  varianzas[,,i]<-temp
  result<-rbind(result,res)
}

hdr.den(result[,1],main='Distribución de la Media de la Longitud')
hdr.den(result[,2],main='Distribución de la Media de la Anchura Máx.')

## 17)
# Niveles
x1<-4
x2<-8
idas.al.baño<-0:8
frecu.x1<-c(400,350,150,50,35,10,4,1, 0)
nivel.segu1<-30

frecu.x2<-c(50,250,450,150,50,35,10,5, 0)
nivel.segu2<-10

Nsim1<-1000
medias1<-apply(matrix(sample(idas.al.baño,Nsim1*nivel.segu1,
                             replace=T,prob=frecu.x1),ncol=nivel.segu1),1,mean)
medias2<-apply(matrix(sample(idas.al.baño,Nsim1*nivel.segu2,
                             replace=T,prob=frecu.x2),ncol=nivel.segu2),1,mean)

library(MASS)
fitdistr(medias1,'gamma')
fitdistr(medias2,'gamma')

Nsim2<-1000
medias11<-rgamma(Nsim2,23.270277,rate=22.476949)
medias22<-rgamma(Nsim2,30.351702,rate=14.5711617)


temporal<-function(medias,n.s){
  y<-c(rpois(n.s[1],medias[1]),rpois(n.s[2],medias[2]))
  x<-c(rep(4,n.s[1]),rep(8,n.s[2]))
  resu<-glm(y~x,family='poisson')$coefficients
  return(resu)
}

resu<-apply(cbind(medias11,medias22),1,temporal,c(nivel.segu1,nivel.segu2))
dim(resu)
colMeans(t(resu))
cov(t(resu))	

colMeans(t(resu))
cov(t(resu))

######################################
## Capítulo 8
######################################
## 18)
# Simulación de la Distribución Predictiva
# Nro. promedio de goles del campeonato colombiano
# Apriori: Normal(2.5, 0.20^2)
# Dist. muestral: Poisson(theta)

dist.predictiva<-function(Nsim=10000,media=2.5,dt=0.20){
  res<-rpois(Nsim,rnorm(Nsim,mean=media,sd=dt))
  tabla<-table(res)
  print(tabla/sum(tabla))
  barplot(tabla)
  print(summary(res))
  print(var(res))
}

dist.predictiva()
title(main='Distribución Predictiva\n Campeonato Colombiano',
      xlab='Número de Goles',
      ylab='Frecuencia')
legend(7,2000,c('Media apriori=2.5','Desv. Est. apriori=0.20'))

res

######################################
## Capítulo 9
######################################
## 19)
#############################################################
# Elicitación de la distribución apriori de una proporción
# Se asume que la apriori es una Beta(alfa,beta)
# El procedimiento consiste en determinar los límites más extremos plausibles
# o sea, donde estamos casi convencidos cae la verdadera proporción y el valor
# de la moda.
# Programador JCCM
# Fecha de última modificación: Agosto 27, 2009
conseguir.beta.apriori<-
  function(){
    # Entrada de la información del elicitador
    LI<-as.numeric(readline("Cuál es el menor valor posible para la proporción? "))
    LS<-as.numeric(readline("Cuál es el mayor valor posible para la proporción? "))
    moda<-as.numeric(readline("Cuál es el valor más posible para la proporción? "))
    
    funcion.a.mini<-function(parametros,moda=0.5,mini=0.3,maxi=0.6){
      a1<-parametros[1]
      b1<-parametros[2]
      temp<-abs((a1-1)/(a1+b1-2)-moda)+abs(pbeta(mini,a1,b1)-0.01)+abs(pbeta(maxi,a1,b1)-0.99)
      return(temp)
    }
    res<-optim(c(2,2),funcion.a.mini,moda=moda,mini=LI,maxi=LS,lower=c(1,1),method='L')
    parametros<-res$par
    alfa<-parametros[1]
    beta<-parametros[2]
    x<-seq(0.01,0.99,length=200)
    y<-dbeta(x,alfa,beta)
    plot(x,y,type='l',ylab='',xlab='')
    abline(v=c(LI,moda,LS),lty=2)
    return(res)
  } 

# Esta función utiliza la aproximación  de Perry y Greig (1975)
# Fin de función principal

conseguir.beta.apriori<-
  function(){
    # Entrada de la información del elicitador
    LI<-as.numeric(readline("Cuál es el menor valor posible para la proporción? "))
    LS<-as.numeric(readline("Cuál es el mayor valor posible para la proporción? "))
    moda<-as.numeric(readline("Cuál es el valor más posible para la proporción? "))
    sigma<-(LS-LI)/3.25
    mu<-(LI+0.95*moda+LS)/2.95
    print(mu)
    print(sigma)
    alfa<-(mu^2*(1-mu))/sigma^2-mu
    print(alfa)
    beta<-alfa*(1-mu)/mu
    print(beta)
    x<-seq(0.01,0.99,length=200)
    y<-dbeta(x,alfa,beta)
    plot(x,y,type='l',ylab='',xlab='')
    abline(v=c(LI,moda,LS),lty=2)
    abline(h=0)
    list(alfa=alfa,beta=beta,media=mu,desvi.tipica=sigma)
  } 

# Fin de función principal

resultados<-conseguir.beta.apriori()

resultados
alfa0<-resultados$par[1]
beta0<-resultados$par[2]
x<-seq(0,1,length=1000)
y<-dbeta(x,alfa0,beta0)
par(mfrow=c(3,1))
plot(x,y,ylab='Densidad',xlab=expression(pi),type='l')
title(main='Distribución Apriori')
# Saco una muestra de 10 estudiantes al azar y obtengo 2 mujeres.
y2<-x^2*(1-x)^8
# Verosimilitud
plot(x,y2,type='l',col='red',ylab='',xlab='')
title(main='Función de Verosimilitud')
# aposteriori
y3<-dbeta(x,alfa0+2,beta0+8)
plot(x,y3,ylab='Densidad',xlab=expression(pi),type='l')
title(main='Distribución Aposteriori')

