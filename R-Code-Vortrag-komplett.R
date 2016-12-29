#Verteilungsfunktion: GEV 
H_definition <- function(chi, mu, psi) {
  H <- function(x) {
    if(chi>0) y = 0
    else y = 1
    if(chi==0) {
      y = exp(-exp(-(x-mu)/psi))
    } else {
      x_new = 1+chi*(x-mu)/psi
      if(x_new > 0) {
        y = exp(-(x_new)^(-1/chi))
      }
    }
    return(y)
  }
  
  return(Vectorize(H))
}

#Frechet mit Parameter alpha
frechet <- function(alpha) {
  chi = 1/alpha
  psi = 1/alpha
  mu = 1
  return(H_definition(chi,mu,psi))
}

#Weibull mit Parameter alpha
weibull <- function(alpha) {
  chi = -1/alpha
  psi = 1/alpha
  mu = -1
  return(H_definition(chi,mu,psi))
}

#Gumbel
gumbel <- function() {
  chi = 0
  psi = 1
  mu = 0
  return(H_definition(chi,mu,psi))
}

#Numerische Ableitung (zur Bestimmung der Dichten):
ableiten <- function(f, h) {
  f_strich <- function(x) {
    y = (f(x+h/2) - f(x-h/2))/h
    return(y)
  }  
  return(Vectorize(f_strich))
}
ableiten_k <- function(f, h, k) {
  if(k<=0) return(f)
  
  else return(ableiten_k(ableiten(f,h),h,k-1))  
}



n = 4
praezision = 5000
alphas = -n:n
for(i in 1:length(alphas)) alphas[i] = 2^alphas[i]
farben = colorRampPalette(c("blue", "black", "red"))(2*n+1)
farben = c(farben[1:n],farben[(n+2):(2*n+1)])

#Plot der GEVs für verschiedene Formparameter:
plot(NaN, NaN, xlim=c(-1,3),ylim=c(0,1))
j = 1
for(i in c(-1:-n,n:1)) {
  curve(ableiten(H_definition(1/i,1,1),.01)(x), add=TRUE, col=farben[j], lw=2 , n=praezision)
  #points(1/i,1,col=farben[j])
  j = j+1
}
curve(ableiten(H_definition(0,1,1),.1)(x), add=TRUE, n=praezision, lw=2)



###########################################################################################################
###########################################################################################################
#Simulationen zum Pickands-Schätzer

m = 100
K = 500
N = K^2
par(mfcol=c(1,1))

#Gleichverteilung hat xi = -1 (Weibull)
xs <- runif(m*N,0,100)

mx <- seq(1,N)
for(i in seq(1,N)){
  mx[i] <- max(xs[((i-1)*m+1):(i*m)]) 
}

hist(mx, probability=TRUE, xlim=c(80,101), ylim=c(0,0.9), main="Histogramm der Maxima", xlab="X_n", ylab="Dichte",breaks=20)

xi <- seq(1:K)
for(k in seq(1:K)) {
  sortmx <- sort(xs[1:k^2], decreasing = TRUE)
  #Pickands-Schätzer
  xi[k] <- 1/log(2) * log((sortmx[k]-sortmx[2*k])/(sortmx[2*k]-sortmx[4*k]))
}

par(mfcol=c(1,2), ask=FALSE)

plot(xi, type='l', main="Pickands-Schätzer", xlab="k=sqrt(n)", ylab="xi")
abline(-1,0,col="blue")

hist(mx, probability=TRUE, xlim=c(90,101), ylim=c(0,0.9), main="Histogramm der Maxima", xlab="X_n", ylab="Dichte",breaks=20)
plot(ableiten(H_definition(-1,99,1),.01), add=TRUE, xlim=c(90,101), col="blue", lwd=2)


#################################################
#################################################
#Normalverteilung hat xi = 0 (Gumbel)
xs <- rnorm(m*N)

mx <- seq(1,N)
for(i in seq(1,N)){
  mx[i] <- max(xs[((i-1)*m+1):(i*m)]) 
}

xi <- seq(1:K)
for(k in seq(1:K)) {
  sortmx <- sort(xs[1:k^2], decreasing = TRUE)
  #Pickands-Schätzer
  xi[k] <- 1/log(2) * log((sortmx[k]-sortmx[2*k])/(sortmx[2*k]-sortmx[4*k]))
}


par(mfcol=c(1,2))

plot(xi, type='l', main="Pickands-Schätzer", xlab="k=sqrt(n)", ylab="xi")
abline(0,0,col="purple")

hist(mx, probability=TRUE, xlim=c(0,6), ylim=c(0,1), main="Histogramm der Maxima von Normalverteilung", xlab="X_n", ylab="Dichte")
plot(ableiten(H_definition(0,2.35,.38),.01), add=TRUE, xlim=c(0,6), col="purple", lwd=2)


#################################################
#################################################
#Cauchy-Verteilung hat xi > 0 (Frechet) (Seite 153-155)
xs <- rcauchy(m*N)

mx <- seq(1,N)
for(i in seq(1,N)){
  mx[i] <- max(xs[((i-1)*m+1):(i*m)]) 
}

xi <- seq(1:K)
for(k in seq(1:K)) {
  sortmx <- sort(xs[1:k^2], decreasing = TRUE)
  #Pickands-Schätzer
  xi[k] <- 1/log(2) * log((sortmx[k]-sortmx[2*k])/(sortmx[2*k]-sortmx[4*k]))
}


plot(xi, type='l', main="Pickands-Schätzer", xlab="k=sqrt(n)", ylab="xi")
abline(1,0,col="red")



###########################################################################################################
###########################################################################################################
#Anwendung auf Daten zu Dänischen Feuerschäden:

daten <- read.table("http://www.macs.hw.ac.uk/~mcneil/ftp/DanishData.txt", header=TRUE)

par(mfcol=c(1,1))
plot(daten$Loss.in.DKM, xlab="Schadenhöhe", ylab="", xaxt="n", type='h')

hist(daten$Loss.in.DKM, 400, probability=TRUE, xlim=c(0,50), xlab="Schadenhöhe <= 50 Mio", ylab="Häufigkeit", main="Histogramm")

plot(hist(daten$Loss.in.DKM, breaks=200, plot=FALSE)$counts, type='h', log='y', lwd=10, lend=2, col="grey", xlab="Schadenhöhe", ylab="Anzahl (log)", main="Histogramm (log)")

plot(ableiten(ecdf(daten$Loss.in.DKM),.5), xlim=c(1,270), log="x",main="Ableitung der empirischen Vtlg. (Log)", xlab="Schadenhöhe (Log)",ylab="Dichte")

xs <- daten$Loss.in.DKM
xssort <- sort(xs, decreasing=TRUE)
xssortlog <- apply(as.array(xssort), MARGIN=1, FUN=log)

#Hill-Schätzer:
alpha <- seq(1:length(xssortlog))
for(k in 1:length(xssortlog)) {
  alpha[k] <- (1/k * sum(xssortlog[1:k])-xssortlog[k])^(-1)  
}

plot(alpha[1:500], ylim=c(1,2.5), type='l', lwd=2,main="Hill-Plot",xlab="k",ylab="Hill-Schätzer")
abline(1.4,0, col="red", lw=2)

######################

library(evd)
plot(ableiten(ecdf(daten$Loss.in.DKM),.5), xlim=c(1,270), log="x",main="Vergleich: Ableitung der emp. Vtlg. mit Frechet-Vtlg.", xlab="Schadenhöhe (Log)",ylab="Dichte")
plot(function(x) dfrechet(x, shape=1.4, loc=.5), xlim=c(1,260), add=TRUE, col="red", n=100000)

qqplot(qfrechet(seq(0.001,0.999,1/2000), shape=1.4),xs,main="QQ-Plot",xlab="",ylab="")
qqline(xs, distribution=function(p) qfrechet(p,shape=1.4))

qqplot(qfrechet(seq(0.001,0.999,1/2000), shape=1.4),xs, xlim=c(0,60),ylim=c(0,60),main="QQ-Plot",xlab="",ylab="")
qqline(xs, distribution=function(p) qfrechet(p,shape=1.4))

########################## Pickands-Schätzer###############################
xe <- seq(1:(length(xssort)/4))
for(k in seq(1:(length(xssort)/4))) {
  #Pickands-Schätzer
  xe[k] <- 1/log(2) * log((xssort[k]-xssort[2*k])/(xssort[2*k]-xssort[4*k]))
}
plot(1/xe, type='l', main="Pickands-Schätzer", xlab="k", ylab="xi", ylim=c(0.5,2))
abline(1.4,0,col="dark blue")
