

p = rpois(1000,500)
b = rbinom(1000, 10, 0.1)
e = rexp(1000,10)
n = rnorm(1000,500)
print(var(p))
print(mean(p))
print(var(b))
print(mean(b))
print(var(e))
print(mean(e))
print(var(n))
print(mean(n))

#Ex 2 si 3 graficele functiilor de masa, densitate si functiile de repartitie

##graficele poisson
p1 = data.frame(DENSITATE=dpois(0:20, 0.1), MASA=ppois(0:20, 0.1))
p2 = data.frame(DENSITATE=dpois(0:20, 2), MASA=ppois(0:20, 2))
p3 = data.frame(DENSITATE=dpois(0:20, 5), MASA=ppois(0:20, 5))
p4 = data.frame(DENSITATE=dpois(0:20, 500), MASA=ppois(0:20, 500))
p5 = data.frame(DENSITATE=dpois(0:20, 0.01), MASA=ppois(0:20, 0.01))

plot(p1$DENSITATE, type="o",col="red", main="Poisson Ddensitate") 
lines(p2$DENSITATE,type="o",col="green")
lines(p3$DENSITATE,type="o",col="blue")
lines(p4$DENSITATE,type="o",col="yellow")
lines(p5$DENSITATE,type="o",col="coral")
legend("topright",c("l=0.1","l=2","l=5","l=500", "l=0.01"), col=c("red","green","blue","yellow","coral"),pch=15)

plot(p1$MASA,type="o",col="red", main="Poisson Masa") 
lines(p2$MASA,type="o",col="green")
lines(p3$MASA,type="o",col="blue")
lines(p4$MASA,type="o",col="yellow")
lines(p5$MASA,type="o",col="coral")
legend("bottomright",c("l=0.1","l=2","l=5","l=500", "l=0.01"), col=c("red","green","blue","yellow","coral"),pch=15)

##graficele binomial
b1 = data.frame(DENSITATE=dbinom(0:20, 20, 0.05), MASA=pbinom(0:20, 20, 0.05))
b2 = data.frame(DENSITATE=dbinom(0:20, 20, 0.1), MASA=pbinom(0:20, 20, 0.1))
b3 = data.frame(DENSITATE=dbinom(0:20, 20, 0.25), MASA=pbinom(0:20, 20, 0.25))
b4 = data.frame(DENSITATE=dbinom(0:20, 20, 0.5), MASA=pbinom(0:20, 20, 0.5))
b5 = data.frame(DENSITATE=dbinom(0:20, 20, 0.75), MASA=pbinom(0:20, 20, 0.75))

plot(b1$DENSITATE, type="o",col="red", main="Binomial Ddensitate")
lines(b2$DENSITATE,type="o",col="green")
lines(b3$DENSITATE,type="o",col="blue")
lines(b4$DENSITATE,type="o",col="yellow")
lines(b5$DENSITATE,type="o",col="coral")
legend("topright",c("p=0.05","p=0.1","p=0.25","p=0.5", "p=0.75"), col=c("red","green","blue","yellow","coral"),pch=15)

plot(b1$MASA,type="o",col="red", main="Binomial Masa")
lines(b2$MASA,type="o",col="green")
lines(b3$MASA,type="o",col="blue")
lines(b4$MASA,type="o",col="yellow")
lines(b5$MASA,type="o",col="coral")
legend("bottomright",c("p=0.05","p=0.1","p=0.25","p=0.5", "p=0.75"), col=c("red","green","blue","yellow","coral"),pch=15)

##graficele exponential
e1 = data.frame(DENSITATE=dexp(0:20, 0.25), MASA=pexp(0:20, 0.25))
e2 = data.frame(DENSITATE=dexp(0:20, 0.75), MASA=pexp(0:20, 0.75))
e3 = data.frame(DENSITATE=dexp(0:20, 1), MASA=pexp(0:20, 1))
e4 = data.frame(DENSITATE=dexp(0:20, 5), MASA=pexp(0:20, 5))
e5 = data.frame(DENSITATE=dexp(0:20, 20), MASA=pexp(0:20, 20))

plot(e1$DENSITATE, type="o",col="red", main="Exponential Ddensitate")
lines(e2$DENSITATE,type="o",col="green")
lines(e3$DENSITATE,type="o",col="blue")
lines(e4$DENSITATE,type="o",col="yellow")
lines(e5$DENSITATE,type="o",col="coral")
legend("topright",c("r=0.25","r=0.75","r=1","r=5", "r=20"), col=c("red","green","blue","yellow","coral"),pch=15)

plot(e1$MASA,type="o",col="red", main="Exponential Masa") 
lines(e2$MASA,type="o",col="green")
lines(e3$MASA,type="o",col="blue")
lines(e4$MASA,type="o",col="yellow")
lines(e5$MASA,type="o",col="coral")
legend("bottomright",c("r=0.25","r=0.75","r=1","r=5", "r=20"), col=c("red","green","blue","yellow","coral"),pch=15)

##graficele normal
n1 = data.frame(DENSITATE=dnorm(0:20, 0, 1), MASA=pnorm(0:20, 0, 1))
n2 = data.frame(DENSITATE=dnorm(0:20, 0, 5), MASA=pnorm(0:20, 0, 5))
n3 = data.frame(DENSITATE=dnorm(0:20, 5, 1), MASA=pnorm(0:20, 5, 1))
n4 = data.frame(DENSITATE=dnorm(0:20, 10, 1), MASA=pnorm(0:20, 10, 1))
n5 = data.frame(DENSITATE=dnorm(0:20, 2, 10), MASA=pnorm(0:20, 2, 10))

plot(n1$DENSITATE, type="o",col="red", main="Normal Ddensitate")
lines(n2$DENSITATE,type="o",col="green")
lines(n3$DENSITATE,type="o",col="blue")
lines(n4$DENSITATE,type="o",col="yellow")
lines(n5$DENSITATE,type="o",col="coral")
legend("topright",c("m=0, sd=1","m=0, sd=5","m=5, sd=1","m=10,sd=1", "m=2,sd=10"), col=c("red","green","blue","yellow","coral"),pch=15)

plot(n1$MASA,type="o",col="red", main="Normal Masa") 
lines(n2$MASA,type="o",col="green")
lines(n3$MASA,type="o",col="blue")
lines(n4$MASA,type="o",col="yellow")
lines(n5$MASA,type="o",col="coral")
legend("bottomright",c("m=0, sd=1","m=0, sd=5","m=5, sd=1","m=10,sd=1", "m=2,sd=10"), col=c("red","green","blue","yellow","coral"),pch=15)

#4
aproxPoisson <- function(n,p,k){
  lambda <- n * p
  sum <- 0
  for(x in 0:k)
    sum <- sum + exp(-lambda) * ((lambda ^ x)/factorial(x))
  return(sum)
}
aproxNormala <- function(n,p,k) {    
  numarator = k -n*p
  numitor = sqrt(n*p*(1-p))
  return(pnorm(numarator/numitor))
}
aproxNormalaFactor <- function(n,p,k){
  numarator = k +0.5 -n*p
  numitor = sqrt(n*p*(1-p))
  return(pnorm(numarator/numitor))
}
aproxCampPaulson <- function(n,p,k){
  b <- 1/ (9*(k+1))
  a <- 1/ (9*(n-k))
  r <- ((k+1)*(1-p)) / (p*(n-k))
  c <- (1-b) * r^(1/3)
  miu <- 1 - a
  sigma <-sqrt( a + b*r^(2/3))
  return(pnorm((c-miu)/sigma))
}
Matrice <- function(n,p) {
  matrice <- matrix(nrow=10,ncol=6)
  for(k in 1:10){
    matrice[k,1] <- k
    matrice[k,2] <- pbinom(k, size = n, prob = p)   
    matrice[k,3] <- aproxPoisson(n,p,k)
    matrice[k,4] <- aproxNormala(n,p,k)
    matrice[k,5] <- aproxNormalaFactor(n,p,k)
    matrice[k,6] <- aproxCampPaulson(n,p,k)
  }
  return(matrice)
}

print(Matrice(25,0.05))
print(Matrice(25,0.1))
print(Matrice(50,0.05))
print(Matrice(50,0.1))
print(Matrice(100,0.05))
print(Matrice(100,0.1))
#5
Kolmogorov <- function(n,p)
{
  Max<- 0
  Max1 <- 0
  Max2 <- 0
  Max3 <- 0

  matrice <- matrix(nrow=10,ncol=6)
  for(k in 1:10){
    matrice[k,1] <- k
    matrice[k,2] <- pbinom(k, size = n, prob = p)   
    matrice[k,3] <- aproxPoisson(n,p,k)
    matrice[k,4] <- aproxNormala(n,p,k)
    matrice[k,5] <- aproxNormalaFactor(n,p,k)
    matrice[k,6] <- aproxCampPaulson(n,p,k)
  }
 
      for(k in 1:10)
      {
        if(abs(matrice[k,2]-matrice[k,3]) > Max)
           Max <- abs(matrice[k,2]-matrice[k,3])
        
          if(abs(matrice[k,4]-matrice[k,2]) > Max1)
            Max1 <- abs(matrice[k,2]-matrice[k,4])
          if(abs(matrice[k,5]-matrice[k,2]) > Max2)
            Max2 <- abs(matrice[k,2]-matrice[k,5])
          if(abs(matrice[k,6]-matrice[k,2]) > Max3)
            Max3 <- abs(matrice[k,2]-matrice[k,6])
      }

  
  return(list(
    MaxP=Max,
    MaxN=Max1,
    MaxNF=Max2,
    MaxCP=Max3
  ))
}

domeniuX <- c(0.01, 0.50)
domeniuY <- c(0, 0.4)

plot(c(),type = "o",xlim=domeniuX,ylim=domeniuY,main="Kolmogorov n = 25")
n <- 25
n1 <- 50
n2 <- 100
for(p in seq(0.01, 0.5, 0.01)) {
  maxime <- Kolmogorov(n,p)
  points(x=p, y=maxime$MaxP, xlim=domeniuX, ylim=domeniuY,col="red")
  points(x=p, y=maxime$MaxN, xlim=domeniuX, ylim=domeniuY,col="green")
  points(x=p, y=maxime$MaxNF, xlim=domeniuX, ylim=domeniuY,col="coral")
  points(x=p, y=maxime$MaxCP, xlim=domeniuX, ylim=domeniuY,col="magenta")
}
plot(c(),type = "o",xlim=domeniuX,ylim=domeniuY,main="Kolmogorov n = 50")

for(p in seq(0.01, 0.5, 0.01)) {
  maxime <- Kolmogorov(n1,p)
  points(x=p, y=maxime$MaxP, xlim=domeniuX, ylim=domeniuY,col="red")
  points(x=p, y=maxime$MaxN, xlim=domeniuX, ylim=domeniuY,col="green")
  points(x=p, y=maxime$MaxNF, xlim=domeniuX, ylim=domeniuY,col="coral")
  points(x=p, y=maxime$MaxCP, xlim=domeniuX, ylim=domeniuY,col="magenta")
}

plot(c(),type = "o",xlim=domeniuX,ylim=domeniuY,main="Kolmogorov n = 100")
for(p in seq(0.01, 0.5, 0.01)) {
  maxime <- Kolmogorov(n2,p)
  points(x=p, y=maxime$MaxP, xlim=domeniuX, ylim=domeniuY,col="red")
  points(x=p, y=maxime$MaxN, xlim=domeniuX, ylim=domeniuY,col="green")
  points(x=p, y=maxime$MaxNF, xlim=domeniuX, ylim=domeniuY,col="coral")
  points(x=p, y=maxime$MaxCP, xlim=domeniuX, ylim=domeniuY,col="magenta")
}
#6

functie <- function(miu,sigma,lambda)
{
  return(function(x){
    result <- 2/sigma * dnorm((x-miu)/sigma) * pnorm(lambda*(x-miu)/sigma)
    return(result)
  })
}

domeniuX <- c(0, 10)
domeniuY <- c(0, 1)

plot(functie(4,5,6),type = "o",xlim=domeniuX,ylim=domeniuY,col="red",main="densitatea repartitiei normal-asimetrice")
plot(functie(1,2,3),type = "o",xlim=domeniuX,ylim=domeniuY,col="green",add=TRUE)
plot(functie(2,5,8),type = "o",xlim=domeniuX,ylim=domeniuY,col="blue",add=TRUE)
plot(functie(1,5,6),type = "o",xlim=domeniuX,ylim=domeniuY,col="yellow",add=TRUE)
plot(functie(4,8,6),type = "o",xlim=domeniuX,ylim=domeniuY,col="coral",add=TRUE)
legend(x=4, y=domeniuY[2],c("l=4,5,6","l=1,2,3","l=2,5,8","l=1,5,6", "l=4,8,6"), col=c("red","green","blue","yellow","coral"),pch=1,cex=0.8)



# 7
lambdaResult <-function(n,p){
  return(function(lambda){
    numarator <- (1 - (2 / 3.14) * (lambda^2/(1+lambda^2)))^3 * (1-2*p)^2
    numitor <- (2/3.14) * (4/3.14 - 1)^2 * (lambda^2/ (1+lambda^2))^3 * (n*p*(1-p))/(1-2*p)^2
    return(numarator-numitor)
  })
}

functieCareCalculeazaChestii <- function(n, p) {
  lambdaPatrat <- uniroot(lambdaResult(n,p),c(0,1000))$root
  lambda <- sign(1-p) * sqrt(lambdaPatrat)
  numarator <- n*p*(1-p)
  numitor <- 1 - (2/3.14) * (lambda/(1+lambda)) 
  sigma <- sqrt(numarator/numitor)
  miu <- n*p - sigma * sqrt((2*3.14) * (lambda/(1+lambda)))
  
  return(list(
    miu=miu,
    sigma=sigma,
    lambda=lambda
  ))
}

sigmaFunction <- function(n,p){
  lambdaPatrat <- uniroot(lambdaResult(n,p),c(0,1000))$root
  lambda <- sign(1-p) * sqrt(lambdaPatrat)
  numarator <- n*p*(1-p)
  numitor <- 1 - (2/3.14) * (lambda/(1+lambda)) 
  sigma <- sqrt(numarator/numitor)
  return(sigma)
}

miuFunction <- function(n,p){
  sigma <- sigmaFunction(n,p)
  lambdaPatrat <- uniroot(lambdaResult(n,p),c(0,1000))$root
  lambda <- sign(1-p) * sqrt(lambdaPatrat)
  miu <- n*p - sigma * sqrt((2*3.14) * (lambda/(1+lambda)))
  return(miu)
}


integ <- Vectorize(function(k,n,p){
  lambdaPatrat <- uniroot(lambdaResult(n,p),c(0,1000))$root
  lambda <- sign(1-p) * sqrt(lambdaPatrat)
  sigma <- sigmaFunction(n,p)
  miu <- miuFunction(n,p)
  return(integrate(function(t) { 2*dnorm(t)*pnorm(lambda * t) },-Inf,(k+0.5 - miu)/ sigma)$value)
})

domeniuX <- c(0, 10)
domeniuY <- c(0, 0.4)
n <- 25
p <- 0.1
p2 <- 0.25
chestii <- functieCareCalculeazaChestii(n, p)
barplot(dbinom(2:10, 25, 0.1), xlim=domeniuX,ylim=domeniuY,col="red")
plot(functie(chestii$miu, chestii$sigma, chestii$lambda), xlim=domeniuX,ylim=domeniuY, col="green",add=TRUE)
chestii2 <- functieCareCalculeazaChestii(n, p2)
barplot(dbinom(1:10, 25, 0.25), xlim=domeniuX,ylim=domeniuY,col="red")
plot(functie(chestii2$miu, chestii2$sigma, chestii2$lambda), xlim=domeniuX,ylim=domeniuY, col="green",add=TRUE)

#8

Kolmogorov_8 <- function(n,p)
{
   Max<- 0
 

  
  matrice <- matrix(nrow=10,ncol=3)
  for(k in 1:10){
    matrice[k,1] <- k
    matrice[k,2] <- pbinom(k, size = n, prob = p)   
    matrice[k,3] <- integ(k,n,p)
  
  }
  
  for(k in 1:10)
  {
    if(abs(matrice[k,2]-matrice[k,3]) > Max)
      Max <- abs(matrice[k,2]-matrice[k,3])
  }
  return(Max)
}
domeniuX <- c(0.05, 0.1)
domeniuY <- c(0, 2)

plot(c(),type = "o",xlim=domeniuX,ylim=domeniuY)
n <- 25
n1 <- 50
n2 <- 100
for(p in seq(0.05, 0.1, 0.001)) {
  points(x=p, y=Kolmogorov_8(n,p), xlim=domeniuX, ylim=domeniuY,col="magenta",main="Kolmogorov")
}
plot(c(),type = "o",xlim=domeniuX,ylim=domeniuY)

for(p in seq(0.05, 0.1, 0.001)) {

  points(x=p,y=Kolmogorov_8(n1,p), xlim=domeniuX, ylim=domeniuY,col="magenta",main="Kolmogorov")
}

plot(c(),type = "o",xlim=domeniuX,ylim=domeniuY)
for(p in seq(0.05, 0.1, 0.001)) {
 
  points(x=p, y=Kolmogorov_8(n2,p), xlim=domeniuX, ylim=domeniuY,col="magenta",main="Kolmogorov")
}


#Problema 2

#Punctul a), calcularea functiei fgam

#integram pentru a < 1
f1 <- function(a)
{
  f2 <- function(x)
  {
    x^(a-1) * exp(-x) 
  }
  integrate(f2, 0, Inf)$value
}

#fgam
fgam <- function(a)
{
  if(a == 1) 1 else
    if(a > 1 && floor(a) == a) (a-1)*fgam(a-1) else f1(a)
}



#Punctul b), calcularea functiei fbet

fbet <- function(a, b)
{
  (fgam(a)*fgam(b))/fgam(a+b)
}




#Punctul c)

fprobgamma1 <- function(a, b)
{
  #functia f dupa distributia gamma
  fgamma <- function(x)
  {
    (1/(b^a * fgam(a))) * x^(a-1) * exp(-x/b)
  } 
  
  integrate(fgamma, 0, 3)$value
  
}

fprobgamma2 <- function(a, b)
{
  #functia f dupa distributia gamma
  fgamma <- function(x)
  {
    (1/(b^a * fgam(a))) * x^(a-1) * exp(-x/b)
  }
  
  nr1 <- integrate(fgamma, 0, 5)
  nr2 <- integrate(fgamma, 0, 2)
  nr1$value - nr2$value
  
}

fprobgamma3 <- function(a, b)
{
  #functia f dupa distributia gamma
  fgamma <- function(x)
  {
    (1/(b^a * fgam(a))) * x^(a-1) * exp(-x/b)
  }
  
  
  (integrate(fgamma, 0, 4)$value - integrate(fgamma, 0, 3)$value) / (1-integrate(fgamma, 0, 2)$value)
  
}

fprobbeta4 <- function(a, b)
{
  #functia f dupa distributia beta
  fbeta <- function(x)
  {
    (1/fbet(a, b)) * x^(a-1) * (1-x)^(b-1)
  }
  
  1 - integrate(fbeta, 0, 2)$value
  
}

fprobbeta5 <- function(a, b)
{
  #functia f dupa distributia beta
  fbeta <- function(x)
  {
    (1/fbet(a, b)) * x^(a-1) * (1-x)^(b-1)
  }
  
  integrate(fbeta, 0, 6)$value - integrate(fbeta, 0, 4)$value
  
}

fprobbeta6 <- function(a, b)
{
  #functia f dupa distributia beta
  fbeta <- function(x)
  {
    (1/fbet(a, b)) * x^(a-1) * (1-x)^(b-1)
  }
  
  (integrate(fbeta, 0, 1)$value - integrate(fbeta, 0, 0)$value) / integrate(fbeta, 0, 7)$value
  
}

v1 <- fprobgamma1(1,2)
v2 <- fprobgamma2(1,0.5)
v3 <- fprobgamma3(1,5)
v4 <- fprobbeta4(5,10)
v5 <- fprobbeta5(1,1)
v6 <- fprobbeta6(2,3)
print(v1)
print(v2)
print(v3)
print(v4)
print(v5)
print(v6)




#Punctul d)   

fprobgamma11 <- function(a, b)
{
  #functia f dupa distributia gamma
  fgamma <- function(x)
  {
    (1/(b^a * gamma(a))) * x^(a-1) * exp(-x/b)
  } 
  
  #integrate(fgamma, 0, 3)$value
  return(pgamma(3, shape=a, scale=b))
  
}

fprobgamma22 <- function(a, b)
{
  #functia f dupa distributia gamma
  fgamma <- function(x)
  {
    (1/(b^a * gamma(a))) * x^(a-1) * exp(-x/b)
  }
  
  #nr1 <- integrate(fgamma, 0, 5)
  #nr2 <- integrate(fgamma, 0, 2)
  #nr1$value - nr2$value
  return(pgamma(5, a, b) - pgamma(2, a, b))
  
}

fprobgamma33 <- function(a, b)
{
  #functia f dupa distributia gamma
  fgamma <- function(x)
  {
    (1/(b^a * gamma(a))) * x^(a-1) * exp(-x/b)
  }
  
  
  #(integrate(fgamma, 0, 4)$value - integrate(fgamma, 0, 3)$value) / (1-integrate(fgamma, 0, 2)$value)
  return((pgamma(4, a, b) - pgamma(3, a, b))/(1 - pgamma(3, a, b)))
}

fprobbeta44 <- function(a, b)
{
  #functia f dupa distributia beta
  fbeta <- function(x)
  {
    (1/beta(a, b)) * x^(a-1) * (1-x)^(b-1)
  }
  
  #1 - integrate(fbeta, 0, 2)$value
  return(1 - pbeta(2, shape1=a, shape2=b))
  
}

fprobbeta55 <- function(a, b)
{
  #functia f dupa distributia beta
  fbeta <- function(x)
  {
    (1/beta(a, b)) * x^(a-1) * (1-x)^(b-1)
  }
  
  #integrate(fbeta, 0, 6)$value - integrate(fbeta, 0, 4)$value
  return(pbeta(6, shape1=a, shape2=b) - pbeta(4, shape1=a, shape2=b))
}

fprobbeta66 <- function(a, b)
{
  #functia f dupa distributia beta
  fbeta <- function(x)
  {
    (1/beta(a, b)) * x^(a-1) * (1-x)^(b-1)
  }
  
  #(integrate(fbeta, 0, 1)$value - integrate(fbeta, 0, 0)$value) / integrate(fbeta, 0, 7)$value
  return((pbeta(1, shape1=a, shape2=b) - pbeta(0, shape1=a, shape2=b)) / pbeta(7, shape1=a, shape2=b))
}

v11 <- fprobgamma11(1,2)
v22 <- fprobgamma22(1,0.5)
v33 <- fprobgamma33(1,5)
v44 <- fprobbeta44(5,10)
v55 <- fprobbeta55(1,1)
v66 <- fprobbeta66(2,3)
print(v11)
print(v22)
print(v33)
print(v44)
print(v55)
print(v66)

vect1 <- c(v1, v2, v3, v4, v5, v6)
vect2 <- c(v11, v22, v33, v44, v55, v66)
mat <- cbind(vect1, vect2)
print(mat)


#3
#a

solve_element <- function(A, r, c, R, C) {
  
  no_elem_line = 0
  sum_line = 0
  no_elem_col = 0
  sum_col = 0
  for (col in 2:(C - 1)) {
    if (A[r, col] != -1) {
      no_elem_line = 1 + no_elem_line
      sum_line = sum_line + A[r, col]
    }
  }
  if (no_elem_line == C - 3) {
    ans = A[r, C] - sum_line
    return(ans)
  }
  for (row in 2:(R - 1)) {
    if (A[row, c] != -1) {
      no_elem_col = 1 + no_elem_col
      sum_col = sum_col + A[row, c]
    }
  }
  if (no_elem_col == R - 3) {
    ans = A[R, c] - sum_col
    return(ans)
  }
  ans = -1
  cov_x
  return(ans)
}


generate_random_joint_distribution <- function(n, m)
{
  x = c(sample(1:(10 * n), n, rep=FALSE))
  y = c(sample(1:(10 * m), m, rep=FALSE))
  x = sort(x)
  y = sort(y)
  aux = matrix(sample(1:(n*m), n * m, replace = TRUE), nrow = n, ncol = m)
  joint_distribution = matrix(0, nrow = n + 2, ncol = m + 2)
  for (i in 1:n) {
    joint_distribution[i + 1, 1] = x[i]
    joint_distribution[i + 1, m + 2] = rowSums(aux)[i]
  }
  for (i in 1:m) {
    joint_distribution[1, i + 1] = y[i]
    joint_distribution[n + 2, i + 1] = colSums(aux)[i]
  }
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      joint_distribution[i, j] = aux[i - 1,j - 1]
    }
  }
  for (i in 2:(n + 2)) {
    for (j in 2:(m + 2)) {
      joint_distribution[i, j] = joint_distribution[i, j] / sum(aux)
    }
  }
  colnames_joint_distribution = c(1:m)
  rownames_joint_distribution = c(1:n)
  for (i in 1:m) {
    colnames_joint_distribution[i] = paste(c("y", colnames_joint_distribution[i]), collapse = "", sep = "")
  }
  for (i in 1:n) {
    rownames_joint_distribution[i] = paste(c("x", rownames_joint_distribution[i]), collapse = "", sep = "")
  }
  rownames(joint_distribution) = c("Y", rownames_joint_distribution, "q")
  colnames(joint_distribution) = c("X", colnames_joint_distribution, "p")
  return(joint_distribution)
}

erase_some_values <- function(joint_distribution) 
{
  
  joint_distribution[2, 2] = -1
  joint_distribution[2, 3] = -1
  mini = min(nrow(joint_distribution), ncol(joint_distribution)) - 1
  for (p in 3:mini) {
    joint_distribution[p, p] = -1
  }
  return(joint_distribution)
}

frepcomgen <- function(n, m)
{
  joint_distribution = generate_random_joint_distribution(n, m)
  print("Repartitia comuna completa")
  print(joint_distribution)
  partial_joint_distribution = erase_some_values(joint_distribution = joint_distribution)
  print("Repartitia comuna incompleta")
  print(partial_joint_distribution)
  return(partial_joint_distribution)
}
partial_joint_distribution = frepcomgen(2,3)
#media variabilei x cu distributia p
get_mean <- function(x, p)
{
  mean_x = 0
  for(i in 1:dim(array(x)))
    mean_x = mean_x + x[i] * p[i]
  return(mean_x)
}

#cov[X, Y] unde X are distributia p, y are distributia q, iar repcom este repartitia comuna
get_cov <- function(x, p, y, q, repcom)
{
  mean_x = get_mean(x, p)
  mean_y = get_mean(y, q)
  new_x = x - mean_x
  new_y = y - mean_y
  xy <- 1:(dim(array(x)) + dim(array(y)))
  pxy <- 1:(dim(array(x)) + dim(array(y)))
  k = 1
  for(i in 1:dim(array(x)))
  {
    for(j in 1:dim(array(y)))
    {
      xy[k] <- new_x[i] * new_y[j]
      pxy[k] <- repcom[i, j]
      k = k + 1
    }
  }
  #cov[X, Y] = E[(X - E[X]) * (Y - E[Y])]
  return(get_mean(xy, pxy))
}

# verifica independenta folosind repartitia comuna
fverind <- function(repcom, p, q)
{
  for(i in 1:dim(repcom)[1])
  {
    for(j in 1:dim(repcom)[2])
      if(repcom[i, j] != p[i] * q[j])
        return(FALSE)
  }
  return(TRUE)
}

# pentru necorelare este de ajuns sa verificam covarianta
fvernecor <- function(x, p, y, q, repcom)
{
  if(get_cov(x, p, y, q, repcom) == 0)
    return(TRUE)
  return(FALSE)
}

complete_prob <- function(p)
{
  for(i in 1:dim(array(p)))
  {
    # gaseste valoarea lipsa din repartitie si completeaza astfel incat suma sa fie 1
    if (is.na(p[i]))
    {
      p[i] <- 0
      p[i] <- 1 - sum(p)
    }
  }
}

fcomplrepcom <- function(x, p, y, q, repcom)
{
  n = dim(array(x))
  m = dim(array(y))
  found = TRUE
  # repeta cat timp inca au fost gasite pozitii ce au putut fi completate
  while(found)
  {
    found = FALSE
    # caut linii in repcom in care lipseste doar o valoare
    for(i in 1:n)
    {
      # daca am toate valorile, dar lipseste p[i] se poate calcula
      if(sum(is.na(repcom[i,])) == 0)
      {
        if(is.na(p[i]))
        {
          p[i] <- sum(repcom[i,])
          found = TRUE
        }
      }
      # daca lipseste o singura valoare, dar avem p[i] se poate calcula valoarea lipsa
      if(sum(is.na(repcom[i,])) == 1)
      {
        if(!is.na(p[i]))
        {
          for(j in 1:m)
          {
            if(is.na(repcom[i, j]))
            {
              repcom[i, j] = 0
              repcom[i, j] = p[i] - sum(repcom[i, ])
            }
          }
          found = TRUE
        }
      }
    }
    
    # repeta pentru coloane
    for(j in 1:m)
    {
      if(sum(is.na(repcom[, j])) == 0)
      {
        if(is.na(q[j]))
        {
          q[j] <- sum(repcom[, j])
          found = TRUE
        }
      }
      if(sum(is.na(repcom[, j])) == 1)
      {
        if(!is.na(q[j]))
        {
          for(i in 1:n)
          {
            if(is.na(repcom[i, j]))
            {
              repcom[i, j] = 0
              repcom[i, j] = q[j] - sum(repcom[, j])
            }
          }
          found = TRUE
        }
      }
    }
    
    # daca lipseste o singura valoare din repartitia lui X
    if(sum(is.na(p)) == 1)
    {
      complete_prob(p)
      found = TRUE
    }
    
    # asemenea pentru y
    if(sum(is.na(q)) == 1)
    {
      complete_prob(q)
      found = TRUE
    }
  }
  return(list(a=repcom, b=p, c=q))
}

# generam un caz de 2x3 in care lipsesc valori din repartitia comuna
x1 <- c(0, 1)
p1 <- c(NA, 3/9)

y1 <- c(0, 1, 2)
q1 <- c(1/9, NA, 1/9)

repcom1 <- matrix(c(NA, 0, 4/9, NA, 1/9, 0), nrow=2, ncol=3)
print(repcom1)
# completam repartitia comuna
res <- fcomplrepcom(x1, p1, y1, q1, repcom1)
repcom1 = res$a
p1 = res$b
q1 = res$c
print(repcom1)

#varificam independenta si necorelarea
print(fverind(repcom1, p1, q1))
print(fvernecor(x1, p1, y1, q1, repcom1))

