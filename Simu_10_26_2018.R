##Part 1 for homogeneous Markov chain with 3 states {1,2,3}

## define a function to get n samples from a given discrete distribution p 
sampleDist = function(n,p) { 
  sample(x = c(1,2,3), n, replace = T, prob = p) 
}
## set the predefine transition matrix
p1=c(0.5,0.25,0.25)
p2=c(1/3,1/3,1/3)
p3=c(0.25,0.5,0.25)

## define a function to record the length of longest consecutive visits of state 1 
## with time length lent
LCV2 <- function(lent,p1,p2,p3){
  currentpos=1
  ccount=0
  hcount=0
  for (i in (1:lent)){
    if (currentpos==1){
      currentpos=sampleDist(1,p1)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else if (currentpos==2){
      currentpos=sampleDist(1,p2)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else{
      currentpos=sampleDist(1,p3)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
      
    }
  }
  return (hcount)
}

##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:5000)){
  tr=c(tr,LCV2(10000,p1,p2,p3))
}
mean(tr)


##derive the mean of L(1,100000) for 500 runs
tr=c()
for (i in (1:500)){
  tr=c(tr,LCV2(100000,p1,p2,p3))
}
mean(tr)

## a function to derive the sum of p_11(l) from l=1 to lent
findsump <- function(lent){
sump=0
init=matrix(c(1,0,0),ncol=3,nrow=1)
tm1=t(matrix(c(0.5,0.25,0.25,1/3,1/3,1/3,0.25,0.5,0.25), ncol = 3, nrow = 3))
tm=t(matrix(c(0.5,0.25,0.25,1/3,1/3,1/3,0.25,0.5,0.25), ncol = 3, nrow = 3))
sump=c(init%*%tm)[1]
for (i in (2:lent)){
  tm=tm%*%tm1
  sump=sump+c(init%*%tm)[1]
}
return (sump)
}

##find the quantity of sum of p_11(l) from l=1 to lent for lent=10000 and lent=100000
findsump(10000)
findsump(100000)


##Pat 2 Simulation for time nonhomogeneous Markov chain with 3 states {1,2,3}
## generate all random number needed for P(t) from t=1 to 100000
set.seed(100)
RG=sample(1:30,8*100000,replace=TRUE)
## define a function to generate the discrete distribution for each time t from sequence RG
getdist<-function(RG,i,s){
  if (s==1){
    x=RG[(8*(i-1)+1):(8*(i-1)+2)]
    a=sum(x)
    x=c(a,x)
    p=x/sum(x)
  }
  else if (s==2){
    x=RG[(8*(i-1)+3):(8*(i-1)+5)]
    p=x/sum(x)
  }
  else {
    x=RG[(8*(i-1)+6):(8*(i-1)+8)]
    p=x/sum(x)
  }
  return (p)
}

LCV3 <- function(lent,RG){
  currentpos=1
  ccount=0
  hcount=0
  for (i in (1:lent)){
    if (currentpos==1){
      p1=getdist(RG,i,1)
      currentpos=sampleDist(1,p1)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else if (currentpos==2){
      p2=getdist(RG,i,2)
      currentpos=sampleDist(1,p2)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else{
      p3=getdist(RG,i,3)
      currentpos=sampleDist(1,p3)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
      
    }
  }
  return (hcount)
}

## a function to derive the sum of p_11(l) from l=1 to lent
findsump <- function(lent,RG){
  sump=0
  init=matrix(c(1,0,0),ncol=3,nrow=1)
  tm=t(matrix(c(getdist(RG,1,1),getdist(RG,1,2),getdist(RG,1,3)), ncol = 3, nrow = 3))
  sump=c(init%*%tm)[1]
  for (i in (2:lent)){
    tm=tm%*%t(matrix(c(getdist(RG,i,1),getdist(RG,i,2),getdist(RG,i,3)), ncol = 3, nrow = 3))
    sump=sump+c(init%*%tm)[1]
  }
  return (sump)
}
findsump(10000,RG)
findsump(100000,RG)

##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:5000)){
  tr=c(tr,LCV3(10000,RG))
}
mean(tr)


##derive the mean of L(1,100000) for 500 runs
tr=c()
for (i in (1:500)){
  tr=c(tr,LCV3(100000,RG))
}
mean(tr)


##part 3 a nonhomogeneous Markov chain with p_11(l)=1/(l+1)

## define a function to generate the desired transition matrix at time t to ganrantee that p_11(t)=(1/(t+1))^0.1
simulateTransition<-function(pi,t){
  if (t==1){
    p1=runif(1,0,1/2)
    p2=1/2-p1
    p3=runif(1,0,1)
    p4=runif(1,0,1-p3)
    p5=1-p3-p4
    p6=runif(1,0,1)
    p7=runif(1,0,1-p6)
    p8=1-p6-p7
    pt=c(1/2,p1,p2,p3,p4,p5,p6,p7,p8)
  }
  else{
    x0=pi[1]
    x1=pi[2]
    x2=pi[3]
    p1=runif(1,0,1/2)
    p2=1/2-p1
    p3=min(((1/(t+1))-(0.5)*x0)/(x1+x2),0.99)
    p4=runif(1,0,1-p3)
    p5=1-p3-p4
    p6=min(((1/(t+1))-(0.5)*x0)/(x1+x2),0.99)
    p7=runif(1,0,1-p6)
    p8=1-p6-p7
    pt=c(1/2,p1,p2,p3,p4,p5,p6,p7,p8)
  }
  return (pt)
}

##define a function to generate all transition matrix from t=1 to lent
generateAll<-function(lent){
  pi=c(1,0,0)
  alltransi=c()
  sump=0
  for (i in (1:lent)){
    pp=simulateTransition(pi,i)
    pi=pi%*%t(matrix(pp, ncol = 3, nrow = 3))
    sump=sump+c(pi)[1]
    alltransi=c(alltransi,pp)
  }
  return (c(sump,alltransi))
}

LCV4<- function(lent,GA){
  GA=GA[2:length(GA)]
  currentpos=1
  ccount=0
  hcount=0
  for (i in (1:lent)){
    if (currentpos==1){
      p1=GA[((i-1)*9+1):((i-1)*9+3)]
      currentpos=sampleDist(1,p1)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else if (currentpos==2){
      p2=GA[((i-1)*9+4):((i-1)*9+6)]
      currentpos=sampleDist(1,p2)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else{
      p3=GA[((i-1)*9+7):((i-1)*9+9)]
      currentpos=sampleDist(1,p3)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
      
    }
  }
  return (hcount)
}
set.seed(100)
GA=generateAll(10000)
##derive the mean of L(1,10000) for 1500 runs
tr=c()
for (i in (1:1500)){
  tr=c(tr,LCV4(10000,GA))
}
GA[1]
mean(tr)

hist(tr, 
     main="Histogram for L(1,10000) with 1500 times", 
     xlab="L(1,10000)", 
     border="blue", 
     col="green", 
     xlim=c(min(tr),max(tr)), 
     breaks=max(tr)-min(tr), 
     prob = TRUE,xaxt="n")
## draw the x-axis with user-defined tick-marks
axis(side=1, at=seq (min(tr),max(tr),1))

##derive the mean of L(1,50000) for 1500 runs
set.seed(100)
GA=generateAll(50000)
tr=c()
for (i in (1:1500)){
  tr=c(tr,LCV4(50000,GA))
}
mean(tr)
