##Simulation for theorem 2.1 and 2.2
sampleDist = function(n,p) { 
  sample(x = c(1,2,3), n, replace = T, prob = p) 
}
## set the predefine transition matrix
p1=c(0.5,0.25,0.25)
p2=c(1/3,1/3,1/3)
p3=c(0.25,0.5,0.25)
## Theorem 2.1
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
powA = function(a,n)
{
  if (n==1)  return (a)
  if (n==2)  return (a%*%a)
  if (n>2) return ( a%*%powA(a,n-1))
}

getfii<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/2
  pp1=t(c(1/4,1/4))
  pp2=rbind(c(1/3,1/3),c(1/2,1/4))
  pp3=c(1/3,1/4)
  s[2]=pp1%*%pp3
  for (i in (3:lent)){
    s[i]=pp1%*%(powA(pp2,i-2))%*%pp3
  }
  return (s)
}

getpii<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/2
  mm=rbind(p1,p2,p3)
 
  for (i in (2:lent)){
    s[i]=powA(mm,i)[1,1]
  }
  return (s)
}
getpiik<-function(fi,pi,lent,r,k){
  s=pi
  s[k]=s[k]-r^k
  for (i in ((k+1):lent)){
    s[i]=s[i]-pi[i-k]*(r^k)
    sum1=0
    for (h in (2:(i-k+1))){
      for (j in (2:h)){
        if (j==h){
          if (i-k-h+1==0){
            sum1=sum1+fi[j]*(r^(k-1))
          }
          else{
            sum1=sum1+fi[j]*pi[i-k-h+1]*(r^(k-1)) 
          }
        }
        else{
          if (i-k-h+1==0){
            sum1=sum1+s[h-j]*fi[j]*(r^(k-1))
          }
          else{
          sum1=sum1+s[h-j]*fi[j]*pi[i-k-h+1]*(r^(k-1)) 
          }
          }
      }
    }
   s[i]=s[i]-sum1
  }
  return (s)
}
getP<-function(fi,pi,piik,lent,r,k){
  s=pi
  s[k]=1-r^k
  for (n in ((k+1):lent)){
    s[n]=1-r^k
    sum1=0
    for (l in (2:(n-k+1))){
      for (h in (0:(l-2))){
        if (h==0){
          sum1=sum1+fi[l-h]*(r^(k-1))
        }
        else{
          sum1=sum1+piik[h]*fi[l-h]*(r^(k-1)) 
        }
      }
    }
    s[n]=s[n]-sum1
  }
  return (s)
  }

fii=getfii(20,p1,p2,p3)
pii=getpii(20,p1,p2,p3)
piik=getpiik(fii,pii,20,0.5,5)
P=getP(fii,pii,piik,20,0.5,5)


##derive the mean of L(1,10000) for 5000 runs
h=pii
set.seed(100)
for (j in (5:20)){
tr=c()
for (i in (1:150000)){
  tr=c(tr,LCV2(j,p1,p2,p3))
}
s1=tr<5
h[j]=sum(s1)/150000
}
P
h

###Theorem 2.2
LCV22 <- function(lent,p1,p2,p3){
  currentpos=1
  ccount=0
  hcount=0
  for (i in (1:lent)){
    if (currentpos==1){
      currentpos=sampleDist(1,p1)
      if (currentpos==2){
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
      if (currentpos==2){
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
      if (currentpos==2){
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

getfjj<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/3
  pp1=t(c(1/3,1/3))
  pp2=rbind(c(1/2,1/4),c(1/4,1/4))
  pp3=c(1/4,1/2)
  s[2]=pp1%*%pp3
  for (i in (3:lent)){
    s[i]=pp1%*%(powA(pp2,i-2))%*%pp3
  }
  return (s)
}


getfij<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/4
  pp1=t(c(1/2,1/4))
  pp2=rbind(c(1/2,1/4),c(1/4,1/4))
  pp3=c(1/4,1/2)
  s[2]=pp1%*%pp3
  for (i in (3:lent)){
    s[i]=pp1%*%(powA(pp2,i-2))%*%pp3
  }
  return (s)
}

getpij<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/4
  mm=rbind(p1,p2,p3)
  
  for (i in (2:lent)){
    s[i]=powA(mm,i)[1,2]
  }
  return (s)
}


getpjj<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/3
  mm=rbind(p1,p2,p3)
  
  for (i in (2:lent)){
    s[i]=powA(mm,i)[2,2]
  }
  return (s)
}


getpijk<-function(fij,fjj,pij,pjj,lent,r,k){
  s=pij
  s[k]=s[k]-(r^(k-1))*pij[1]
  s[k+1]=s[k+1]-(r^k)*pij[1]-(r^(k-1))*fij[2]
  for (i in ((k+2):lent)){
    sum1=0
    for (t in (3:(i-k+1))){
      for (z in (1:(t-2))){
        if (i-k-t+1==0){
          if (z<k){
            sum1=sum1+pij[z]*fjj[t-z]*r^(k-1)
          }
          else{
            sum1=sum1+s[z]*fjj[t-z]*r^(k-1) 
          }
        }
        else{
          if (z<k){
            sum1=sum1+pij[z]*fjj[t-z]*pjj[i-k-t+1]*r^(k-1) 
          }
          else{
          sum1=sum1+s[z]*fjj[t-z]*pjj[i-k-t+1]*r^(k-1) 
          }
        }
      }
    }
    for (w in (1:(i-k+1))){
      if (i-w-k+1==0){
        sum1=sum1+fij[w]*(r^(k-1))
      }
      else{
      sum1=sum1+fij[w]*pjj[i-w-k+1]*(r^(k-1))
      }
    }
    s[i]=s[i]-sum1
  }
  return (s)
}

getP<-function(fij,fjj,pijk,lent,r,k){
  s=pijk
  s[k]=1-fij[1]*(r^(k-1))
  s[k+1]=1-fij[1]*(r^(k-1))-fij[2]*(r^(k-1))
  for (n in ((k+2):lent)){
    sum1=0
    for (l in (3:(n-k+1))){
      for (t in (1:(l-2))){
        hold1=pijk[t]*fjj[l-t]*(r^(k-1))
        sum1=sum1+hold1
      }
    }
    s[n]=1-sum1
    sum2=0
    for (z in (1:(n-k+1))){
      hold=fij[z]*(r^(k-1))
      sum2=sum2+hold
    }
    s[n]=s[n]-sum2
  }
  return (s)
}



fjj=getfjj(20,p1,p2,p3)
fij=getfij(20,p1,p2,p3)
pij=getpij(20,p1,p2,p3)
pjj=getpjj(20,p1,p2,p3)
pijk=getpijk(fij,fjj,pij,pjj,20,1/3,5)
P=getP(fij,fjj,pijk,20,1/3,5)


##derive the mean of L(1,10000) for 5000 runs
h=pjj
set.seed(100)
for (j in (5:20)){
  tr=c()
  for (i in (1:150000)){
    tr=c(tr,LCV22(j,p1,p2,p3))
  }
  s1=tr<5
  h[j]=sum(s1)/150000
}
P
h






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