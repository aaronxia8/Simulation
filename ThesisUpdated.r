library(ggplot2)
##simulation for bernoulli case
## a function to derive the array of longest head for each time
B1<-function(num,p){
  c1=rbinom(num, 1,p)
  curmax=0
  hismax=0
  for (i in (1:num)){
    if (c1[i]==1){
      curmax=curmax+1
      hismax=max(curmax,hismax)
    }
    else{
      curmax=0
    }
    c1[i]=hismax
  }
  return (c1)
}
y1=B1(1000,0.5)
x1=seq(3,1000)
for (i in (1:9999)){
  ss=B1(1000,0.5)
  y1=c(y1,ss)
}
s1=t(matrix(y1,ncol=10000))
meany=colMeans(s1,dims=1)
s2=c(log2(1))
s4=c(log2(1)+(-digamma(1)/log(2))-(3/2))
s3=c(log2(0.5))
for (i in (1:999)){
  s2=c(s2,log2(i))
  s3=c(s3,log2(0.5*i))
  s4=c(s4,log2(i)+(-digamma(1)/log(2))-(3/2))
}

##plot figure 1 in thesis (comparsion of 3 estimates of E(L(n)) and simulations)
df1 <- data.frame("Length" = c(x1,x1,x1,x1), "Value" = c(meany[3:1000],s2[3:1000],s3[3:1000],s4[3:1000]),"Type"=c(rep("Simulation",998),rep("log2(n)",998),rep("log2(0.5n)",998),rep("E(L(n)) by Flajolet et al.",998)))
d <- ggplot(df1, aes(Length, Value))
d  + geom_point(aes(colour = Type),size=2,alpha=0.8,position=position_jitter(h=0.05, w=0.05)) +  xlab("Length n") + ylab("Estimation of E(L(n))")+ theme_grey(base_size = 30)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(size=25))













#########Simulation for theorem 2.1 and 2.2
#########

sampleDist = function(n,p) { 
  sample(x = c(1,2,3), n, replace = T, prob = p) 
}

## set the predefined transition matrix to derive the sequence f_jj, p_ij,f_ij,p_jj
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
##a helper function to do matrix multiplication
powA = function(a,n)
{
  if (n==1)  return (a)
  if (n==2)  return (a%*%a)
  if (n>2) return ( a%*%powA(a,n-1))
}

##a function to derive f_ii
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

##a function to derive p_ii
getpii<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/2
  mm=rbind(p1,p2,p3)
  
  for (i in (2:lent)){
    s[i]=powA(mm,i)[1,1]
  }
  return (s)
}

##a function to derive p_iik
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

#a function to derive exact probability
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
## define a function to record the length of longest consecutive visits of state 1 
## with time length lent
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


##a function to find f_jj
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

##a function to find f_ij
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

##a function to find p_jj
getpjj<-function(lent,p1,p2,p3){
  s=seq(1,lent)
  s[1]=1/3
  mm=rbind(p1,p2,p3)
  
  for (i in (2:lent)){
    s[i]=powA(mm,i)[2,2]
  }
  return (s)
}

##a function to find p_ijk
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

##a function to derive exact probability
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









###Appendix A: Illustration of Lemma 2.6
##A.1 for homogeneous Markov chain with 3 states {1,2,3}
set.seed(100)
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
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))

##derive the mean of L(1,100000) for 500 runs
tr=c()
for (i in (1:1000)){
  tr=c(tr,LCV2(100000,p1,p2,p3))
}
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))

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






##plot
##A.1 for homogeneous Markov chain with 3 states {1,2,3}
set.seed(100)
## define a function to get n samples from a given discrete distribution p 
sampleDist = function(n,p) { 
  sample(x = c(1,2,3), n, replace = T, prob = p) 
}
## set the predefine transition matrix
p1=c(0.5,0.25,0.25)
p2=c(1/3,1/3,1/3)
p3=c(0.25,0.5,0.25)

## a function to derive the sum of p_11(l) from l=1 to lent

sump=c()
sump1=0
init=matrix(c(1,0,0),ncol=3,nrow=1)
tm1=t(matrix(c(0.5,0.25,0.25,1/3,1/3,1/3,0.25,0.5,0.25), ncol = 3, nrow = 3))
tm=t(matrix(c(0.5,0.25,0.25,1/3,1/3,1/3,0.25,0.5,0.25), ncol = 3, nrow = 3))
sump1=c(init%*%tm)[1]
sump=c(sump,log2(sump1))
for (i in (2:1000)){
  tm=tm%*%tm1
  sump1=sump1+c(init%*%tm)[1]
  sump=c(sump,log2(sump1))
}



## define a function to record the length of longest consecutive visits of state 1 
## with time length lent
LCV2 <- function(lent,p1,p2,p3){
  currentpos=1
  ccount=0
  hcount=0
  s1=c()
  for (i in (1:lent)){
    if (currentpos==1){
      currentpos=sampleDist(1,p1)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
      }
    }
    else if (currentpos==2){
      currentpos=sampleDist(1,p2)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
      }
    }
    else{
      currentpos=sampleDist(1,p3)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
      }
      
    }
  }
  return (s1)
}

##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:10000)){
  tr=c(tr,LCV2(1000,p1,p2,p3))
}
s1=t(matrix(tr,ncol=10000))
meany=colMeans(s1,dims=1)
x1=seq(3,1000)
##zz=c()
##for (i in (1:1000)){
##  zz=c(zz,log(i)/log(2))
##}
##plot figure 1 in thesis (comparsion of 3 estimates of E(L(n)) and simulations)
##df1 <- data.frame("Length" = c(x1,x1,x1), "Value" = c(sump[3:1000],meany[3:1000],zz[3:1000]),"Type"=c(rep("Our estimation",998),rep("Simulation",998),rep("log2(n)",998)))

df1 <- data.frame("Length" = c(x1,x1), "Value" = c(sump[3:1000],meany[3:1000]),"Type"=c(rep("Our estimation",998),rep("Simulation",998)))
d <- ggplot(df1, aes(Length, Value))
d  + geom_point(aes(colour = Type),size=3,alpha=0.8,position=position_jitter(h=0.005, w=0.005)) +  xlab("Length n") + ylab("Estimation of E(L(1,n))")+ theme_grey(base_size = 24)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(size=25))

























##A.2 Simulation for time nonhomogeneous Markov chain with 3 states {1,2,3}
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
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))

##derive the mean of L(1,100000) for 1500 runs
tr=c()
for (i in (1:1500)){
  tr=c(tr,LCV3(100000,RG))
}
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))



##plot
##A.2 Simulation for time nonhomogeneous Markov chain with 3 states {1,2,3}
## generate all random number needed for P(t) from t=1 to 1000
sampleDist = function(n,p) { 
  sample(x = c(1,2,3), n, replace = T, prob = p) 
}

set.seed(100)
RG=sample(1:30,8*1000,replace=TRUE)
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
## a function to derive the sum of p_11(l) from l=1 to lent
sump=C()
sump1=0
init=matrix(c(1,0,0),ncol=3,nrow=1)
tm=t(matrix(c(getdist(RG,1,1),getdist(RG,1,2),getdist(RG,1,3)), ncol = 3, nrow = 3))
sump=c(init%*%tm)[1]
for (i in (2:1000)){
  tm=tm%*%t(matrix(c(getdist(RG,i,1),getdist(RG,i,2),getdist(RG,i,3)), ncol = 3, nrow = 3))
  sump1=sump1+c(init%*%tm)[1]
  sump=c(sump,log2(sump1))
}


LCV3 <- function(lent,RG){
  currentpos=1
  ccount=0
  hcount=0
  s1=c()
  for (i in (1:lent)){
    if (currentpos==1){
      p1=getdist(RG,i,1)
      currentpos=sampleDist(1,p1)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
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
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
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
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
      }
      
    }
  }
  return (s1)
}


##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:10000)){
  tr=c(tr,LCV3(1000,RG))
}
s1=t(matrix(tr,ncol=10000))
meany=colMeans(s1,dims=1)
x1=seq(3,1000)
##generate log2(n) sequence
zz=c()
for (i in (1:1000)){
  zz=c(zz,log(i)/log(2))
}
##plot figure 1 in thesis (comparsion of 3 estimates of E(L(n)) and simulations)
##df1 <- data.frame("Length" = c(x1,x1,x1), "Value" = c(sump[3:1000],meany[3:1000],zz[3:1000]),"Type"=c(rep("Our estimation",998),rep("Simulation",998),rep("log2(n)",998)))
df1 <- data.frame("Length" = c(x1,x1), "Value" = c(sump[3:1000],meany[3:1000]),"Type"=c(rep("Our estimation",998),rep("Simulation",998)))
d <- ggplot(df1, aes(Length, Value))
d  + geom_point(aes(colour = Type),size=3,alpha=0.8,position=position_jitter(h=0.005, w=0.005)) +  xlab("Length n") + ylab("Estimation of E(L(1,n))")+ theme_grey(base_size = 24)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(size=25))
















##A.2.2 a nonhomogeneous Markov chain with p_11(l)=1/(l+1)

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
##derive the mean of L(1,10000) for 2500 runs
tr=c()
for (i in (1:5000)){
  tr=c(tr,LCV4(10000,GA))
}
GA[1]
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))


##derive the mean of L(1,100000) for 2500 runs
set.seed(100)
GA=generateAll(100000)
tr=c()
for (i in (1:2500)){
  tr=c(tr,LCV4(100000,GA))
}
GA[1]
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))


##plot
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

pi=c(1,0,0)
alltransi=c()
sum1=0
sump=c()
for (i in (1:1000)){
  pp=simulateTransition(pi,i)
  pi=pi%*%t(matrix(pp, ncol = 3, nrow = 3))
  sum1=sum1+c(pi)[1]
  sump=c(sump,log2(sum1))
  alltransi=c(alltransi,pp)
}


LCV4<- function(lent,GA){
  currentpos=1
  ccount=0
  hcount=0
  ss1=c()
  for (i in (1:lent)){
    if (currentpos==1){
      p1=GA[((i-1)*9+1):((i-1)*9+3)]
      currentpos=sampleDist(1,p1)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        ss1=c(ss1,hcount)
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
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        ss1=c(ss1,hcount)
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
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        ss1=c(ss1,hcount)
      }
      
    }
  }
  return (ss1)
}
set.seed(100)
##derive the mean of L(1,10000) for 2500 runs
tr=c()
for (i in (1:10000)){
  tr=c(tr,LCV4(1000,alltransi))
}
s1=t(matrix(tr,ncol=10000))
meany=colMeans(s1,dims=1)
x1=seq(3,1000)

##plot figure 1 in thesis (comparsion of 3 estimates of E(L(n)) and simulations)
df1 <- data.frame("Length" = c(x1,x1), "Value" = c(sump[3:1000],meany[3:1000]),"Type"=c(rep("Our estimation",998),rep("Simulation",998)))
d <- ggplot(df1, aes(Length, Value))
d  + geom_point(aes(colour = Type),size=3,alpha=0.8,position=position_jitter(h=0.005, w=0.005)) +  xlab("Length n") + ylab("Estimation of E(L(1,n))")+ theme_grey(base_size = 24)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(size=25))




##B.3 a counterexample where the state is transient
## define a function to generate the desired transition matrix at time t 
simulateTransition<-function(t){
  if (t==1){
    p2=0.5
    p3=0.5
    pt=c(1/2,0.5,p2,p3)
  }
  else{
    pt=c(1/2,0.5,(1/(t^2)),1-(1/(t^2)))
  }
  return (pt)
}

sampleDist = function(n,p) { 
  sample(x = c(1,2), n, replace = T, prob = p) 
}
##define a function to generate all transition matrix from t=1 to lent
generateAll<-function(lent){
  pi=c(1,0)
  alltransi=c()
  sump=0
  for (i in (1:lent)){
    pp=simulateTransition(i)
    pi=pi%*%t(matrix(pp, ncol = 2, nrow = 2))
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
      p1=GA[((i-1)*4+1):((i-1)*4+2)]
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
    else{
      p2=GA[((i-1)*4+3):((i-1)*4+4)]
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
  }
  return (hcount)
}
set.seed(100)
GA=generateAll(10000)
##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:5000)){
  tr=c(tr,LCV4(10000,GA))
}
GA[1]
mean(tr)

GA=generateAll(5000)
##derive the mean of L(1,5000) for 5000 runs
tr=c()
for (i in (1:5000)){
  tr=c(tr,LCV4(5000,GA))
}
GA[1]
mean(tr)


##plot
##B.3 for homogeneous Markov chain with 2 states {1,2}
set.seed(100)
## define a function to get n samples from a given discrete distribution p 
sampleDist = function(n,p) { 
  sample(x = c(1,2), n, replace = T, prob = p) 
}

simulateTransition<-function(t){
  if (t==1){
    p2=0.5
    p3=0.5
    pt=c(1/2,0.5,p2,p3)
  }
  else{
    pt=c(1/2,0.5,(1/(t^2)),1-(1/(t^2)))
  }
  return (pt)
}
## a function to derive the sum of p_11(l) from l=1 to lent

##define a function to generate all transition matrix from t=1 to lent
generateAll1<-function(lent){
  pi=c(1,0)
  alltransi=c()
  sump=c()
  sump1=0
  for (i in (1:lent)){
    pp=simulateTransition(i)
    pi=pi%*%t(matrix(pp, ncol = 2, nrow = 2))
    sump1=sump1+pi[1]
    sump=c(sump,log2(sump1))
  }
  return(sump)
}

##define a function to generate all transition matrix from t=1 to lent
generateAll<-function(lent){
  pi=c(1,0)
  alltransi=c()
  sump=0
  for (i in (1:lent)){
    pp=simulateTransition(i)
    pi=pi%*%t(matrix(pp, ncol = 2, nrow = 2))
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
  s1=c()
  for (i in (1:lent)){
    if (currentpos==1){
      p1=GA[((i-1)*4+1):((i-1)*4+2)]
      currentpos=sampleDist(1,p1)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
      }
    }
    else{
      p2=GA[((i-1)*4+3):((i-1)*4+4)]
      currentpos=sampleDist(1,p2)
      if (currentpos==1){
        ccount=ccount+1
        if (ccount>=hcount){
          hcount=ccount
        }
        s1=c(s1,hcount)
      }
      else{
        ccount=0
        s1=c(s1,hcount)
      }
    }
  }
  return (s1)
}

GA=generateAll(1000)
##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:10000)){
  tr=c(tr,LCV4(1000,GA))
}
s1=t(matrix(tr,ncol=10000))
meany=colMeans(s1,dims=1)
x1=seq(3,1000)
sump=generateAll1(1000)
library(ggplot2)
##plot figure 1 in thesis (comparsion of 3 estimates of E(L(n)) and simulations)
df1 <- data.frame("Length" = c(x1,x1), "Value" = c(sump[3:1000],meany[3:1000]),"Type"=c(rep("Our estimation",998),rep("Simulation",998)))
d <- ggplot(df1, aes(Length, Value))
d  + geom_point(aes(colour = Type),size=3,alpha=0.8,position=position_jitter(h=0.005, w=0.005)) +  xlab("Length n") + ylab("Estimation of E(L(1,n))")+ theme_grey(base_size = 24)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(size=25))







## Illustration for Lemma 3.2
##study the length of longest consecutive visits of a pattern in time nonhomogeneous Markov chain


sampleDist1 = function(n,p) { 
  sample(x = c(1,2,3,4), n, replace = T, prob = p) 
}

simulateTransition1<-function(t){
  p1=1/(4*(t^(0.1)))
  p5=1/(4*(t^(0.1)))
  p6=1/(4*(t^(0.1)))
  p7=0.5-(1/(4*(t^(0.1))))
  p8=0.5-(1/(4*(t^(0.1))))
  pt=c(0.2,0.1,0.35,0.35,0.3,0.2,0.25,0.25,p1,0.5,(1-p1-0.5)/2,(1-p1-0.5)/2,p5,p6,p7,p8)
  return (pt)
}



generateAll1<-function(lent){
  pi=c(1,0,0,0)
  alltransi=c()
  sump=0
  for (i in (1:lent)){
    pp=simulateTransition1(i)
    pi=pi%*%t(matrix(pp, ncol = 4))
    sump=sump+c(pi)[2]
    alltransi=c(alltransi,pp)
  }
  return (c(sump,alltransi))
}



LCV6<- function(lent,GA){
  GA=GA[2:length(GA)]
  currentpos=1
  ccount=0
  hcount=0
  step1=0
  for (i in (1:lent)){
    if (currentpos==1){
      p1=GA[((i-1)*16+1):((i-1)*16+4)]
      currentpos=sampleDist1(1,p1)
      if ((currentpos==3)&(step1=2)){
        step1=3
        ccount=ccount+1
        hcount=max(ccount,hcount)
      }
      else if (currentpos==2){
        step1=1
        ccount=0
      }
      else{
        step1=0
        ccount=0
      }
    }
    else if (currentpos==2){
      p2=GA[((i-1)*16+5):((i-1)*16+8)]
      currentpos=sampleDist1(1,p2)
      if ((step1=1)&(currentpos==1)){
        step1=2
      }
      else if(currentpos==2){
        step1=1
        ccount=0
      }
      else{
        ccount=0
        step1=0
      }
    }
    
    else if (currentpos==3){
      p3=GA[((i-1)*16+9):((i-1)*16+12)]
      currentpos=sampleDist1(1,p3)
      if ((step1=3)&(currentpos==2)){
        step1=1
      }
      else{
        step1=0
        ccount=0
      }
    }
    
    
    else {
      p4=GA[((i-1)*16+13):((i-1)*16+16)]
      currentpos=sampleDist1(1,p4)
      ccount=0
      step1=0
    }
    
  }
  return (hcount)
}

##derive the mean of L(1,1000) for 5000 runs
set.seed(100)
GA=generateAll1(10000)
tr=c()
for (i in (1:5000)){
  tr=c(tr,LCV6(10000,GA))
}
GA[1]
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))



set.seed(100)
GA=generateAll1(100000)
##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:2000)){
  tr=c(tr,LCV6(100000,GA))
}
GA[1]
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))

##plot simulation results from n=1 to 1000
##define a function to generate all transition matrix from t=1 to lent and also the sum of probability moving to set A 
sampleDist1 = function(n,p) { 
  sample(x = c(1,2,3,4), n, replace = T, prob = p) 
}

simulateTransition1<-function(t){
  p1=1/(4*(t^(0.1)))
  p5=1/(4*(t^(0.1)))
  p6=1/(4*(t^(0.1)))
  p7=0.5-(1/(4*(t^(0.1))))
  p8=0.5-(1/(4*(t^(0.1))))
  pt=c(0.2,0.1,0.35,0.35,0.3,0.2,0.25,0.25,p1,0.5,(1-p1-0.5)/2,(1-p1-0.5)/2,p5,p6,p7,p8)
  return (pt)
}


pi=c(1,0,0,0)
alltransi=c()
sum1=0
sump=c()
for (i in (1:1000)){
  pp=simulateTransition1(i)
  pi=pi%*%t(matrix(pp, ncol = 4))
  sum1=sum1+c(pi)[2]
  sump=c(sump,(log(sum1)/log(19.0476)))
  alltransi=c(alltransi,pp)
}



LCV6<- function(lent,GA){
  currentpos=1
  ccount=0
  hcount=0
  step1=0
  ss1=c()
  for (i in (1:lent)){
    if (currentpos==1){
      p1=GA[((i-1)*16+1):((i-1)*16+4)]
      currentpos=sampleDist1(1,p1)
      if ((currentpos==3)&(step1=2)){
        step1=3
        ccount=ccount+1
        hcount=max(ccount,hcount)
        ss1=c(ss1,hcount)
      }
      else if (currentpos==2){
        step1=1
        ccount=0
        ss1=c(ss1,hcount)
      }
      else{
        step1=0
        ccount=0
        ss1=c(ss1,hcount)
      }
    }
    else if (currentpos==2){
      p2=GA[((i-1)*16+5):((i-1)*16+8)]
      currentpos=sampleDist1(1,p2)
      if ((step1=1)&(currentpos==1)){
        step1=2
        ss1=c(ss1,hcount)
      }
      else if(currentpos==2){
        step1=1
        ccount=0
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        step1=0
        ss1=c(ss1,hcount)
      }
    }
    
    else if (currentpos==3){
      p3=GA[((i-1)*16+9):((i-1)*16+12)]
      currentpos=sampleDist1(1,p3)
      if ((step1=3)&(currentpos==2)){
        step1=1
        ss1=c(ss1,hcount)
      }
      else{
        step1=0
        ccount=0
        ss1=c(ss1,hcount) 
      }
    }
    
    
    else {
      p4=GA[((i-1)*16+13):((i-1)*16+16)]
      currentpos=sampleDist1(1,p4)
      ccount=0
      step1=0
      ss1=c(ss1,hcount)
    }
    
  }
  return (ss1)
}

##derive the mean of L(1,1000) for 4000 runs
tr=c()
for (i in (1:10000)){
  tr=c(tr,LCV6(1000,alltransi))
}  
s1=t(matrix(tr,ncol=10000))
meany=colMeans(s1,dims=1)
x1=seq(6,1000)

##plot figure 1 in thesis (comparsion of 3 estimates of E(L(n)) and simulations)
df1 <- data.frame("Length" = c(x1,x1), "Value" = c(sump[6:1000],meany[6:1000]),"Type"=c(rep("Our estimation",995),rep("Simulation",995)))
d <- ggplot(df1, aes(Length, Value))
d  + geom_point(aes(colour = Type),size=3,alpha=0.8,position=position_jitter(h=0.005, w=0.005)) +  xlab("Length n") + ylab("Estimation of E(L(PT,n))")+ theme_grey(base_size = 24)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(size=25))










## Illustration for Lemma 4.2
##study the length of longest consecutive visits of a subset of states in time nonhomogeneous Markov chain

sampleDist1 = function(n,p) { 
  sample(x = c(1,2,3,4), n, replace = T, prob = p) 
}

simulateTransition1<-function(t){
  p1=1/(4*(t^(0.1)))
  p2=1/(4*(t^(0.1)))
  p3=0.5-(1/(4*(t^(0.1))))
  p4=0.5-(1/(4*(t^(0.1))))
  p5=1/(4*(t^(0.1)))
  p6=1/(4*(t^(0.1)))
  p7=0.5-(1/(4*(t^(0.1))))
  p8=0.5-(1/(4*(t^(0.1))))
  pt=c(0.2,0.1,0.35,0.35,0.3,0.2,0.25,0.25,p1,p2,p3,p4,p5,p6,p7,p8)
  return (pt)
}


##define a function to generate all transition matrix from t=1 to lent and also the sum of probability moving to set A 
generateAll1<-function(lent){
  pi=c(1,0,0,0)
  alltransi=c()
  sump=0
  for (i in (1:lent)){
    pp=simulateTransition1(i)
    pi=pi%*%t(matrix(pp, ncol = 4))
    sump=sump+c(pi)[1]+c(pi)[2]
    alltransi=c(alltransi,pp)
  }
  return (c(sump,alltransi))
}

LCV5<- function(lent,GA){
  GA=GA[2:length(GA)]
  currentpos=1
  ccount=0
  hcount=0
  for (i in (1:lent)){
    if (currentpos==1){
      p1=GA[((i-1)*16+1):((i-1)*16+4)]
      currentpos=sampleDist1(1,p1)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else if (currentpos==2){
      p2=GA[((i-1)*16+5):((i-1)*16+8)]
      currentpos=sampleDist1(1,p2)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
    }
    else if (currentpos==3){
      p3=GA[((i-1)*16+9):((i-1)*16+12)]
      currentpos=sampleDist1(1,p3)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
      }
      else{
        ccount=0
      }
      
    }
    else {
      p4=GA[((i-1)*16+13):((i-1)*16+16)]
      currentpos=sampleDist1(1,p4)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
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

##derive the mean of L(1,1000) for 5000 runs
set.seed(100)
GA=generateAll1(10000)
tr=c()
for (i in (1:5000)){
  tr=c(tr,LCV5(10000,GA))
}
GA[1]
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))



set.seed(100)
GA=generateAll1(100000)
##derive the mean of L(1,10000) for 5000 runs
tr=c()
for (i in (1:2000)){
  tr=c(tr,LCV5(100000,GA))
}
GA[1]
mean(tr)
quantile(tr)
quantile(tr, probs = c(0.05, 0.95))


##plot simulation results from n=1 to 1000
##define a function to generate all transition matrix from t=1 to lent and also the sum of probability moving to set A 
sampleDist1 = function(n,p) { 
  sample(x = c(1,2,3,4), n, replace = T, prob = p) 
}

simulateTransition1<-function(t){
  p1=1/(4*(t^(0.1)))
  p2=1/(4*(t^(0.1)))
  p3=0.5-(1/(4*(t^(0.1))))
  p4=0.5-(1/(4*(t^(0.1))))
  p5=1/(4*(t^(0.1)))
  p6=1/(4*(t^(0.1)))
  p7=0.5-(1/(4*(t^(0.1))))
  p8=0.5-(1/(4*(t^(0.1))))
  pt=c(0.2,0.1,0.35,0.35,0.3,0.2,0.25,0.25,p1,p2,p3,p4,p5,p6,p7,p8)
  return (pt)
}


pi=c(1,0,0,0)
alltransi=c()
sum1=0
sump=c()
for (i in (1:1000)){
  pp=simulateTransition1(i)
  pi=pi%*%t(matrix(pp, ncol = 4))
  sum1=sum1+c(pi)[1]+c(pi)[2]
  sump=c(sump,(log(sum1)/log(2.688)))
  alltransi=c(alltransi,pp)
}



LCV6<- function(lent,GA){
  currentpos=1
  ccount=0
  hcount=0
  ss1=c()
  for (i in (1:lent)){
    if (currentpos==1){
      p1=GA[((i-1)*16+1):((i-1)*16+4)]
      currentpos=sampleDist1(1,p1)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        ss1=c(ss1,hcount)
      }
    }
    else if (currentpos==2){
      p2=GA[((i-1)*16+5):((i-1)*16+8)]
      currentpos=sampleDist1(1,p2)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        ss1=c(ss1,hcount)
      }
    }
    else if (currentpos==3){
      p3=GA[((i-1)*16+9):((i-1)*16+12)]
      currentpos=sampleDist1(1,p3)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        ss1=c(ss1,hcount)
      }
      
    }
    else {
      p4=GA[((i-1)*16+13):((i-1)*16+16)]
      currentpos=sampleDist1(1,p4)
      if ((currentpos==1)|(currentpos==2)){
        ccount=ccount+1
        if (ccount>hcount){
          hcount=ccount
        }
        ss1=c(ss1,hcount)
      }
      else{
        ccount=0
        ss1=c(ss1,hcount)
      }
      
    }
  }
  return (ss1)
}

##derive the mean of L(1,5000) for 5000 runs
tr=c()
for (i in (1:10000)){
  tr=c(tr,LCV6(1000,alltransi))
}  
s1=t(matrix(tr,ncol=10000))
meany=colMeans(s1,dims=1)
x1=seq(3,1000)

##plot figure 1 in thesis (comparsion of 3 estimates of E(L(n)) and simulations)
df1 <- data.frame("Length" = c(x1,x1), "Value" = c(sump[3:1000],meany[3:1000]),"Type"=c(rep("Our estimation",998),rep("Simulation",998)))
d <- ggplot(df1, aes(Length, Value))
d  + geom_point(aes(colour = Type),size=3,alpha=0.8,position=position_jitter(h=0.005, w=0.005)) +  xlab("Length n") + ylab("Estimation of E(L(A,n))")+ theme_grey(base_size = 24)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(size=25))







