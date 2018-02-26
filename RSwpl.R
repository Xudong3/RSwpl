
#setting: notation
N1=10 ## number of strata 
N2=10 ##number of elements in each strata (population level)
latitude<-1:N2
longitude<-1:N1
population<-expand.grid(lat=latitude,long=longitude)

population$strata<-population$long
overlap=ceiling(N2*3/5)

model_cluster<-function(population, overlap){
   population$cluster<-numeric(nrow(population))
   
   id<-ifelse(population$lat<=overlap, 
              population$long, 
              ((population$long+population$lat-overlap) %% N1)+1
   )
   population$cluster<-id
   population	
}

population<-model_cluster(population, overlap)
T=length(unique(population$cluster))

##check
table(population$cluster==population$strata)
table(table(population$cluster))
table(table(population$strata))

#Model: parameter from random slope model  model
truebeta1=1
truebeta2=3
truesigma2=4
truetau2_11=1.5
truetau2_12=2
truetau2_22=3
truevalue<-c(truebeta1,truebeta2, truesigma2, truetau2_11, truetau2_12, truetau2_22)
names(truevalue)<-c("beta1", "beta2", "sigma2", "tau2_11", "tau2_12", "tau2_22")

##Population data
#install.packages("MASS")
library("MASS")
#install.packages("rockchalk")
library(rockchalk)
PairCov=matrix(c(truetau2_11, truetau2_12, truetau2_12, truetau2_22), nrow=2, ncol=2, byrow=TRUE)
re=mvrnorm(n=T, mu = c(0,0), Sigma = PairCov)
population$a<-re[,1][population$cluster]
population$b<-re[,2][population$cluster]

population$x<-rnorm(N1*N2)+rnorm(T)[population$cluster]
population$y<-with(population, truebeta1+a+truebeta2*x+b*x+rnorm(N1*N2,s=sqrt(truesigma2)))
population$r=with(population, x*(y-truebeta1-truebeta2*x))

#uninformative stratify sampling design (SRSWOR)
##uninformative means number of elements (sample level) in each strata is fixed
n2=ceiling(N2/10) ##number of elements in each strata (sample level)

#information stratify sampling design (SRSWOR)
##Informative means the number of elments in each stata (sample level) depends on the model variable y,x.
param=c(0.05, 3.5)
n2informative= function(y,x,r, sc, param, N2){
   a=rep(NA, length=length(unique(population$strata)))
   b=rep(NA, length=length(unique(population$strata)))
   for (i in 1: length(unique(sc))){
      a[i]=mean(r[sc==i])
      b[i]=2*ceiling((param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))*N2/2)
   }
   b
}
n2is=n2informative(population$y, population$x, population$r,population$strata, param ,N2)

# Using sampling package for stratify sampling (SRSWOR) 
#install.packages("sampling")
library("sampling")

##uninformative stratify sampling design (SRSWOR) and extracts the observed data
s<-strata(data=population, stratanames="strata", size=rep(n2, N1),
          method="srswor",description=FALSE)
StrSRSWORSample=getdata(population, s)


##informative stratify sampling design (SRSWOR) and extracts the observed data
sis<-strata(data=population, stratanames="strata", size=n2is,
            method="srswor",description=FALSE)
StrSRSWORSampleis=getdata(population, sis)


# Estimation: full-likelihood
#install.packages("lme4")
library(lme4)

##Census full-liklihood estimator (random slope without random intercept)
lmer(y~(1+x|cluster)+x,data=population)



# Estimation: pairwise likelihood (without weight)
l2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22 #pairwise covariance for 12

   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,pc12,0)
   det<-pc11*pc12-iftau^2
    (-(r1*r1*pc-2*r1*r2*iftau+r2*r2*pc22)/2/det-log(det)/2)
}	

dalpha<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22 #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,pc12,0)
   det<-pc11*pc22-iftau^2
   ((r1*pc22-r2*iftau-r1*iftau+r2*pc11)/det)
}

dbeta<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2){
   
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22 #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,pc12,0)
   det<-pc11*pc22-iftau^2
   ( (r1*pc22*x1-r2*iftau*x1-r1*iftau*x2+r2*pc11*x2)/det )
}	

dsigma2<-function(y1,y2, g1,g2, x1,x2,alpha, beta, sigma2,tau2){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22 #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,pc12,0)
   det<-pc11*pc22-iftau^2
   (-( (r1*r1+r2*r2)/2/det - (r1*r1*pc22-2*r1*r2*iftau+r2*r2*pc11)*(pc11+pc22)/2/det/det)- 1/2*(pc11+pc22)/det)
}	

dtau2_11<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22 #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   iftau<-ifelse(g1==g2,tau2,0)
   ifistau<-1*(g1==g2)
   det<-st^2-iftau^2
   ddet<-2*(st-iftau)
   (-( (r1*r1-2*r1*r2*ifistau+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(ddet)/2/det/det) -ddet/det/2)
}	


dtau2_12<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22 #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,tau2,0)
   ifistau<-1*(g1==g2)
   det<-st^2-iftau^2
   ddet<-2*(st-iftau)
   (-( (r1*r1-2*r1*r2*ifistau+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(ddet)/2/det/det) -ddet/det/2)
}	


dtau2_22<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22 #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   ifistau<-1*(g1==g2)
   det<-st^2-iftau^2
   ddet<-2*(st-iftau)
   (-( (r1*r1-2*r1*r2*ifistau+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(ddet)/2/det/det) -ddet/det/2)
}	
