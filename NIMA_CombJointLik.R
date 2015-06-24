################################################################################
#                     Combined Likelihood for family and twin data             #
#                          Brunilda - Sophia Balliu                            #
#                       Leiden University Medical Center                       #
#                                 July 2011                                    #
################################################################################

library(mvtnorm)
library(ltm)

################################################################################
#           Functions To Generate Offsprings'  Phenotypes                      #
################################################################################

  #######################################################################################
  #             Step 1 :    Function to calculate the cases of  NIMA                    #
  #    This function will produce a vector of dim=Numb.Offspring with TRUE/FALSE        #
  #    in case of NIMA effect of not. As NIMA effect we define the case in which        #
  #    the mather but NOT the offspring has the protective allele "d"                   #
  #######################################################################################

NIMA<-function(genotypes) {
    Nof<-dim(genotypes)[2]-2
    sapply(1:Nof,function(i){
    ((genotypes[,1]==1))&(genotypes[,(2+i)]==0)
                            })
                          }
#nima=NIMA(genotypes)

OAE=function(genotypes){
      Nof<-dim(genotypes)[2]-2
      sapply(1:Nof,function(i){
		  t2<-genotypes[,(2+i)]==2
      t1<-genotypes[,(2+i)]==1
      t0<-!(t1|t2)
      1*t2+1*t1+0*t0
       })
       }

#of=OAE(genotypes)



################################################################################
#           Step 2: Penetrance function for one individual                     #
# This function calculates the penetrance  P(Yij|Goij,Gpij) for                #
# individual  j from family i                                                  #
################################################################################

Penetrance_1<-function(u,dbr,oae,nima.par,genotypes) {
  #dbr is the Disease Baseline risk
  #oae is the offspring's allelic effect
  #nima.par is the maternal NIMA effect
    pr=dbr+(oae*OAE(genotypes))+(nima.par*NIMA(genotypes))+u
    return(plogis(pr))
		}
#(Penetrance_1(u,dbr,oae,nima.par,genotypes))
################################################################################
#                                                                              #
#                      Functions to Calculate the Likelihood                   #
#                                                                              #
################################################################################


################################################################################
#                                                                              #
#                       Numerator:  P(Y|Gp,Gc)*Pr(Gc|Gp)P(Gp)                  #
#                                                                              #
################################################################################



###############################################################
#             Step 1 : Penetrance P(Y|Gp,Gc)                  #
###############################################################

Penetrance_fam<-function(u,tu,dbr,oae,nima.par,genotypes,phenotypes) {
        p<-Penetrance_1(u,dbr,oae,nima.par,genotypes)
        B=ifelse(phenotypes,p,1-p)
        lik<-rowSums(log(B))
        f.u<-log(dnorm(u,0,sqrt(tu)))
        return(exp(lik+f.u))
        }

#Penetrance_fam(u,tu,dbr,oae,nima.par,genotypes,phenotypes)

GHQ=function(tu,dbr,oae,nima.par,genotypes,phenotypes) {
npoints=20
N<-dim(genotypes)[1]
GH <- ltm:::gauher(npoints)
wGH<-sqrt(tu)*sqrt(2)*GH$w*exp(GH$x^2)
x=sqrt(tu)*sqrt(2)*GH$x
wGH=matrix(rep(wGH,N),ncol=npoints,byrow=TRUE)
x=matrix(rep(x,N),ncol=npoints,byrow=TRUE)
l=sapply(1:npoints, function(i) Penetrance_fam(x[,i],tu,dbr,oae,nima.par,genotypes,phenotypes)*wGH[,i] )
return(rowSums(l))
}

#GHQ(tu,dbr,oae,nima.par,genotypes,phenotypes)



 ###############################################################
#     Step 2 : Transmission probabilities  Pr(Go|Gp)          #
###############################################################


a<-cbind(rep("2",3),c("2","1","0"))
b<-cbind(rep("1",3),c("2","1","0"))
e<-cbind(rep("0",3),c("2","1","0"))

rows<-rbind(a,b,e)
rows<-paste(rows[,1],rows[,2],sep=" ")
columns<-c("2","1","0")
l<-matrix(rep(NA,27),ncol=3)
colnames(l)<-columns
rownames(l)<-rows

l[1,]<-c(1,0,0)
l[2,]<-c(1/2,1/2,0)
l[3,]<-c(0,1,0)
l[4,]<-c(1/2,1/2,0)
l[5,]<-c(1/4,1/2,1/4)
l[6,]<-c(0,1/2,1/2)
l[7,]<-c(0,1,0)
l[8,]<-c(0,1/2,1/2)
l[9,]<-c(0,0,1)

tra.prob<-function(genotypes) {
      Nof<-dim(genotypes)[2]-2
      N<-dim(genotypes)[1]
      par.gen<-paste(genotypes[,1],genotypes[,2],sep=" ")
      off.gen<-genotypes[,3:(2+Nof)]
      p<-sapply(1:Nof, function(j) {
        a<-(matrix(rep(rownames(l),N),nrow=N,byrow=TRUE))==par.gen
        b<-(matrix(rep(colnames(l),N),nrow=N,byrow=TRUE))==as.character(off.gen[,j])
        proba=sapply(1:N, function(i){
            rows=a[i,]
            cols=b[i,]
            l[rows,cols]
                                      })
        return(proba)
                                   })
  return(rowSums(log(p)))
}

#tra.prob(genotypes)

tra.prob.twins<-function(genotypes) {
      Nof<-dim(genotypes)[2]-2
      N<-dim(genotypes)[1]
      par.gen<-paste(genotypes[,1],genotypes[,2],sep=" ")
      off.gen<-genotypes[,3:(2+Nof)]
      p<-sapply(1:Nof, function(j) {
        a<-(matrix(rep(rownames(l),N),nrow=N,byrow=TRUE))==par.gen
        b<-(matrix(rep(colnames(l),N),nrow=N,byrow=TRUE))==as.character(off.gen[,j])
        proba=sapply(1:N, function(i){
            rows=a[i,]
            cols=b[i,]
            l[rows,cols]
                                      })
        return(proba)
                                   })
  return(apply(p,1,prod))
}

#tra.prob.twins(genotypes)






###############################################################
#        Step 3: Calculate the P(Gp) probability           #
###############################################################

P.Gp<-function(a.f,genotypes){
ind.gen<-genotypes[,1:2]
I1<-ind.gen==2
I2<-ind.gen==1
I3<-ind.gen==0
q1<-a.f^2
q2<-2*a.f*(1-a.f)
q3<-(1-a.f)^2
lik<-(q1^I1)*(q2^I2)*(q3^I3)
return(rowSums(log(lik)))
}

#P.Gp(a.f,genotypes)

###############################################################
#       step 4:  Numerator=log(P(Y|Go,Gp)*P(Go|Gp)*P(Gp))     #
###############################################################

Numerator<-function(a.f,tu,dbr,oae,nima.par,genotypes,phenotypes) {
  penetrance<-log(GHQ(tu,dbr,oae,nima.par,genotypes,phenotypes))
  p.tran<-tra.prob(genotypes)
  p.gp<-P.Gp(a.f,genotypes)
  return(penetrance+p.tran+p.gp)
  }

#Numerator(a.f,tu,dbr,oae,nima.par,genotypes,phenotypes)

Numerator.twins<-function(a.f, tu, dbr, oae, nima.par,p.gs,genotypes,phenotypes) {
  N<-dim(genotypes)[1]
  num=sapply(1:9, function(k){
  p.gss=p.gs[k,1:2]
  weight=p.gs[k,3]
  genotype=cbind(matrix(rep(p.gss,N),ncol=2,byrow=TRUE),genotypes)
  penetrance<-GHQ(tu,dbr,oae,nima.par,genotype,phenotypes)
  p.tran<-tra.prob.twins(genotype)
  exp(log(penetrance)+log(p.tran)+2*log(weight))
  })
  return(log(rowSums(num))+log(sum(p.gs[,3]^2)))
  }

#Numerator.twins(a.f, tu, dbr, oae, nima.par,p.gs,genotypes,phenotypes)


########################################################################
#                                                                      #
#                          Denominator=                                #
#  1-logSum_{Gc,Gp}(P(Go|Gp)*P(Gp)*{P(SumY=1|Gc,Gp)+ P(SumY=0|Gc,Gp)}) #
#                                                                      #
########################################################################

###############################################################
#        Step 1: Calculate all possible combinations          #
###############################################################

fam.comb.f<-function(Nof) {
n=Nof+2
to<-vector("list",n)
ti<-lapply(to,function(x) c(0,1,2))
gts<-expand.grid(ti)
return(gts)
}


#fam.comb=fam.comb.f(Nof)
###############################################################
#        Step 2: Exclude impossible combinations              #
#   Imposible according to Mendelian Inheritance              #
###############################################################

exclude.f<-function(fam.comb) {
  p.tran=exp(tra.prob(fam.comb))
  t<-p.tran!=0
  p.tran=p.tran[t]
  fam.comb<-list(genotypes=fam.comb[t,],p.tran=p.tran)
  return(fam.comb)
}

#exclude=exclude.f(fam.comb)




################################################################################
#            Step 3: Correct for Ascertainment                                 #
################################################################################

#  Probability that all offspring will be unaffected
all.unaffected<-function(tu,dbr,oae,nima.par,genotypes) {
Nof<-dim(genotypes)[2]-2
N=dim(genotypes)[1]
phenotypes<-matrix(rep(0,Nof*N),N,byrow=TRUE)
return(GHQ(tu,dbr,oae,nima.par,genotypes,phenotypes))
}
#all.unaffected(tu,dbr,oae,nima.par,genotypes)



#  Probability that exactly one offspring will be affected
one.affected<-function(tu,dbr,oae,nima.par,genotypes) {
Nof<-dim(genotypes)[2]-2
N=dim(genotypes)[1]
pheno<-diag(1,Nof,Nof)
p<-sapply(1:Nof,function(i) {
    phenotypes<-matrix(rep(pheno[,i],N),N,byrow=TRUE)
    GHQ(tu,dbr,oae,nima.par,genotypes,phenotypes)
    })
return(rowSums(p))
}
#one.affected(tu,dbr,oae,nima.par,genotypes)




###############################################################
#       step 4:  Denominator                                  #
###############################################################

Denominator<-function(a.f,tu,dbr,oae,nima.par,Nof) {
exclude=exclude.f(fam.comb.f(Nof))
genotypes=exclude$genotypes
p.0<-all.unaffected(tu,dbr,oae,nima.par,genotypes)
p.1<-one.affected(tu,dbr,oae,nima.par,genotypes)
p.tran<-exp(tra.prob(genotypes))
p.gp<-exp(P.Gp(a.f,genotypes))
denom=log(1-sum(p.tran*p.gp*(p.1+p.0)))
return(denom)
}

#Denominator(a.f,tu,dbr,oae,nima.par,Nof)



Denominator.twins<-function(a.f,tu,dbr,oae,nima.par,Nof) {
exclude=exclude.f(fam.comb.f(Nof))
genotypes=exclude$genotypes
p.0<-all.unaffected(tu,dbr,oae,nima.par,genotypes)
p.tran<-exclude$p.tran
p.gp<-exp(P.Gp(a.f,genotypes))
denom=log(1-sum(p.tran*p.gp*p.0))
return(denom)
}

#Denominator.twins(a.f,tu,dbr,oae,nima.par,Nof)

##################################################################################
#####################    Likelihood for all families           ###################
##################################################################################


Log.likelihood<-function(a.f,tu,dbr,oae,nima.par,genotypes,phenotypes){
Nof<-dim(genotypes)[2]-2
N<-dim(genotypes)[1]
-sum(Numerator(a.f,tu,dbr,oae,nima.par,genotypes,phenotypes))+N*Denominator(a.f,tu,dbr,oae,nima.par,Nof)
}

Log.likelihood.twins<-function(a.f,tu ,dbr ,oae ,nima.par ,p.gs,genotypes.t,phenotypes.t){
N<-dim(genotypes.t)[1]
-sum(Numerator.twins(a.f,tu,dbr,oae,nima.par,p.gs,genotypes.t,phenotypes.t))+N*Denominator.twins(a.f,tu,dbr,oae,nima.par,Nof=2)
}


P.Gp.num<-function(a.f,genotypes){
I1<-genotypes==2
I2<-genotypes==1
I3<-genotypes==0
q1<-a.f^2
q2<-2*a.f*(1-a.f)
q3<-(1-a.f)^2
lik<-(q1^I1)*(q2^I2)*(q3^I3)
return(apply(lik,1,prod))
}

Log.likelihood2<-function(p,genotypes,phenotypes,genotypes.t,phenotypes.t){
a.f=inv.logit(p[1])
tu=exp(p[2])
dbr=p[3]
oae=p[4]
nima.par=p[5]

a<-cbind(rep(2,3),c(2,1,0))
b<-cbind(rep(1,3),c(2,1,0))
e<-cbind(rep(0,3),c(2,1,0))
p.gs<-rbind(a,b,e)
weights=P.Gp.num(a.f,p.gs)
p.gs=cbind(p.gs,weights)
Log.likelihood(a.f=a.f,tu=tu,dbr=dbr,oae=oae,nima.par=nima.par,genotypes,phenotypes)+Log.likelihood.twins(a.f=a.f,tu=tu,dbr=dbr,oae=oae,nima.par=nima.par,p.gs,genotypes.t,phenotypes.t)
}




