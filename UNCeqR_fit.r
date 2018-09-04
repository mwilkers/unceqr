library(VGAM)

#description:


args=commandArgs()
print(args)
startInd=6

file=args[startInd]
dnaOnly=args[startInd+1]

pbetabin.ab=pbetabinom.ab

#fit distribution
logLik.bb <- function(ab){
  a = ab[1]
  b = ab[2]
  logL = sum(lchoose(ns,ks) + lbeta(ks+a,ns-ks+b))
  logL = logL - length(ns)*lbeta(a,b)
  -logL
}

pinch=function(x){
  #pinch p-value range
  if(!is.na(x) & x<=0){
    return(1e-16)
  }
  if(!is.na(x) & x >= 1){
    return( 1 - 1e-16)
  }
  return(x)
}

#check duplicates.

d=read.csv(pipe(paste("grep  '^#INIT' ",file, " | cut -d, -f 9,13,23,27,54,55",sep="")),skip=1,header=F)#header=T)
head(d)
for(i in 1:4){
  d[,i] = as.numeric(as.character(d[,i]))
}

#test
#sn=10000
#s=sample(1:100,sn,replace=T)
#v=rbetabinom.ab(sn, s, 0.2, 200, .dontuse.prob = NULL)

#d=cbind(s,v)

#DNA
	ind=which(d[,5]==1)
	ks=d[ind,2]
	ns=apply(d[ind,1:2],1,sum,na.rm=T)
	if(sum(d[ind,2]>ns,na.rm=T)>0){
	  message("dna error, k > n");
	  q();
	}
	summary(ks)
	summary(ns)
	op = optim(par=c(1,10), fn=logLik.bb,control=list(maxit=5000))
	op
	a_d = op$par[1]
	b_d = op$par[2]
	print("DNA")
	print(op)
	print(as.character(a_d))
	print(as.character(b_d))

#RNA
a_r=NA
b_r=NA
if(dnaOnly != "1"){
	ind=which(d[,6]==1)
	ks=d[ind,4]
	ns=apply(d[ind,3:4],1,sum,na.rm=T)
	if(sum(d[ind,4]>ns,na.rm=T)>0){
	  message("rna error, k > n");
	  q();
	}
	summary(ks)
	summary(ns)
	op = optim(par=c(1,10), fn=logLik.bb,control=list(maxit=5000))
	a_r = op$par[1]
	b_r = op$par[2]
	print("RNA")
	print(op)
	print(as.character(a_r))
	print(as.character(b_r))
}

m=rbind(
	c("dna",a_d,b_d),
	c("rna",a_r,b_r))


#write out.
mf=paste(file,".fit",sep="")
write.table(m,sep=",",row.names=F,col.names=F,quote=F,file=mf)



