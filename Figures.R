#######File that generates Manuscript Figures###################################
#######Mark Albrecht, Ansu Chatterjee        ###################################

################################################################################
###Set directory path and load libraries
################################################################################
setwd("K:/xover/Haplotype/Imputation_SHF")
rm(list=ls(all=TRUE)) 
gc()
graphics.off()
dev.off()
library(lattice)
library(car)
library(alr3)
library(MASS)
library(faraway)
library(quantreg)
library(latticeExtra)

#######################################
#######################################
###Support Example:  normal hyperparameter, multinomial sample, non-cartesian support
#######################################
#######################################
mean=0
sd=3
lim<-8
from=-10
to=10
by=0.5
grid<-data.frame("X"=seq(from=from,to=to,by=by),"Y"=seq(from=from,to=to,by=by))
mesh<-expand.grid(grid)
mesh<-mesh[which(apply(mesh,1,sum)<=lim & apply(mesh,1,sum)>=-lim),]
plot(mesh)
dmesh<-t(apply(mesh,1,dnorm,mean=mean,sd=sd)) ###calculate marginals
dmesh<-apply(dmesh,1,prod) ###calculate joint prob
dmesh<-(1/sum(dmesh))*dmesh ###normalize to sum to 1
sample<-rmultinom(n=1,size=100000,prob=dmesh)
pmf_sample<-sample/100000
true_marginal<-matrix(nrow=nrow(grid),ncol=ncol(grid))
sample_marginal<-matrix(nrow=nrow(grid),ncol=ncol(grid))
rownames(true_marginal)<-grid$X
colnames(true_marginal)<-colnames(grid)
rownames(sample_marginal)<-grid$X
colnames(sample_marginal)<-colnames(grid)
for (i in 1:nrow(grid)){
true_marginal[i,1]<-sum(dmesh[which(mesh$X==grid$X[i])])
true_marginal[i,2]<-sum(dmesh[which(mesh$Y==grid$Y[i])])
sample_marginal[i,1]<-sum(pmf_sample[which(mesh$X==grid$X[i])])
sample_marginal[i,2]<-sum(pmf_sample[which(mesh$Y==grid$Y[i])])
}
true_marginal
sample_marginal

Z_model_true<-matrix(nrow=nrow(mesh),ncol=5)
Z_model_sample<-matrix(nrow=nrow(mesh),ncol=5)
colnames(Z_model_true)<-c("X","Y","[XY]","[X]","[Y]")
colnames(Z_model_sample)<-c("X","Y","[XY]","[X]","[Y]")
for (i in 1:nrow(mesh)){
Z_model_true[i,]<-cbind(mesh[i,1],mesh[i,2],dmesh[i],
true_marginal[which(mesh[i,1]==grid$X),1],
true_marginal[which(mesh[i,2]==grid$Y),1])
Z_model_sample[i,]<-cbind(mesh[i,1],mesh[i,2],pmf_sample[i],
sample_marginal[which(mesh[i,1]==grid$X),1],
sample_marginal[which(mesh[i,2]==grid$Y),1])
}
Z_model_true<-as.data.frame(Z_model_true)
Z_model_sample<-as.data.frame(Z_model_sample)

######### Fit Copula Models ########
m_true<-lm(`[XY]`~`[X]`*`[Y]`,data=Z_model_true)
summary(m_true)
m_sample<-lm(`[XY]`~`[X]`*`[Y]`,data=Z_model_sample)
summary(m_sample)

#####Proper Estimates#######
trellis.device(color=F)
contourplot(`[XY]`~X*Y,data=Z_model_true,main="True Density",
region=TRUE,pretty=T,cuts=10,labels=TRUE)
contourplot(`[XY]`~X*Y,data=Z_model_sample,main="Empirical Density (n=100,000)",
region=TRUE,pretty=T,cuts=10,labels=TRUE)


#####Regression Model Prediction Over Propper Support#####
dat_true<-data.frame("Y_hat"=predict(m_true),Z_model_true)
dat_sample<-data.frame("Y_hat"=predict(m_sample),Z_model_true)
trellis.device(color=F)
contourplot(Y_hat~X*Y,data=dat_true,main="Model Fit on True Density",
region=TRUE,pretty=T,cuts=10,labels=TRUE)
contourplot(Y_hat~X*Y,data=dat_sample,main="Model Fit on Empirical Density",
region=TRUE,pretty=T,cuts=10,labels=TRUE)

#####Impropper Extrapolation of Regression Model######
impropper_support<-expand.grid(grid)
Z_model_t<-matrix(nrow=nrow(impropper_support),ncol=4)
Z_model_s<-matrix(nrow=nrow(impropper_support),ncol=4)
colnames(Z_model_t)<-c("X","Y","[X]","[Y]")
colnames(Z_model_s)<-c("X","Y","[X]","[Y]")
for (i in 1:nrow(impropper_support)){
Z_model_t[i,]<-cbind(impropper_support[i,1],impropper_support[i,2],
true_marginal[which(impropper_support[i,1]==grid$X),1],
true_marginal[which(impropper_support[i,2]==grid$Y),1])
Z_model_s[i,]<-cbind(impropper_support[i,1],impropper_support[i,2],
sample_marginal[which(impropper_support[i,1]==grid$X),1],
sample_marginal[which(impropper_support[i,2]==grid$Y),1])
}
b_t<-data.frame("[X]"=Z_model_t[,3],"[Y]"=Z_model_t[,4])
colnames(b_t)<-c("[X]","[Y]")
dat_imp_t<-data.frame("Y_hat"=predict(m_true,newdata=b_t),Z_model_t)
b_s<-data.frame("[X]"=Z_model_s[,3],"[Y]"=Z_model_s[,4])
colnames(b_s)<-c("[X]","[Y]")
dat_imp_s<-data.frame("Y_hat"=predict(m_sample,newdata=b_s),Z_model_s)

trellis.device(color=F)
contourplot(Y_hat~X*Y,data=dat_imp_t,main="Impropper Extrapolation on True Density",
region=TRUE,pretty=T,cuts=10,labels=TRUE)
contourplot(Y_hat~X*Y,data=dat_imp_s,main="Impropper Extrapolation on Empirical Density",
region=TRUE,pretty=T,cuts=10,labels=TRUE)




#######################################
#######################################
###Figure 1
#######################################
#######################################

################################################################################
###Set directory path and load libraries
################################################################################
setwd("K:/xover/Haplotype/Imputation_SHF")
rm(list=ls(all=TRUE)) 
gc()
graphics.off()
dev.off()
library(lattice)
library(car)
library(alr3)
library(MASS)
library(faraway)
library(quantreg)
library(latticeExtra)


###################################################################################
###Load coefficient data
################################################################################
graphics.off()
filename=paste0(getwd(),"/output/2013-01-18Coefficients_OLS.csv")
COEFF<-read.csv(file=filename, check.names=FALSE)
COEFF_names<-gsub("`","",colnames(COEFF))
COEFF_names<-gsub("~",",",COEFF_names)
COEFF_names<-gsub(":","][",COEFF_names)
COEFF_names<-c(COEFF_names[c(1,2)],paste("log([",COEFF_names[-c(1,2)],"])",sep=""))
colnames(COEFF)<-COEFF_names
head(COEFF)
race<-as.character(unique(COEFF$Race))
trellis.device(color=T)
###################################################################################
###Generate Coefficient Graph for CAU
################################################################################
coeff<-COEFF[which(COEFF$Race=="NAM"),]
rownames(coeff)<-coeff$Level
coeff<-t(coeff[,-c(1,2)])
labels<-rep(rownames(coeff),4)
plot_data<-cbind(labels,stack(as.data.frame(coeff)))
plot_data<-plot_data[which(!is.na(plot_data$values)),]
plot_data$level="tmp"
plot_data$level[which(plot_data$ind=="4!")]="r=4!"
plot_data$level[which(plot_data$ind=="3!")]="r=3!"
plot_data$level[which(plot_data$ind=="2!")]="r=2!"
plot_data$level[which(plot_data$ind=="1!")]="r=1!"
head(plot_data)
xx<-data.frame()
xx=rbind(xx,plot_data[order(abs(plot_data$value[which(plot_data$ind=="4!")])),])
xx=rbind(xx,plot_data[order(abs(plot_data$value[which(plot_data$ind=="3!")]))+nrow(xx),])
xx=rbind(xx,plot_data[order(abs(plot_data$value[which(plot_data$ind=="2!")]))+nrow(xx),])
xx=rbind(xx,plot_data[order(abs(plot_data$value[which(plot_data$ind=="1!")]))+nrow(xx),])
xx$index=1:nrow(xx)
names<-list(p1=xx$labels[which(xx$ind=="1!")],p2=xx$labels[which(xx$ind=="2!")],
p3=xx$labels[which(xx$ind=="3!")],p4=xx$labels[which(xx$ind=="4!")])


barchart(index~values|level,data=xx,horiz=T,layout=c(1,4),origin=0,ylab="Model Components",
xlab=expression("Coefficient Values:"~~beta*"'s"),
scales=list(y=list(relation="free",labels=names,cex=0.6)))
resizePanels()


#####Make individual bar charts#####

barchart(index~values|level,data=xx[which(xx$ind=="4!"),],horiz=T,layout=c(1,1),origin=0,ylab="Model Components",
         xlab=expression("Coefficient Values:"~~beta*"'s"),strip=FALSE,
         scales=list(y=list(relation="free",labels=names$p4,cex=1)))
resizePanels()

barchart(index~values|level,data=xx[which(xx$ind=="3!"),],horiz=T,layout=c(1,1),origin=0,ylab="Model Components",
         xlab=expression("Coefficient Values:"~~beta*"'s"),strip=FALSE,
         scales=list(y=list(relation="free",labels=names$p3,cex=1)))
resizePanels()

barchart(index~values|level,data=xx[which(xx$ind=="2!"),],horiz=T,layout=c(1,1),origin=0,ylab="Model Components",
         xlab=expression("Coefficient Values:"~~beta*"'s"),strip=FALSE,
         scales=list(y=list(relation="free",labels=names$p2,cex=1)))
resizePanels()

barchart(index~values|level,data=xx[which(xx$ind=="1!"),],horiz=T,layout=c(1,1),origin=0,ylab="Model Components",
         xlab=expression("Coefficient Values:"~~beta*"'s"),strip=FALSE,
         scales=list(y=list(relation="free",labels=names$p1,cex=1)))
resizePanels()



#######################################
#######################################
###Figure 2 CAU imputed coeff
#######################################
#######################################

################################################################################
###Set directory path and load libraries
################################################################################
setwd("K:/xover/Haplotype/Imputation_SHF")
rm(list=ls(all=TRUE)) 
gc()
graphics.off()
dev.off()
library(lattice)
library(car)
library(alr3)
library(MASS)
library(faraway)
library(quantreg)

###################################################################################
###Load coefficient data to subset data matrix (will not be needed in code module)
################################################################################
filename=paste0(getwd(),"/output/2013-01-18Coefficients_OLS.csv")
filename_train=paste0(getwd(),"/output/2013-01-18Index_Train_OLS.csv")
filename_valid=paste0(getwd(),"/output/2013-01-18Index_Validation_OLS.csv")
COEFF<-read.csv(file=filename, check.names=FALSE)
COEFF_names<-gsub("`","",colnames(COEFF))
colnames(COEFF)<-COEFF_names
head(COEFF)
race<-as.character(unique(COEFF$Race))
r_ind<-12
filename=paste0(getwd(),"/nemo_est5loc_parts/",race[r_ind],".ARS.h5_component.freqs.csv")
DATA<-read.csv(file=filename, check.names=FALSE, nrows=50000,row.names=1)
index_train=read.csv(file=filename_train, check.names=FALSE ,row.names=1)
index_validation=read.csv(file=filename_valid, check.names=FALSE ,row.names=1)
pull<-index_validation[which(rownames(index_validation)==race[r_ind]),]
pull_t<-index_train[which(rownames(index_train)==race[r_ind]),]
data<-DATA[as.numeric(pull),]
data_t<-DATA[as.numeric(pull_t),]
density=c("5!",rep("4!",5),rep("3!",10),rep("2!",10),rep("1!",5))
names(density)<-colnames(data)
MSEp<-numeric()
b=1 #set to 4! when b=1

#######################################################################
###Strip out Data to test Imputation Approach
#######################################################################
strip=c(0.2,0.2,0.2,0.2) ###precent of data to strip from (4!,3!,2!,1!)
data_s<-as.matrix(data)
indexing<-unique(density)[-1]
i<-1
for (i in 1:length(strip)) {
act<-which(density==indexing[i])
pull<-sample(1:length(data_s[,act]),length(data_s[,act])*strip[i])
tmp<-data_s[,act]
tmp[pull]=NA
data_s[,act]=tmp  
}
head(data_s)


#######################################################################
###Start prediction loop over levels 4!, 3!, 2!, 1!
#######################################################################

levels<-unique(density)[-1]
FREQUENCY_hat<-matrix(nrow=nrow(data_s),ncol=length(levels))
PERCENT_miss<-matrix(nrow=nrow(data_s),ncol=length(levels))
colnames(FREQUENCY_hat)<-levels
rownames(FREQUENCY_hat)<-rownames(data_s)
colnames(PERCENT_miss)<-levels
rownames(PERCENT_miss)<-rownames(data_s)
coeff<-COEFF[which(COEFF$Race==race[r_ind] & COEFF$Level==levels[b]),]
coeff=coeff[which(is.na(coeff)==FALSE)]
coeff_names<-colnames(coeff)
int<-0
if (any(coeff_names=="(Intercept)")){
int<-1
coeff=coeff[-which(coeff_names=="(Intercept)")]
coeff_names=coeff_names[-which(coeff_names=="(Intercept)")]
}


#######################################################################
###This is where the code starts for imputation of missing values using MVN density assumption
#######################################################################
X<-log(data[,-1])
XX<-log(data_s[,-1])
XT<-log(data_t[,-1])
XXX<-XX
V<-rbind(XX,XT)
mu<-colMeans(V,na.rm=T)
sigma<-var(V,na.rm=T)
print(mu)
print(sigma)
p<-ncol(XX)
for (i in 1:nrow(XX)){
index<-which(is.na(XX[i,]))
q<-length(index)
if (q>=1 & p!=q){
T<-c(XX[i,-index],XX[i,index])
Mu<-c(mu[-index],mu[index])
Si<-rbind(sigma[-index,],sigma[index,])
Si<-cbind(Si[,-index],Si[,index])
m_miss<-matrix(Mu[(p+1-q):p],nrow=q,ncol=1)
var_prod<-Si[(p+1-q):p,1:(p-q)]%*%solve(Si[1:(p-q),1:(p-q)])
shift<-matrix(T[1:(p-q)]-Mu[1:(p-q)],nrow=p-q,ncol=1)
new<-m_miss+var_prod%*%shift
XX[i,index]=new}else{
if (p==q){
XX[i,index]=mu
}}
}
head(XX)


#######################################################################
###Plot of Imputed Values Versus Observed
#######################################################################
labels<-colnames(XX)
labels<-gsub("~",",",labels)
labels<-gsub(":","][",labels)
labels<-paste("[",labels,"]",sep="")
par(mfrow=c(6,5),mar=c(3,3.5,3,1)+0.1,mgp=c(2,1,0))
for (i in 1:ncol(XX)){
Y<-X[is.na(XXX[,i]),i]
Y_hat<-XX[is.na(XXX[,i]),i]
if (length(Y) > 0) {
plot(Y_hat~Y,main=list(paste0(labels[i]),cex=1),xlab=expression("log("*X*")"),ylab=expression("log("*hat(X)*")"),cex=0.15)
maxx<-max(Y)
minx<-min(Y)
maxy<-max(Y_hat)
miny<-min(Y_hat)
R<-(cor(Y,Y_hat))^2
R<-floor(R*100)
R<-R/100
tmpcmd<-paste("R^2 ==",R)
tmpcmd<-paste("text(minx+(maxx-minx)*0.85,miny+(maxy-miny)*0.2,expression(",tmpcmd,"),cex=1.25)",sep="")
eval(parse(text = tmpcmd)) 
m1<-lm(Y_hat~Y)
abline(m1)
}
}



#######################################
#######################################
###Figure 3
#######################################
#######################################
################################################################################
###Set directory path and load libraries
################################################################################
setwd("K:/xover/Haplotype/Imputation_SHF")
rm(list=ls(all=TRUE)) 
gc()
graphics.off()
dev.off()
library(lattice)
library(car)
library(alr3)
library(MASS)
library(faraway)
library(quantreg)

###################################################################################
###Load coefficient data to subset data matrix (will not be needed in code module)
################################################################################
filename=paste0(getwd(),"/output/2013-01-18Coefficients_OLS.csv")
filename_train=paste0(getwd(),"/output/2013-01-18Index_Train_OLS.csv")
filename_valid=paste0(getwd(),"/output/2013-01-18Index_Validation_OLS.csv")
COEFF<-read.csv(file=filename, check.names=FALSE)
COEFF_names<-gsub("`","",colnames(COEFF))
colnames(COEFF)<-COEFF_names
head(COEFF)
par(mfrow=c(6,5),mar=c(3,3.5,3,1)+0.1,mgp=c(2,1,0))
#par(mfrow=c(5,6))
race<-as.character(unique(COEFF$Race))
race<-race[-which(race=="UNK")]
race<-race[-which(race=="DEC")]
race<-race[-which(race=="OTH")]



r_ind<-2
for (r_ind in 1:length(race)){###Start Race Loop 
filename=paste0(getwd(),"/nemo_est5loc_parts/",race[r_ind],".ARS.h5_component.freqs.csv")
DATA<-read.csv(file=filename, check.names=FALSE, nrows=50000,row.names=1)
###index_train=read.csv(file=filename_train, check.names=FALSE ,row.names=1)
index_validation=read.csv(file=filename_valid, check.names=FALSE ,row.names=1)
pull<-index_validation[which(rownames(index_validation)==race[r_ind]),]
data<-DATA[as.numeric(pull),]
data<-DATA
head(data)
density=c("5!",rep("4!",5),rep("3!",10),rep("2!",10),rep("1!",5))
names(density)<-colnames(data)


#######################################################################
###Start prediction loop over levels 4!
#######################################################################

levels<-unique(density)[-1]
b<-4 #set level to 4! equals 1 (reverse order)
coeff<-COEFF[which(COEFF$Race==race[r_ind] & COEFF$Level==levels[b]),]
coeff=coeff[which(is.na(coeff)==FALSE)]
coeff_names<-colnames(coeff)
int<-0
if (any(coeff_names=="(Intercept)")){
int<-1
coeff=coeff[-which(coeff_names=="(Intercept)")]
coeff_names=coeff_names[-which(coeff_names=="(Intercept)")]
}

#########################################################################
###Generate Data Matrix from densities for unstripped data
#########################################################################

X<-matrix(nrow=nrow(data),ncol=(length(coeff)-2))
for (i in 3:length(coeff)) {
q=strsplit(coeff_names[i],":")
tmp<-do.call(rbind, q)
Z=matrix(nrow=nrow(data),ncol=length(tmp))
for (ii in 1:length(tmp)){
Z[,ii]=data[,which(colnames(data)==tmp[ii])]
}
Z=apply(Z,1,prod)
X[,i-2]=Z
}
X<-as.matrix(X)
colnames(X)<-coeff_names[-c(1,2)]
X<-log(X)
head(X)

#######################################################################
###Use Regression Model Results to get Y versus Y_hat for data
#######################################################################
if (int==1){
X=cbind(1,X)
coeff<-COEFF[which(COEFF$Race==race[r_ind] & COEFF$Level==levels[b]),]
coeff=coeff[which(is.na(coeff)==FALSE)]
coeff_names<-colnames(coeff)
}
Freq_hat<-X%*%t(coeff[3:length(coeff)])
Freq<-log(data[,1])
MSEp<-sum((Freq-Freq_hat)^2)/nrow(data)
R<-(cor(Freq,Freq_hat))^2
plot(Freq_hat~Freq,main=list(paste0(race[r_ind]), cex=1),xlab=expression("log("*Y*")"),ylab=expression("log("*hat(Y)*")"),cex=0.1)
maxx<-max(Freq)
minx<-min(Freq)
maxy<-max(Freq_hat)
miny<-min(Freq_hat)
R<-floor(R*100)
R<-R/100
tmpcmd<-paste("R^2 ==",R)
tmpcmd<-paste("text(minx+(maxx-minx)*0.85,miny+(maxy-miny)*0.2,expression(",tmpcmd,"),cex=1.25)",sep="")
eval(parse(text = tmpcmd)) 
###plot(exp(Freq_hat)~exp(Freq),main=list(paste0(race[r_ind]), cex=1),xlab=expression("Y"),ylab=expression(""*hat(Y)*""))
m1<-lm(Freq_hat~Freq)
abline(m1)
}#####end race loop












par(mfrow=c(1,2),mar=c(3,3.25,3,1)+0.1,mgp=c(2,1,0))
plot(Freq_hat~Freq,main=list(paste0("Log(Y):",race[r_ind]), cex=1),xlab=list(expression(Y),cex=1),ylab=list(expression(hat(Y)),cex=1))
plot(exp(Freq_hat)~exp(Freq),main=list(paste0("Y:",race[r_ind]), cex=1),xlab=list(expression(Y),cex=1),ylab=list(expression(hat(Y)),cex=1))
MSEp<-sum((exp(Freq)-exp(Freq_hat))^2)/nrow(data)
print(MSEp)
tmp_out<-cbind(Freq,Freq_hat,exp(Freq),exp(Freq_hat))
rownames(tmp_out)<-rownames(data)
colnames(tmp_out)<-c("log(Y)","log(Y_hat)","Y","Y_hat")
head(tmp_out)
write.csv(tmp_out, file=paste0(getwd(),"/validaiton_joel.csv"),row.names=TRUE,col.names=TRUE)



#######################################
#######################################
###Figure 4
#######################################
#######################################

################################################################################
###Set directory path and load libraries
################################################################################
setwd("K:/xover/Haplotype/Imputation_SHF")
rm(list=ls(all=TRUE)) 
gc()
graphics.off()
dev.off()
library(lattice)
library(car)
library(alr3)
library(MASS)
library(faraway)
library(quantreg)

###################################################################################
###Load coefficient data to subset data matrix (will not be needed in code module)
################################################################################
filename=paste0(getwd(),"/output/2013-01-18Coefficients_OLS.csv")
filename_train=paste0(getwd(),"/output/2013-01-18Index_Train_OLS.csv")
filename_valid=paste0(getwd(),"/output/2013-01-18Index_Validation_OLS.csv")
COEFF<-read.csv(file=filename, check.names=FALSE)
COEFF_names<-gsub("`","",colnames(COEFF))
colnames(COEFF)<-COEFF_names
head(COEFF)
race<-as.character(unique(COEFF$Race))
r_ind<-12
filename=paste0(getwd(),"/nemo_est5loc_parts/",race[r_ind],".ARS.h5_component.freqs.csv")
DATA<-read.csv(file=filename, check.names=FALSE, nrows=50000,row.names=1)
index_train=read.csv(file=filename_train, check.names=FALSE ,row.names=1)
index_validation=read.csv(file=filename_valid, check.names=FALSE ,row.names=1)
pull<-index_validation[which(rownames(index_validation)==race[r_ind]),]
pull_t<-index_train[which(rownames(index_train)==race[r_ind]),]
data<-DATA[as.numeric(pull),]
data_t<-DATA[as.numeric(pull_t),]
density=c("5!",rep("4!",5),rep("3!",10),rep("2!",10),rep("1!",5))
names(density)<-colnames(data)
strip_all<-matrix(nrow=0,ncol=5)
strip_step<-seq(0,1,length=50)
for (i in 1:4){
tmp<-matrix(1,nrow=length(strip_step),ncol=4)
tmp[,i:4]=strip_step
tmp<-cbind(tmp,i)
strip_all<-rbind(strip_all,tmp)
}
MSEp<-numeric()
b=1
c=2
for (b in 1:4){### start levels loop
st<-strip_all[which(strip_all[,5]==b),]
for (c in 1:nrow(st)){ ###start loop for plotting MSEp by percent missing

#######################################################################
###Strip out Data to test Imputation Approach
#######################################################################
strip=st[c,1:4] ###precent of data to strip from (4!,3!,2!,1!)
data_s<-as.matrix(data)
indexing<-unique(density)[-1]
i<-1
for (i in 1:length(strip)) {
act<-which(density==indexing[i])
pull<-sample(1:length(data_s[,act]),length(data_s[,act])*strip[i])
tmp<-data_s[,act]
tmp[pull]=NA
data_s[,act]=tmp  
}
head(data_s)


#######################################################################
###Start prediction loop over levels 4!, 3!, 2!, 1!
#######################################################################

levels<-unique(density)[-1]
FREQUENCY_hat<-matrix(nrow=nrow(data_s),ncol=length(levels))
PERCENT_miss<-matrix(nrow=nrow(data_s),ncol=length(levels))
colnames(FREQUENCY_hat)<-levels
rownames(FREQUENCY_hat)<-rownames(data_s)
colnames(PERCENT_miss)<-levels
rownames(PERCENT_miss)<-rownames(data_s)
coeff<-COEFF[which(COEFF$Race==race[r_ind] & COEFF$Level==levels[b]),]
coeff=coeff[which(is.na(coeff)==FALSE)]
coeff_names<-colnames(coeff)
int<-0
if (any(coeff_names=="(Intercept)")){
int<-1
coeff=coeff[-which(coeff_names=="(Intercept)")]
coeff_names=coeff_names[-which(coeff_names=="(Intercept)")]
}


#######################################################################
###This is where the code starts for imputation of missing values using MVN density assumption
#######################################################################
ZZ<-log(data_s[,-1])
ZT<-log(data_t[,-1])
ZZZ<-ZZ
V<-rbind(ZZ,ZT)
mu<-colMeans(V,na.rm=T)
sigma<-var(V,na.rm=T)
print(mu)
print(sigma)
p<-ncol(ZZ)
for (i in 1:nrow(ZZ)){
index<-which(is.na(ZZ[i,]))
q<-length(index)
if (q>=1 & p!=q){
T<-c(ZZ[i,-index],ZZ[i,index])
Mu<-c(mu[-index],mu[index])
Si<-rbind(sigma[-index,],sigma[index,])
Si<-cbind(Si[,-index],Si[,index])
m_miss<-matrix(Mu[(p+1-q):p],nrow=q,ncol=1)
var_prod<-Si[(p+1-q):p,1:(p-q)]%*%solve(Si[1:(p-q),1:(p-q)])
shift<-matrix(T[1:(p-q)]-Mu[1:(p-q)],nrow=p-q,ncol=1)
new<-m_miss+var_prod%*%shift
ZZ[i,index]=new}else{
if (p==q){
ZZ[i,index]=mu
}}
}
ZZ<-exp(ZZ)
head(ZZ)

#########################################################################
###Generate Data Matrix from densities for unstripped data
#########################################################################

X<-matrix(nrow=nrow(data),ncol=(length(coeff)-2))
for (i in 3:length(coeff)) {
q=strsplit(coeff_names[i],":")
tmp<-do.call(rbind, q)
Z=matrix(nrow=nrow(data),ncol=length(tmp))
for (ii in 1:length(tmp)){
Z[,ii]=data[,which(colnames(data)==tmp[ii])]
}
Z=apply(Z,1,prod)
X[,i-2]=Z
}
X<-as.matrix(X)
colnames(X)<-coeff_names[-c(1,2)]
X<-log(X)
head(X)


#########################################################################
###Generate Data Matrix from densities for stripped data
#########################################################################

XX<-matrix(nrow=nrow(ZZ),ncol=(length(coeff)-2))
for (i in 3:length(coeff)) {
q=strsplit(coeff_names[i],":")
tmp<-do.call(rbind, q)
Z=matrix(nrow=nrow(ZZ),ncol=length(tmp))
for (ii in 1:length(tmp)){
Z[,ii]=ZZ[,which(colnames(ZZ)==tmp[ii])]
}
Z=apply(Z,1,prod)
XX[,i-2]=Z
}
XX<-as.matrix(XX)
colnames(XX)<-coeff_names[-c(1,2)]
XX<-log(XX)
head(XX)



#######################################################################
###Figure out MSEp by % missing level
#######################################################################
if (int==1){
XX=cbind(1,XX)
coeff<-COEFF[which(COEFF$Race==race[r_ind] & COEFF$Level==levels[b]),]
coeff=coeff[which(is.na(coeff)==FALSE)]
coeff_names<-colnames(coeff)
}
Freq_hat<-XX%*%t(coeff[3:length(coeff)])
Freq<-log(data[,1])
MSEp<-c(MSEp,sum((Freq-Freq_hat)^2)/nrow(data))
#	#x11()
#	#plot(Freq_hat~Freq,main=paste(levels[b],"MSEp=",round(MSEp[length(MSEp)],2),"g(Z) ave miss=",round(PERCENT_miss,2)))

}#### end msep loop
}#### end levels loop



MSEp<-cbind(strip_all[,c(4,5)],MSEp)
colnames(MSEp)<-c("Percent Missing","Levels","MSEp")
MSEp<-as.data.frame(MSEp)
MSEp$Levels<-levels[MSEp$Levels]
combs<-combn(levels,2)
mesh<-list()
for (i in 1:ncol(combs)){
indx<-combs[,i]
mesh<-rbind(mesh,expand.grid(m1=strip_step,m2=strip_step,label_1=indx[1],label_2=indx[2]))
}


x<-numeric()
y<-numeric()
for (i in 1:nrow(mesh)){
x<-c(x,MSEp$MSEp[which(MSEp$Levels==mesh$label_1[i] & MSEp$Per==mesh$m1[i])])
y<-c(y,MSEp$MSEp[which(MSEp$Levels==mesh$label_2[i] & MSEp$Per==mesh$m2[i])])
}

obj<-x/y
mesh$label_1<-paste("i=",mesh$label_1)
mesh$label_2<-paste("j=",mesh$label_2)
mesh$label_1<-gsub("!","",mesh$label_1)
mesh$label_2<-gsub("!","",mesh$label_2)
group<-paste(mesh$label_1,"vs",mesh$label_2)
plot_dat<-cbind(mesh,group,x,y,obj)
trellis.device(color=T)
contourplot(obj~m2*m1|group,data=plot_dat,layout=c(2,3),xlab=expression("% Missing of "*Z["#S"<="j"]),ylab=expression("% Missing of "*Z["#S"<="i"]),
scales=list(relation="free"),region=TRUE,pretty=T,cuts=10,labels=TRUE)

1/(MSEp$MSEp[which(MSEp$P==0.75 & MSEp$Levels=="4!")]/MSEp$MSEp[which(MSEp$P==1.0 & MSEp$Levels=="4!")])
1/(MSEp$MSEp[which(MSEp$P==0.75 & MSEp$Levels=="3!")]/MSEp$MSEp[which(MSEp$P==1.0 & MSEp$Levels=="3!")])
1/(MSEp$MSEp[which(MSEp$P==0.75 & MSEp$Levels=="2!")]/MSEp$MSEp[which(MSEp$P==1.0 & MSEp$Levels=="2!")])
1/(MSEp$MSEp[which(MSEp$P==0.75 & MSEp$Levels=="1!")]/MSEp$MSEp[which(MSEp$P==1.0 & MSEp$Levels=="1!")])





#######################################
#######################################
###Imputation Test for Joel
#######################################
#######################################

################################################################################
###Set directory path and load libraries
################################################################################
setwd("K:/xover/Haplotype/Imputation_SHF")
rm(list=ls(all=TRUE)) 
gc()
graphics.off()
dev.off()
library(lattice)
library(car)
library(alr3)
library(MASS)
library(faraway)
library(quantreg)

###################################################################################
###Load coefficient data to subset data matrix (will not be needed in code module)
################################################################################
filename=paste0(getwd(),"/output/2013-01-18Coefficients_OLS.csv")
filename_train=paste0(getwd(),"/output/2013-01-18Index_Train_OLS.csv")
filename_valid=paste0(getwd(),"/output/2013-01-18Index_Validation_OLS.csv")
COEFF<-read.csv(file=filename, check.names=FALSE)
COEFF_names<-gsub("`","",colnames(COEFF))
colnames(COEFF)<-COEFF_names
head(COEFF)
race<-as.character(unique(COEFF$Race))
r_ind<-12
filename=paste0(getwd(),"/nemo_est5loc_parts/",race[r_ind],".ARS.h5_component.freqs.csv")
DATA<-read.csv(file=filename, check.names=FALSE, nrows=50000,row.names=1)
data<-DATA[c(1:31),]##test top 50 values##
data_t<-DATA
density=c("5!",rep("4!",5),rep("3!",10),rep("2!",10),rep("1!",5))
names(density)<-colnames(data)
MSEp<-numeric()
b=1 #set to 4! when b=1

#######################################################################
###Strip out Data to test Imputation Approach
#######################################################################
data_s<-as.matrix(data)
for (i in 2:nrow(data)){
data_s[i,c(2:i)]<-NA
}

#######################################################################
###Start prediction loop over levels 4!, 3!, 2!, 1!
#######################################################################

levels<-unique(density)[-1]
FREQUENCY_hat<-matrix(nrow=nrow(data_s),ncol=length(levels))
PERCENT_miss<-matrix(nrow=nrow(data_s),ncol=length(levels))
colnames(FREQUENCY_hat)<-levels
rownames(FREQUENCY_hat)<-rownames(data_s)
colnames(PERCENT_miss)<-levels
rownames(PERCENT_miss)<-rownames(data_s)
coeff<-COEFF[which(COEFF$Race==race[r_ind] & COEFF$Level==levels[b]),]
coeff=coeff[which(is.na(coeff)==FALSE)]
coeff_names<-colnames(coeff)
int<-0
if (any(coeff_names=="(Intercept)")){
int<-1
coeff=coeff[-which(coeff_names=="(Intercept)")]
coeff_names=coeff_names[-which(coeff_names=="(Intercept)")]
}


#######################################################################
###This is where the code starts for imputation of missing values using MVN density assumption
#######################################################################
X<-log(data[,-1])
XX<-log(data_s[,-1])
XT<-log(data_t[,-1])
XXX<-XX
V<-XT
mu<-colMeans(V,na.rm=TRUE)
sigma<-var(V,na.rm=TRUE)
print(mu)
print(sigma)
p<-ncol(XX)
for (i in 1:nrow(XX)){
index<-which(is.na(XX[i,]))
q<-length(index)
if (q>=1 & p!=q){
T<-c(XX[i,-index],XX[i,index])
Mu<-c(mu[-index],mu[index])
Si<-rbind(sigma[-index,],sigma[index,])
Si<-cbind(Si[,-index],Si[,index])
m_miss<-matrix(Mu[(p+1-q):p],nrow=q,ncol=1)
var_prod<-Si[(p+1-q):p,1:(p-q)]%*%solve(Si[1:(p-q),1:(p-q)])
shift<-matrix(T[1:(p-q)]-Mu[1:(p-q)],nrow=p-q,ncol=1)
new<-m_miss+var_prod%*%shift
XX[i,index]=new}else{
if (p==q){
XX[i,index]=mu
}}
}
head(XX)


#########################################################################
###Generate Data Matrix from densities for unstripped data
#########################################################################
data[,-1]<-exp(XX)
X<-matrix(nrow=nrow(data),ncol=(length(coeff)-2))
for (i in 3:length(coeff)) {
q=strsplit(coeff_names[i],":")
tmp<-do.call(rbind, q)
Z=matrix(nrow=nrow(data),ncol=length(tmp))
for (ii in 1:length(tmp)){
Z[,ii]=data[,which(colnames(data)==tmp[ii])]
}
Z=apply(Z,1,prod)
X[,i-2]=Z
}
X<-as.matrix(X)
colnames(X)<-coeff_names[-c(1,2)]
X<-log(X)
head(X)

#######################################################################
###Use Regression Model Results to get Y versus Y_hat for data
#######################################################################

if (int==1){
X=cbind(1,X)
coeff<-COEFF[which(COEFF$Race==race[r_ind] & COEFF$Level==levels[b]),]
coeff=coeff[which(is.na(coeff)==FALSE)]
coeff_names<-colnames(coeff)
}
Freq_hat<-X%*%t(coeff[3:length(coeff)])
Freq<-log(data[,1])
MSEp<-sum((Freq-Freq_hat)^2)/nrow(data)
plot(Freq_hat~Freq,main=list(paste0(race[r_ind]), cex=1),xlab=expression("log("*Y*")"),ylab=expression("log("*hat(Y)*")"))




par(mfrow=c(1,2),mar=c(3,3.25,3,1)+0.1,mgp=c(2,1,0))
plot(Freq_hat~Freq,main=list(paste0("Log(Y):",race[r_ind]), cex=1),xlab=list(expression(Y),cex=1),ylab=list(expression(hat(Y)),cex=1))
plot(exp(Freq_hat)~exp(Freq),main=list(paste0("Y:",race[r_ind]), cex=1),xlab=list(expression(Y),cex=1),ylab=list(expression(hat(Y)),cex=1))
MSEp<-sum((exp(Freq)-exp(Freq_hat))^2)/nrow(data)
print(MSEp)
tmp_out<-cbind(Freq,Freq_hat,exp(Freq),exp(Freq_hat))
rownames(tmp_out)<-rownames(data)
colnames(tmp_out)<-c("log(Y)","log(Y_hat)","Y","Y_hat")
head(tmp_out)
write.csv(tmp_out, file=paste0(getwd(),"/validaiton_joel.csv"),row.names=TRUE,col.names=TRUE)





