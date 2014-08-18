########File for Synthetic Haplotype Estimation Based Upon 	####################
########Lower level marginal distributions.  This is based 	####################
########upon a regression approach.                        	####################
########Code was updated to use independent densities		####################
########Creation date:  12/07/2012                         	####################
########Author:  Mark Albrecht, Ansu Chatterjee            	####################
########Revision:  4.0	(added an intercept)		     	####################



################################################################################
###Set directory path and load libraries
################################################################################
##setwd("/netware/grpshare/Haplotype/Imputation_SHF")
setwd("K:/xover/Haplotype/Imputation_SHF")
rm(list=ls(all=TRUE)) 
gc()
graphics.off()
#	#dev.off()
#	#library(lattice)
#	#library(car)
#	#library(alr3)
#	#library(MASS)
#	#library(faraway)
#	#library(quantreg)

Regression=c("OLS")###  Specifiy one of :  c("OLS","IRWLS","QR")
mod_select=c("both") ###  Specifiy one of :  c("backward","both","forward")
race=c("AAFA","AFA","AFB","AINDI","AISC","ALANAM","AMIND","API","CARB","CARHIS","CARIBI","CAU","DEC","FILII","HAWI","HIS","JAPI","KORI","MENAFC",
"MLT","MSWHIS","NAM","NAMER","NCHI","OTH","SCAHIS","SCAMB","SCSEAI","UNK","VIET")
race=c("AFA","CAU","HIS","API")
sample_sizes<-c(rep(50000,length(race)))
r_ind<-1
for (r_ind in 1:length(race)){ 
print(race[r_ind])

###################################################################################
###Load data
################################################################################
filename=paste0(getwd(),"/nemo_est5loc_parts/",race[r_ind],".ARS.h5_component.freqs.csv")
DATA<-read.csv(file=filename, check.names=FALSE, row.names=1,nrows=sample_sizes[r_ind])
colnames(DATA)[1]<-"Freq"
data<-DATA
################################################################################
###format data to generate model matrix having all terms with order<=5
################################################################################
density=c("5!",rep("4!",5),rep("3!",10),rep("2!",10),rep("1!",5))
names(density)<-colnames(data)
subset<-density[-1]
ord<-c(rep(4,5),rep(3,10),rep(2,10),rep(1,5))
content=matrix(0,nrow=5,ncol=30)
rownames(content)<-c("A","C","B","DRB1","DQB1")
colnames(content)<-names(subset)
content[1,c(grep("A",names(subset)))]=1
content[2,c(grep("C",names(subset)))]=1
content[3,c(grep("~B",names(subset)))]=1
content[3,c(grep("B~",names(subset)))]=1
content[3,28]=1
content[4,c(grep("DRB1",names(subset)))]=1
content[5,c(grep("DQB1",names(subset)))]=1
print(content)




order_1<-5
comb_1<-1
OUT=matrix(nrow=5)
DEN=numeric()
j=1
for (j in 1:5){
comb<-combn(1:30,j)
o<-ord[comb]
o<-matrix(o,nrow=j,ncol=ncol(comb))
den<-apply(o,2,max)
o<-apply(o,2,sum)
comb<-comb[,which(o<=order_1)]
comb<-matrix(comb,nrow=j)
den<-den[which(o<=order_1)]
cont<-matrix(nrow=5,ncol=ncol(comb))
rownames(cont)<-rownames(content)
ii=1
for (ii in 1:5){
tmp<-matrix(content[ii,comb],nrow=j,ncol=ncol(comb))
cont[ii,]<-apply(tmp,2,sum)
}
comb=comb[,which(apply(cont,2,max)<=comb_1)]
comb<-matrix(comb,nrow=j)
den=den[which(apply(cont,2,max)<=comb_1)]
values<-matrix(names(subset)[comb],nrow=nrow(comb),ncol=ncol(comb))
tt=matrix(NA,nrow=5,ncol=ncol(comb))
tt[1:nrow(comb),]=comb
OUT<-cbind(OUT,tt)
DEN<-c(DEN,den)
}
OUT<-OUT[,-1]
NAMES<-matrix(names(subset)[OUT],nrow=nrow(OUT),ncol=ncol(OUT))

sub_dat<-data[,-1]
X<-matrix(nrow=nrow(sub_dat))
j=1
for (j in 1:ncol(OUT)){
tmp<-as.matrix(sub_dat[,OUT[!is.na(OUT[,j]),j]])
tmp<-apply(tmp,1,prod)
tmp<-as.matrix(tmp)
colnames(tmp)<-paste(NAMES[!is.na(NAMES[,j]),j],collapse=":")
X=cbind(X,tmp)
}
X<-X[,-1]







################################################################################
###Partition and log-transform the data
################################################################################
DATA<-as.data.frame(cbind(DATA[,1],X))
DEN<-c(5,DEN)###add response density
colnames(DATA)[1]<-"Freq"
DATA<-log(DATA)
index_train<-sample(1:sample_sizes[r_ind],sample_sizes[r_ind]*0.8)
index_validation<-1:sample_sizes[r_ind]
index_validation<-index_validation[-index_train]
data<-DATA[index_train,]#code for top haplotypes assuming pre-sorted file


################################################################################
###Regression Modeling
###Partition by class of distributions
###4!,3!,2!,1!
################################################################################
y_hat<-matrix(NA,nrow(data),4)
print(race[r_ind])
UB.Fit<-lm(Freq~.,data=data)
model_coefficients<-matrix(NA,4,length(UB.Fit$coef))
colnames(model_coefficients)<-names(UB.Fit$coef)
identifier<-names(UB.Fit$coef)

loop<-c(4,3,2,1)
b<-1
i<-1
for (i in loop){  

#################################################################
##### Add iterative reweighting while loop for WLS estimation####
#################################################################


if (Regression=="OLS") {###Select correct regression method:  c("OLS","IRWLS","QR")
print(race[r_ind])
UB.Fit<-lm(Freq~.,data=data[,c(1,which(DEN<=i))])
LB.Fit<-lm(Freq~1,data=data[,c(1,which(DEN<=i))])
Step.Linear.Fit=step(LB.Fit,scope=list(upper=formula(UB.Fit),lower=formula(LB.Fit)),direction = mod_select,k=log(nrow(data)))
}else
if (Regression=="IRWLS"){
obj<-1 #starting value for objective function
wts<-rep(1,nrow(data))  #starting weights
q<-1
MSE<-1
while (obj>=0.0001) {
UB.Fit<-lm(Freq~.,data=data[,c(1,which(DEN<=i))],weights=wts)
LB.Fit<-lm(Freq~1,data=data[,c(1,which(DEN<=i))],weights=wts)
Step.Linear.Fit=step(LB.Fit,scope=list(upper=formula(UB.Fit),lower=formula(LB.Fit)),direction = mod_select,k=log(nrow(data)))
y_hat[,b]=predict(Step.Linear.Fit)
yy<-as.numeric((residuals(Step.Linear.Fit)^2))
xx<-as.numeric(y_hat[,b])
lo.fit<-loess(log(yy)~xx,degree=1,span=0.5)
var_est<-predict(lo.fit,data.frame(xx=y_hat[,b]))
wts<-1/exp(var_est)
tmp<-sum(residuals(Step.Linear.Fit)^2)/Step.Linear.Fit$d
tmp<-tmp[length(tmp)]
MSE<-cbind(MSE,tmp)
obj<-(MSE[q+1]-MSE[q])^2 ###Stopping criteria
q<-q+1
labels<-c("4!","3!","2!","1!")
plot(log(yy)~xx,main=paste(race[r_ind],labels[b],"MSE=",round(MSE[q],digits=2)),ylab="Log(Residual^2)",xlab="Predicted Log(Frequency)")
points(var_est~y_hat[,b],pch=20)
} 


}else
if (Regression=="QR") {

print(race[r_ind])
UB.Fit<-lm(Freq~.,data=data[,c(1,which(DEN<=i))])
LB.Fit<-lm(Freq~1,data=data[,c(1,which(DEN<=i))])
Step.Linear.Fit=step(LB.Fit,scope=list(upper=formula(UB.Fit),lower=formula(LB.Fit)),direction = mod_select,k=log(nrow(data)))

UB.Fit<-rq(formula(Step.Linear.Fit),data=data[,c(1,which(DEN<=i))])
Step.Linear.Fit=step(UB.Fit,scope=list(upper=formula(UB.Fit),lower=Freq~1),direction = mod_select,k=log(nrow(data)))
}

assign(paste("m", b, sep = ""), Step.Linear.Fit) 
	ii=1
	for (ii in 1:length(Step.Linear.Fit$coef)) {
	indx<-which(colnames(model_coefficients)==names(Step.Linear.Fit$coef[ii]))
	model_coefficients[b,indx]=Step.Linear.Fit$coef[ii]
	}

b<-b+1
}


print(model_coefficients)

levels<-rbind("4!","3!","2!","1!")
thisCoefficients=cbind(race[r_ind],levels,model_coefficients)
colnames(thisCoefficients)=c("Race","Level",colnames(model_coefficients))
thisIndexValidation=matrix(index_validation,nrow=1)
rownames(thisIndexValidation)<-race[r_ind]
thisIndexTrain=matrix(index_train,nrow=1)
rownames(thisIndexTrain)<-race[r_ind]
if (!exists("allCoefficients")){
allCoefficients=thisCoefficients
IndexValidation=thisIndexValidation
IndexTrain=thisIndexTrain
} else {
allCoefficients=rbind(allCoefficients,thisCoefficients)
IndexValidation=rbind(IndexValidation,thisIndexValidation)
IndexTrain=rbind(IndexTrain,thisIndexTrain)
}



################################################################################
###Validation Loop
################################################################################
validation<-DATA[index_validation,]#code for top haplotypes assuming pre-sorted file
levels<-c("4!","3!","2!","1!")
i=1
msep=as.numeric()

plotit<-FALSE  ###change to false to supress plotting
if (plotit==TRUE){
x11()
par(mfrow=c(2,2))}
for (i in 1:length(levels)){
plot_data<-validation
tmpCmd = paste("y_pred=predict(m",i,",plot_data)", sep = "")     
eval(parse(text = tmpCmd)) 
msep=c(msep,sum((y_pred-validation[,1])^2)/length(y_pred))

if (plotit==TRUE){
plot(y_pred,validation[,1])
handle<-lm(validation[,1]~y_pred-1)
abline(handle,lwd=2)
title(main=paste0("Validation Y vs Y_hat",levels[i]," r=",race[r_ind]))
handle<-lm(validation[,1]~y_pred-1)
}

tmpCmd = paste("y_pred",i,"=y_pred", sep = "")     
eval(parse(text = tmpCmd)) 
}


if (!exists("all_y_pred1")){
for (i in 1:length(levels)){
tmpCmd = paste0("all_y_pred",i,"=y_pred",i)
eval(parse(text = tmpCmd)) 
}
all_validation=validation[,1]
all_msep=msep
} else {
for (i in 1:length(levels)){
tmpCmd = paste0("all_y_pred",i,"=","cbind(all_y_pred",i,",y_pred",i,")")
eval(parse(text = tmpCmd)) 
}
all_validation=cbind(all_validation,validation[,1])
all_msep=cbind(all_msep,msep)
}




} #### end race loop




for (y in 1:4) {
name<-4:1
pdf(paste0(getwd(),"/output/",Sys.Date(),"validation_plots_",Regression,name[y],".pdf"))


par(mfrow=c(5,6),mar=c(1,1,1,1)+0.1,mgp=c(3,0.3,0),cex.axis=.7)

for (r in 1:ncol(all_msep)){

if (y==1){
y_pred=all_y_pred1[,r]
} else if (y==2){
y_pred=all_y_pred2[,r]
} else if (y==3){
y_pred=all_y_pred3[,r]
} else if (y==4){
y_pred=all_y_pred4[,r]
}

msep=signif(all_msep[y,r],digits=3)
validation=all_validation[,r]

plot(y_pred,validation)
left=min(y_pred)
right=max(y_pred)
top=max(validation)
bottom=min(validation)
text(left-0.1*(right-left),top-0.05*(top-bottom),paste0("MSEp=",msep),pos=4,col="red",cex=1,font=2)
#	#text(left-0.1*(right-left),top-0.25*(top-bottom),paste0("MSE=",mse),pos=4,col="red",cex=1,font=2)
#	#handle<-lm(validation~y_pred)
#	#abline(handle,lwd=2)
title(main=paste(race[r],levels[y]))


}# end plot race loop

dev.off()


}# end plot yfactor loop
colnames(all_msep)<-race
rownames(all_msep)<-levels
write.csv(allCoefficients, file=paste0(getwd(),paste0("/output/",Sys.Date(),"Coefficients_",Regression,".csv")),row.names=FALSE)
write.csv(all_msep, file=paste0(getwd(),paste0("/output/",Sys.Date(),"MSE_predicted_",Regression,".csv")),row.names=TRUE)
write.csv(IndexTrain, file=paste0(getwd(),paste0("/output/",Sys.Date(),"Index_Train_",Regression,".csv")),row.names=TRUE)
write.csv(IndexValidation, file=paste0(getwd(),paste0("/output/",Sys.Date(),"Index_Validation_",Regression,".csv")),row.names=TRUE)
