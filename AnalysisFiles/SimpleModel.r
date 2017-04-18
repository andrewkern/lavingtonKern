#!/usr/bin/Rscript 


load("NoHets")

Newdata<-NoHets


tempdata<-NULL
tempdata<-as.data.frame(tempdata)
for(m in 1:(ncol(Newdata)-5)){

	
	placeholder<-1#for keeping track of the output column
	
	fit.model<-anova(lm(Newdata[[m+5]]~Sex + In2Lt + In3RMo + Line%in%In2Lt, data=Newdata))#fit the 'good' model

	for (i in 1:5){
		tempdata[m,placeholder]<-fit.model$'Sum Sq'[i]#drop one effect from the model
		placeholder<-placeholder+1
	}#to get the sum of squares
		
	for (i in 1:4){
		tempdata[m,placeholder]<-fit.model$'F value'[i]#drop one effect from the model
		placeholder<-placeholder+1
		tempdata[m,placeholder]<-fit.model$'Pr(>F)'[i]#drop one effect from the model
		placeholder<-placeholder+1
		
	}#to get the F and p values for each effect 
#calculate the R^2 of the model and eta^2 for each effect

for(i in 1:4){
	tempdata[m,placeholder]<-tempdata[m,i]/sum(tempdata[m,(1:5)])
	placeholder<-placeholder+1
}#to get the eta^2 's

tempdata[m,placeholder]<-1-tempdata[m,5]/sum(tempdata[m,(1:5)])#to get the R^2 of the model
placeholder<-placeholder+1

tempdata[m,placeholder]<-colnames(Newdata)[m+5]#get the probe name from current column of the Newfile data frame 

}
colnames(tempdata)[1]<-"ss_Sex"
colnames(tempdata)[2]<-"ss_In2Lt"
colnames(tempdata)[3]<-"ss_In3RMo"
colnames(tempdata)[4]<-"ss_In2Lt.Line"
colnames(tempdata)[5]<-"ss_Residuals"
colnames(tempdata)[6]<-"F_Sex"
colnames(tempdata)[7]<-"p_Sex"
colnames(tempdata)[8]<-"F_In2Lt"
colnames(tempdata)[9]<-"p_In2Lt"
colnames(tempdata)[10]<-"F_In3RMo"
colnames(tempdata)[11]<-"p_In3RMo"
colnames(tempdata)[12]<-"F_In2Lt.Line"
colnames(tempdata)[13]<-"p_In2Lt.Line"
colnames(tempdata)[14]<-"eta_Sex"
colnames(tempdata)[15]<-"eta_In2Lt"
colnames(tempdata)[16]<-"eta_In3RMo"
colnames(tempdata)[17]<-"eta_In2Lt.Line"
colnames(tempdata)[18]<-"R^2"
colnames(tempdata)[19]<-"Probe"


write.table(tempdata, file="/san/personal/lavington/SimpleModel.txt")



