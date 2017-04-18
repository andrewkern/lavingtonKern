#!/usr/bin/Rscript --vanilla

args<-commandArgs(trailingOnly=TRUE)
iteration<-as.numeric(args[1])
out.fileF<-paste("/san/personal/lavington/SimpleModel/SimpleModel_F_",iteration, ".txt", sep="")
out.fileP<-paste("/san/personal/lavington/SimpleModel/SimpleModel_p_",iteration, ".txt", sep="")

#read in data table
load("NoHets") #data with polymorphic lines removed (for both In2Lt and In3RMo)

Newdata<-NoHets #adapt old code for new data

#find number of columns and define number of not probe columns probes explicitly
numCOLS<-ncol(Newdata)


numNOTprobes<-5 #or find them, they will be the NOT numeric columns

#set number of replications, permutations, explicitly

m<-1
#permute the Inversion state for each line, for each inversion-Inversion state is independent of line, one line can have multiple states

RandomLt<-sample(Newdata$In2Lt)
RandomMO<-sample(Newdata$In3RMo)


tempF<-matrix(data=NA, ncol=(numCOLS-numNOTprobes), nrow=4)
tempP<-matrix(data=NA, ncol=(numCOLS-numNOTprobes), nrow=4)

#colnames(allFs)<-probeNAMES #;add column name as we go in the for() loop

#RUN REPLCATES OF RANDOMIZED DATA; Randomize the inversion state of a line and retain relationships between samples in a line and expression in a sample

	
	for(m in 1:(numCOLS-numNOTprobes)){
	
		tempAov <- anova(lm(Newdata[[m+5]]~Sex + RandomLt + RandomMO + Line%in%RandomLt,data=Newdata))
		
	for (j in 1:4){
		tempF[j,m]<-tempAov$'F value'[j]
		tempP[j,m]<-tempAov$'Pr(>F)'[j]
		}#add the effect F and p values to the dataframe
  
	}# ANOVA test for SINGLE strREPLICATE RANDOM inversion state, two inversions, all probes



write.table(tempF, file=out.fileF, append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(tempP, file=out.fileP, append=TRUE, row.names=FALSE, col.names=FALSE)
#-----------#Output rows, from 1 to 12:
#1:Sex
#2:In2Lt
#3:In3RMo
#4:Line in In2Lt