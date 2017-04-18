#!/usr/bin/Rscript --vanilla

#Input file names
#NoHets_0_Sex_F
#NoHets_1_In2Lt_F
#NoHets_2_In3RMo_F
#NoHets_3_Line.In2Lt_F



NoHetData<-read.table("SimpleModel.txt", header=TRUE)

#initialize the target matrix
Percs<-matrix(data=NA, nrow=18952, ncol=4)

#-----------------------------

Fs<-read.table("/san/personal/lavington/SimpleModel/SimpleModel_0",header=FALSE)#
Fs<-as.matrix(Fs)
for( i in 1:ncol(Fs)){
	tempfunc<-ecdf(Fs[,i])
	Percs[i,1]<-tempfunc(NoHetData[i,6])

}#Sex

Fs<-read.table("/san/personal/lavington/SimpleModel/SimpleModel_1")
Fs<-as.matrix(Fs)
for( i in 1:ncol(Fs)){
	tempfunc<-ecdf(Fs[,i])
	Percs[i,2]<-tempfunc(NoHetData[i,8])

}#In2Lt

Fs<-read.table("/san/personal/lavington/SimpleModel/SimpleModel_2")
Fs<-as.matrix(Fs)
for( i in 1:ncol(Fs)){
	tempfunc<-ecdf(Fs[,i])
	Percs[i,3]<-tempfunc(NoHetData[i,10])

}#In3RMo

Fs<-read.table("/san/personal/lavington/SimpleModel/SimpleModel_3")
Fs<-as.matrix(Fs)
for( i in 1:ncol(Fs)){
	tempfunc<-ecdf(Fs[,i])
	Percs[i,4]<-tempfunc(NoHetData[i,12])

}#Sex*In2Lt

rm(Fs)

write.table(Percs, file="/san/personal/lavington/SimpleModel/SMpercentiles")
