chr2L_56<-matrix(nrow=100000,ncol=3);

for(j in 1:100000){
	y <-all_2L[sort(sample(nrow(all_2L),56,),decreasing=FALSE),]#sample n rows from data.frame(all_2Lt) and sort from first to last:rows of all_2Lt are all genes in order from proximal to distal
	
	for(i in 2:nrow(y)){y$DistBetween[i]<-y$r5.49Start[i]-y$r5.49Stop[(i-1)]}#calculate distance from gene in row i to gene in previous row (nearest neighbor on the left side)
	for(i in 2:nrow(y)){if(y$DistBetween[i] < 0){y$DistBetween[i]<-0;}}	#overlapping genes will result in a negative distance: find and set to 0
	chr2L_56[j,1]<-mean(y$DistBetween[(2:nrow(y))])#average(mean) of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr2L_56[j,2]<-sd(y$DistBetween[(2:nrow(y))])#standard deviation of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr2L_56[j,3]<-chr2L_56[j,2]/chr2L_56[j,1]#coefficient of variation (CV) of distance between genes in this bootstrap, ignores row 1 as uninformative implicitly from previous calculations
}		

chr2R_70<-matrix(nrow=100000,ncol=3);

for(j in 1:100000){
	y <-all_2R[sort(sample(nrow(all_2R),70,),decreasing=FALSE),]#sample n rows from data.frame(all_2Lt) and sort from first to last:rows of all_2Lt are all genes in order from proximal to distal
	
	for(i in 2:nrow(y)){y$DistBetween[i]<-y$r5.49Start[i]-y$r5.49Stop[(i-1)]}#calculate distance from gene in row i to gene in previous row (nearest neighbor on the left side)
	for(i in 2:nrow(y)){if(y$DistBetween[i] < 0){y$DistBetween[i]<-0;}}	#overlapping genes will result in a negative distance: find and set to 0
	chr2R_70[j,1]<-mean(y$DistBetween[(2:nrow(y))])#average(mean) of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr2R_70[j,2]<-sd(y$DistBetween[(2:nrow(y))])#standard deviation of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr2R_70[j,3]<-chr2R_70[j,2]/chr2R_70[j,1]#coefficient of variation (CV) of distance between genes in this bootstrap, ignores row 1 as uninformative implicitly from previous calculations
}

chr3L_56<-matrix(nrow=100000,ncol=3);

for(j in 1:100000){
	y <-all_3L[sort(sample(nrow(all_3L),56,),decreasing=FALSE),]#sample n rows from data.frame(all_2Lt) and sort from first to last:rows of all_2Lt are all genes in order from proximal to distal
	
	for(i in 2:nrow(y)){y$DistBetween[i]<-y$r5.49Start[i]-y$r5.49Stop[(i-1)]}#calculate distance from gene in row i to gene in previous row (nearest neighbor on the left side)
	for(i in 2:nrow(y)){if(y$DistBetween[i] < 0){y$DistBetween[i]<-0;}}	#overlapping genes will result in a negative distance: find and set to 0
	chr3L_56[j,1]<-mean(y$DistBetween[(2:nrow(y))])#average(mean) of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr3L_56[j,2]<-sd(y$DistBetween[(2:nrow(y))])#standard deviation of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr3L_56[j,3]<-chr3L_56[j,2]/chr3L_56[j,1]#coefficient of variation (CV) of distance between genes in this bootstrap, ignores row 1 as uninformative implicitly from previous calculations
}

chr3R_271<-matrix(nrow=100000,ncol=3);

for(j in 1:100000){
	y <-all_3R[sort(sample(nrow(all_3R),271,),decreasing=FALSE),]#sample n rows from data.frame(all_2Lt) and sort from first to last:rows of all_2Lt are all genes in order from proximal to distal
	
	for(i in 2:nrow(y)){y$DistBetween[i]<-y$r5.49Start[i]-y$r5.49Stop[(i-1)]}#calculate distance from gene in row i to gene in previous row (nearest neighbor on the left side)
	for(i in 2:nrow(y)){if(y$DistBetween[i] < 0){y$DistBetween[i]<-0;}}	#overlapping genes will result in a negative distance: find and set to 0
	chr3R_271[j,1]<-mean(y$DistBetween[(2:nrow(y))])#average(mean) of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr3R_271[j,2]<-sd(y$DistBetween[(2:nrow(y))])#standard deviation of distances between genes in this bootstrap, ignore row 1 as uninformative
	chr3R_271[j,3]<-chr3R_271[j,2]/chr3R_271[j,1]#coefficient of variation (CV) of distance between genes in this bootstrap, ignores row 1 as uninformative implicitly from previous calculations
}

chrX_32<-matrix(nrow=100000,ncol=3);

for(j in 1:100000){
	y <-all_X[sort(sample(nrow(all_X),32,),decreasing=FALSE),]#sample n rows from data.frame(all_2Lt) and sort from first to last:rows of all_2Lt are all genes in order from proximal to distal
	
	for(i in 2:nrow(y)){y$DistBetween[i]<-y$r5.49Start[i]-y$r5.49Stop[(i-1)]}#calculate distance from gene in row i to gene in previous row (nearest neighbor on the left side)
	for(i in 2:nrow(y)){if(y$DistBetween[i] < 0){y$DistBetween[i]<-0;}}	#overlapping genes will result in a negative distance: find and set to 0
	chrX_32[j,1]<-mean(y$DistBetween[(2:nrow(y))])#average(mean) of distances between genes in this bootstrap, ignore row 1 as uninformative
	chrX_32[j,2]<-sd(y$DistBetween[(2:nrow(y))])#standard deviation of distances between genes in this bootstrap, ignore row 1 as uninformative
	chrX_32[j,3]<-chrX_32[j,2]/chrX_32[j,1]#coefficient of variation (CV) of distance between genes in this bootstrap, ignores row 1 as uninformative implicitly from previous calculations
}