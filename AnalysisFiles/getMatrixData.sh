
#datafile1 needs to be the large file to extract data from

#datafile2 needs to ba a list of matrix column/row names to extract columns/rows from the main data file 

#specifically designed for a matrix of between-probe r-square values for Ayroles et al 2009 expression data. 

#Main data file was output from R, so probe names begin with X, seperated by a single space, and surounded by quotes
#NA was used to fill the first field of the first row (header row), without quotes
#first field of each row is probe label as described for header row
#data values, including NA's are seperated by a single space, no quotes

awk 'FILENAME=="AyrolesMatrix"{
	if(NR==1){
		z=split($0,names," "); #split the header row into an array, keep array size as z
	}
	else{
		split($0,temp," ");
		for(i;i<=z;i++){
			mat[$1,names[i]]=$i;
		}

	}#
}
FILENAME=="2LtNames"{
	IAL[FNR]=$1;
	LAI[FNR]=$1;
	sizeIAL=FNR;
}
{
	
	for(j=1;j<=sizeIAL;j++){
		for(k=j+1;j<=sizeIAL;k++){
			print mat[IAL[j],LAI[k]] >> "out";
		}
	}	
}' AyrolesMatrix 2LtNames