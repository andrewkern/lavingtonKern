#!/usr/bin/awk -f


# awk [ 2Lt , 3RMo , 2RNs , etc. ] to get a list of inversions: awk '{ if ($1 ~ "In") print $1}' snpReportInvState.chr[n] where [n] = 2L, 2R, 3L, 3R, or X
#command line variables: Inver="Inversion pattern" outfile="output filename" 

#SNP and Inversion state input file (snpReport...) header:
#	1	2			3		4		5		6	7		8:n-1	n
#	loc	sampleSize	nStates	state1	state2	MAF	RAL-100	...	RAL-93



#Execution Example:
# 	(to get the LD between each SNP on chr2L with In(2L)t)
#where the data file with chr2L SNP data and Inversion states for the DGRP is named: "snpReportInvState.chr2L" 
#and this awk script file is named: "SNPandINVcorr.sh"
#The command should look like:
#SNPandINVcorr.sh -v inv="2Lt" -v outfile="YourOutputFilename" snpReportInvState.chr2L

{ if($1 ~ inv){
			split($i,inversions)	
		}else if($1 ~ /^[0-9]/ ){
			split($i,snps);
		 	tl=0;
			tr=0;
			bl=0;
			br=0;
			total=0;
			trash=0;
		
			for (x=7;x<=211;x++){
				if (inversions[x] ~ "ST" && snps[x] ~ snps[4]){
					++tl;
				} else if(inversions[x] ~ "ST" && snps[x] ~ snps[5]){
					++tr;
				} else if (inversions[x] ~ "INV" && snps[x] ~ snps[4]){
					++bl;
				} else if (inversions[x] ~ "INV" && snps[x] ~ snps[5]){
					++br;
				} else{
					++trash;
				}
			}
			
			left=tl+bl;
			right=tr+br;
			top=tl+tr;
			bottom=bl+br;
			total=tl+tr+bl+br;
			if(tl > 0 && tr > 0 && bl >0 && br > 0){
				rsquare=((((tl/total)*(br/total))-((tr/total)*(bl/total)))^2)/((left/total)*(right/total)*(top/total)*(bottom/total));
				chisquare=rsquare*total;
			} else if(top == 0){
				rsquare = "0";
				chisquare=0;
			} else if(bottom == 0){
				rsquare = "0";
				chisquare=0;
			} else if(left == 0){
				rsquare = "0";
				chisquare=0;
			} else if(right == 0){
				rsquare = "0";
				chisquare=0;
			} else if((((bl+tr)>0)&&((tl+br)==0)) || (((bl+tr) == 0)&&((tl+br)>0))){
				rsquare=((((tl/total)*(br/total))-((tr/total)*(bl/total)))^2)/((left/total)*(right/total)*(top/total)*(bottom/total));
				chisquare=rsquare*total;
			} else{
				rsquare=((((tl/total)*(br/total))-((tr/total)*(bl/total)))^2)/((left/total)*(right/total)*(top/total)*(bottom/total));
				chisquare=rsquare*total;
			}
			print NR,snps[1],tl,tr,bl,br,left,top,total,rsquare,chisquare >> outfile;
		}
}