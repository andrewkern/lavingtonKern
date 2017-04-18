#!/bin/bash
# first line of snpReports is a header; 
#read in line two as the first 'OldSnps' array
#Start loop on line 3: 
		#run correlation between NewSnps and OldSnps, 
		#then rename NewSnps as OldSnps,
		#output data line, 
		#repeat to input EOF
#command line variables: $1: Inversion pattern $2 output filename 
#Top=State1 (Major allele) OldSnps
#Left=State1 (Major allele) NewSnps
#Input file (snpReports) header:
#	1	2			3		4		5		6	7		8:n-1	n
#	loc	sampleSize	nStates	state1	state2	MAF	RAL-100	...	RAL-93

awk ' { if($NR ~ 1){
			print $0;
		}
		else if($NR ~ 2){
			split($i,OldSnps);
		}
		else if($1 ~ /^[0-9]/ ){
			n=split($i,NewSnps);
		 	tl=0;
			tr=0;
			bl=0;
			br=0;
			total=0;
			trash=0;
		
			for (x=7;x<=211;x++){
				if (OldSnps[4] ~ OldSnps[x] && NewSnps[4] ~ NewSnps[x]){
					++tl;
				} else if(OldSnps[4] ~  OldSnps[x] && NewSnps[5] ~ NewSnps[x]){
					++tr;
				} else if (OldSnps[5] ~  OldSnps[x] && NewSnps[4] ~ NewSnps[x]){
					++bl;
				} else if (OldSnps[5] ~  OldSnps[x] && NewSnps[5] ~ NewSnps[x]){
					++br;
				} else{
					++trash;
				}
			
			}
		}
			left=tl+bl;
			right=tr+br;
			top=tl+tr;
			bottom=bl+br;
			total=tl+tr+bl+br;
			dist=NewSnps[1]-OldSnps[1];
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
			
			print NewSnps[1],total,rsquare,NewSnps[6],dist,tl,tr,bl,br,left,top,chisquare >> "/san/personal/lavington/corr2LSnps";
			
			delete OldSnps;
			
			for(y=1;y<=n;y++){
			
				OldSnps[y]=NewSnps[y];
				delete NewSnps[y];
			}
			}' /san/data/melanogaster/drosophilaNexus/dgrp/snpReports/allDGRP.chr2L.fa.snpReport