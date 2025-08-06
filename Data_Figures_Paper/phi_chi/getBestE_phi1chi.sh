#!/bin/bash



FILE=$1
HAND=$2

if [ HAND == "" ]
then 
	HAND=1
fi

for i in {0..2}
do
	cat $FILE | awk '{print ($1+(-1.0)*'$HAND'*$2+230)%360-230,($2+10)%20 -30 +'$i'*20,$3,$4,$5}' | awk 'BEGIN{phi1=-1;chi=-12345;eBest=0}{
		if($2 != chi) {
			if(chi != -12345) printf("%+4d %+4d %12.8f %12.8f %12.8e\n",phi1 ,chi, zBest,rBest,eBest)
			phi1=$1
			chi=$2
			zBest=$3
			rBest=$4
			eBest=$5
		}
		else if ($5<eBest){
			zBest=$3
			rBest=$4
			eBest=$5
		}
	}' 
done > tmp

for i in {-230..130..1} ; do  
	#~ echo $i ;
	cat tmp | awk '{if($1=='$i') print}' | sort -n -k 2 ;
	echo ;
done

rm tmp
