#!/bin/bash



FILE=$1
HAND=$2

if [ HAND == "" ]
then 
	HAND=1
fi

for i in {0..17}
do
	cat $FILE | awk '{print ($1+'$i'*20) % 360 ,($2+'$HAND'*'$i'*20+360) % 360,$3,$4,$5}' | awk 'BEGIN{phi1=-1;phi2=-1;eBest=0}{
		if($2 != phi2) {
			if(phi2 != -1) printf("%+4d %+4d %12.8f %12.8f %12.8e\n",(phi1+180) %360 -180 ,(phi2+180) %360 -180, zBest,rBest,eBest)
			phi1=$1
			phi2=$2
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

for i in {-180..180..1} ; do  echo $i ; cat tmp | awk '{if($1=='$i') print}' | sort -n -k 2 ;echo ; done

rm tmp
