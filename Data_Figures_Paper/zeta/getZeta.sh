

FILE=$1

if [ -z "$FILE" ]
then
	echo "Defin potential file"
	exit
fi



cat $FILE | awk 'BEGIN{
	phi1=-1
	phi2=-1
	zeta=-1
	Rbest=0
	Ebest=0
}{
	if(zeta!=$3){
		if (zeta != -1) printf("%5d %5d %12.8f %12.8f %12.8f\n", phi1,phi2,zeta, Rbest,Ebest)
		phi1=$1
		phi2=$2
		zeta=$3
		Rbest=$4
		Ebest=$5
	}

	if($5<Ebest){
		Rbest=$4
		Ebest=$5
	} 

}END{
}'
