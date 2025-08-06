

FILE=$1



if [ -z "$FILE" ]
then
	echo "Defin potential file"
	exit
fi



sort -g -k 3 $FILE | awk 'BEGIN{
	phi1=-1
	phi2=-1
	zeta=12345
	Rbest=0
	Ebest=10000
}{
	if(zeta!=$3){
		if (zeta != 12345) printf("%5d %5d %12.8f %12.8f %12.8f\n", phi1,phi2,zeta, Rbest,Ebest)
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
