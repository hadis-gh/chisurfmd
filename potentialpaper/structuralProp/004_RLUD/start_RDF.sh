

DIR=.


###############GET RDF

python3 ../getRDF.py -init 7.5 25.0 0.1 -out RDF_all.dat
for file in  $DIR/best_geos/bcluster_*.dat
do 
	echo $file
   python3 ../getRDF.py -read RDF_all.dat -file  $file  -out tmp
   mv tmp RDF_all.dat
done
