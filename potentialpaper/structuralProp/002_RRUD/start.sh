

DIR=../../003_FILMS/002_RRUD

##############GET RONN
cp /dev/null data_all.dat
for file in  $DIR/best_geos/bcluster_*.dat
do 
	echo $file
 	python3 ../getRONN.py $file >> data_all.dat
done


###############GET RDF

#python3 ../getRDF.py -init 7.5 25.0 0.1 -out RDF_all.dat
#for i in {0..99}
#do 
#   ls best_geos/bcluster_$i.dat
#   python3 ../getRDF.py -read RDF_all.dat -file  best_geos/bcluster_$i.dat -out tmp
#   mv tmp RDF_all.dat
#done
