

DIR=.

##############GET RONN
cp /dev/null data_all.dat
for file in  $DIR/best_geos/bcluster_*.dat
do 
	echo $file
 	python3 ../getRONN.py $file >> data_all.dat
done

