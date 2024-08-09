#!/usr/bin/env bash

particleNum=20
particleInit="DLA"
areaL=20
time=50
dt=0.01
thermoInterval=0.5
exclusionRadius=0.8
cutoff=10.0
seed=15
relaxationTime=40.0
temperature=0.1
collisionFr=$(( 1 / 60 ))
startCutoff=10
finalCutoff=50
stepCutoff=2
runTime="N_runnnigTime.dat"

> B_cutoffList.dat
> B_runTimeList.dat

make testNParticleMD ||  { echo "compiling failed!"; exit 1; }

for cf in $(seq "$startCutoff" "$stepCutoff" "$finalCutoff")
do
    echo "$cf" >> B_cutoffList.dat
    
    ./testNParticleMD --seed "$seed" -t "$time" --particleInit "$particleInit" -T "$temperature" -n "$particleNum" --areaL "$areaL" --cutoff "$cutoff" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval" ||  { echo "run failed!"; exit 1; }

    time=$(tail -n 1 $runTime)
    echo "$time" >> B_runTimeList.dat

    echo "done for cutoff: $cf"
done