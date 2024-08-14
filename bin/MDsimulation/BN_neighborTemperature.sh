#!/usr/bin/env bash

particleNum=49
particlesInit="RANDOM"
areaL=20
particlesDensity=70.0
time=100
dt=0.001
thermoInterval=0.1
exclusionRadius=0.8
cutoff=10.0
seed=15
relaxationTime=40.0
collisionFr=$(echo "scale=4; 1 / 60" | bc -l)
lowTemperature=0.02
highTemperature=1.00
stepTemperature=0.02
saveParticles="N_Config"
neighborFile="N_particle_NeighborsMD.dat"

calNeighbors() {
    local startT="$1"
    local stepT="$2"
    local endT="$3"
    local fileName="$4"
    prevT=$startT
    for tp in $(seq $(echo "$startT + $stepT" | bc) "$stepT" "$endT"); do
        echo "$tp" | tr '\n' '\t' >> "$fileName"

        testNParticleMD -T "$tp" --particlesInit "${saveParticles}_${prevT}.dat" --saveParticles "${saveParticles}_${tp}.dat" --appendLog --dt "$dt" -n "$particleNum" --seed "$seed" -t "$time" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"

        NearestNeighbor=$(tail -n 7 "$neighborFile" | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }')
        echo "$NearestNeighbor" >> "$fileName"

        echo "done for temperature: $tp"
        prevT=$tp
    done
}

> B_NeighborsCount.dat
> B_NeighborsCount2.dat

testNParticleMD -T "$highTemperature" --particlesInit "$particlesInit" --saveParticles "${saveParticles}_${highTemperature}.dat" --dt "$dt" -n "$particleNum" --seed "$seed" -t "$time" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"
echo "done for temperature: $highTemperature"

calNeighbors "$highTemperature" "-$stepTemperature" "$lowTemperature" "B_NeighborsCount.dat"
calNeighbors "$lowTemperature" "$stepTemperature" "$highTemperature" "B_NeighborsCount2.dat"
