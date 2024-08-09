#!/usr/bin/env bash

particleNum=100
particleInit="RANDOM"
areaL=40
time=20
dt=0.01
thermoInterval=0.5
exclusionRadius=0.8
cutoff=10.0
seed=15
relaxationTime=40.0
collisionFr=$(( 1 / 60 ))
lowTemperature=0.1
highTemperature=10
stepTemperature=0.2
saveParticles="N_Config"
neighborFile="N_particle_NeighborsMD.dat"

calNeighbors() {
    local startT="$1"
    local stepT="$2"
    local endT="$3"
    local fileName="$4"

    for tp in $(seq "$((startT + stepT))" "$stepT" "$endT"); do
        echo "$tp" | tr '\n' '\t' >> "$fileName"

        prevT=$((tp - stepT))

        ./testNParticleMD --seed "$seed" -t "$time" --particleInit "${saveParticles}_${prevT}.dat" -T "$tp" --saveParticles "${saveParticles}_${tp}.dat" -n "$particleNum" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"

        NearestNeighbor=$(tail -n 7 "$neighborFile" | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }')
        echo "$NearestNeighbor" >> "$fileName"

        echo "done for temperature: $tp"
    done
}

> B_NeighborsCount.dat
> B_NeighborsCount2.dat

./testNParticleMD --seed "$seed" -t "$time" --particleInit "$particleInit" -T "$highTemperature" --saveParticles "${saveParticles}_${highTemperature}.dat" -n "$particleNum" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"

calNeighbors "$highTemperature" "-$stepTemperature" "$lowTemperature" "B_NeighborsCount.dat"

calNeighbors "$lowTemperature" "$stepTemperature" "$highTemperature" "B_NeighborsCount2.dat"
