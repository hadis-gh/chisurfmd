#!/usr/bin/env bash

particleNum=49
particlesInit="RANDOM"
areaL=20
particlesDensity=70.0
time=30
dt=0.005
thermoInterval=0.5
exclusionRadius=0.8
cutoff=10.0
seed=15
relaxationTime=40.0
collisionFr=$(echo "scale=4; 1 / 60" | bc -l)
lowTemperature=0.10
highTemperature=1.00
stepTemperature=0.05
saveParticles="N_Config"
neighborFile="N_particle_NeighborsMD.dat"

calNeighbors() {
    local startT="$1"
    local stepT="$2"
    local endT="$3"
    local fileName="$4"

    for tp in $(seq $(echo "$startT + $stepT" | bc) "$stepT" "$endT"); do
        formattedTp=$(printf "%.2f" "$tp")
        echo "$formattedTp" | tr '\n' '\t' >> "$fileName"

        prevT=$(echo "$tp - $stepT" | bc)
        formattedPrevT=$(printf "%.2f" "$prevT")

        testNParticleMD -T "$formattedTp" --particlesInit "${saveParticles}_${formattedPrevT}.dat" --saveParticles "${saveParticles}_${formattedTp}.dat" --appendLog --dt "$dt" -n "$particleNum" --seed "$seed" -t "$time" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"

        NearestNeighbor=$(tail -n 7 "$neighborFile" | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }')
        echo "$NearestNeighbor" >> "$fileName"

        echo "done for temperature: $formattedTp"
    done
}

> B_NeighborsCount.dat
> B_NeighborsCount2.dat

formattedStartT=$(printf "%.2f" "$highTemperature")
testNParticleMD -T "$formattedStartT" --particlesInit "$particlesInit" --saveParticles "${saveParticles}_${formattedStartT}.dat" --dt "$dt" -n "$particleNum" --seed "$seed" -t "$time" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"
echo "done for temperature: $formattedStartT"

calNeighbors "$highTemperature" "-$stepTemperature" "$lowTemperature" "B_NeighborsCount.dat"
calNeighbors "$lowTemperature" "$stepTemperature" "$highTemperature" "B_NeighborsCount2.dat"
