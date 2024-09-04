#!/usr/bin/env bash

if [ -f config.sh ]; then
    . config.sh
fi

particleNum=49
particlesInit="RANDOM"
areaL=20
particlesDensity=70.0
timeCooling=${timeCooling:-100}
timeHeating=${timeHeating:-50}
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

MD_EXE=../test/testNParticleMD

calNeighbors() {
    local startT="$1"
    local stepT="$2"
    local endT="$3"
    local fileName="$4"
    local time="$5"
    prevT=$startT
    for tp in $(seq $(echo "$startT + $stepT" | bc) "$stepT" "$endT"); do
        echo -n -e "${tp}\t" >> "$fileName"

        $MD_EXE -T "$tp" --particlesInit "${saveParticles}_${prevT}.dat" --saveParticles "${saveParticles}_${tp}.dat" --appendLog --dt "$dt" -n "$particleNum" --seed "$seed" -t "$time" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"

        NearestNeighbors=$(tail -n 7 "$neighborFile" | awk '
            {
                for (i = 2; i <= NF; i++) {
                    sum[i] += $i;
                }
            }
            END {
                for (i = 2; i <= length(sum) + 1; i++) {
                    avg = sum[i] / NR;
                    printf "%.6f\t", avg;
                }
            }
        ')

        echo -e "$NearestNeighbors" >> "$fileName"
        echo "done for temperature: $tp"
        prevT=$tp
    done
}

> B_NeighborsCount.dat
> B_NeighborsCount2.dat

$MD_EXE -T "$highTemperature" --particlesInit "$particlesInit" --saveParticles "${saveParticles}_${highTemperature}.dat" --dt "$dt" -n "$particleNum" --seed "$seed" -t "$timeCooling" --areaL "$areaL" --exclusionRadius "$exclusionRadius" --thermoInterval "$thermoInterval"
echo "done for temperature: $highTemperature"

echo "Start Cooling Process"
calNeighbors "$highTemperature" "-$stepTemperature" "$lowTemperature" "B_NeighborsCount.dat" "$timeCooling"
echo "Start Heating Process"
calNeighbors "$lowTemperature" "$stepTemperature" "$highTemperature" "B_NeighborsCount2.dat" "$timeHeating"
