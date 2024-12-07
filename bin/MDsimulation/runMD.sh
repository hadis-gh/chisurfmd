#!/usr/bin/env bash

if [ -f config.sh ]; then
    . config.sh
fi

particleNum=${particleNum:-49}
areaL=${areaL:-20}
particlesDensity=${particlesDensity:-70.0}
particlesInit=${particlesInit:-"RANDOM"}
seed=${seed:-15}

timeCooling=${timeCooling:-100}
timeHeating=${timeHeating:-50}
dt=${dt:-0.001}
writeStateInterval=${writeStateInterval:-0.05}
writeEnergyInterval=${writeEnergyInterval:-0.5}
thermoInterval=${thermoInterval:-0.01}

cutoff=${cutoff:-10.0}
exclusionRadius=${exclusionRadius:-0.8}
relaxationTime=${relaxationTime:-40.0}
collisionFr=$(echo "scale=4; 1 / 60" | bc -l)

integration=${integration:-"VelocityVerlet"}

lowTemperature=${lowTemperature:-0.02}
highTemperature=${highTemperature:-1.00}
stepTemperature=${stepTemperature:-0.01}

saveParticles="N_Config"
coolingDir="cooling_configs"
heatingDir="heating_configs"
neighborFile="N_particle_NeighborsMD.dat"

LJPhiOrder=${LJPhiOrder:-2}
momentI=${momentI:-1}

MD_EXE=../test/testNParticleMD

mkdir -p "$coolingDir"
mkdir -p "$heatingDir"

calNeighbors() {
    local startT="$1"
    local stepT="$2"
    local endT="$3"
    local fileName="$4"
    local time="$5"
    local configDir="$6"
    prevT=$startT
    for tp in $(seq $(echo "$startT + $stepT" | bc) "$stepT" "$endT"); do
        echo -n -e "${tp}\t" >> "$fileName"

        $MD_EXE -T "$tp" --particlesInit "${configDir}/${saveParticles}_${prevT}.dat" \
            --saveParticles "${configDir}/${saveParticles}_${tp}.dat" --appendLog --dt "$dt" \
            -n "$particleNum" --seed "$seed" -t "$time" --areaL "$areaL" --exclusionRadius "$exclusionRadius" \
            --thermoInterval "$thermoInterval" --writeStateInterval "$writeStateInterval" \
            --writeEnergyInterval "$writeEnergyInterval" --integration "$integration" \
            --LJPhiOrder "$LJPhiOrder" --momentI "$momentI" 

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

$MD_EXE -T "$highTemperature" --particlesInit "$particlesInit" --saveParticles "${coolingDir}/${saveParticles}_${highTemperature}.dat" \
    --dt "$dt" -n "$particleNum" --seed "$seed" -t "$timeCooling" --areaL "$areaL" --exclusionRadius "$exclusionRadius" \
    --thermoInterval "$thermoInterval"  --writeStateInterval "$writeStateInterval" --writeEnergyInterval "$writeEnergyInterval"  \
    --integration "$integration" --LJPhiOrder "$LJPhiOrder" --momentI "$momentI" 

echo "done for temperature: $highTemperature"

echo "Start Cooling Process"
calNeighbors "$highTemperature" "-$stepTemperature" "$lowTemperature" "B_NeighborsCount.dat" "$timeCooling" "$coolingDir"  --integration "$integration"

lastCoolingConfig="${coolingDir}/${saveParticles}_${lowTemperature}.dat"
echo "Start Heating Process from last cooling configuration: $lastCoolingConfig"

$MD_EXE -T "$lowTemperature" --particlesInit "$lastCoolingConfig" --saveParticles "${heatingDir}/${saveParticles}_${lowTemperature}.dat" \
    --appendLog --dt "$dt" -n "$particleNum" --seed "$seed" -t "$timeHeating" --areaL "$areaL" --exclusionRadius "$exclusionRadius" \
    --thermoInterval "$thermoInterval"  --writeStateInterval "$writeStateInterval" --writeEnergyInterval "$writeEnergyInterval"  \
    --integration "$integration"--LJPhiOrder "$LJPhiOrder" --momentI "$momentI" 

calNeighbors "$lowTemperature" "$stepTemperature" "$highTemperature" "B_NeighborsCount2.dat" "$timeHeating" "$heatingDir"
