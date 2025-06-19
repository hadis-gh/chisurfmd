#!/usr/bin/env bash

if [ -f config.sh ]; then
    source config.sh
fi

particleNum=${particleNum:-49}
areaL=${areaL:-20}
particlesDensity=${particlesDensity:-70.0}
particlesInit=${particlesInit:-"RANDOM"}
seed=${seed:-15}

timeCooling=${timeCooling:-100}
timeHeating=${timeHeating:-100}
dt=${dt:-0.001}
writeStateInterval=${writeStateInterval:-0.05}
writeEnergyInterval=${writeEnergyInterval:-0.5}
thermoInterval=${thermoInterval:-0.01}

cutoff=${cutoff:-10.0}
exclusionRadius=${exclusionRadius:-0.8}
relaxationTime=${relaxationTime:-40.0}
collisionFr=${collisionFr:-0.016}

integration=${integration:-"VelocityVerlet"}

LJPhiOrder=${LJPhiOrder:-"LJPhiOrder"}
momentI=${momentI:-"momentI"}
LJangularScale=${LJangularScale:-"LJangularScale"}

lowTemperature=${lowTemperature:-0.02}
highTemperature=${highTemperature:-1.00}
stepTemperature=${stepTemperature:-0.01}

numCoolingRuns=$(echo "($highTemperature - $lowTemperature) / $stepTemperature" | bc)
numHeatingRuns=$(echo "($highTemperature - $lowTemperature) / $stepTemperature" | bc)
totalRuns=$((numCoolingRuns + numHeatingRuns))

MD_EXE=${MD_EXE:-../test/testNParticleMD}
outputDir=${outputDir:-"./outputs"}
mkdir -p "$outputDir"



echo "Starting simulation #${runIndex} | Temperature: ${temperature}"

$MD_EXE --runIndex 0 \
    --particlesInit "input.pb" \
    --temperature 0.02 \
    --time "$timeCooling" --dt "$dt" \
    --particlesNum "$particleNum" \
    --seed "$seed" \
    --areaL "$areaL" \
    --exclusionRadius "$exclusionRadius" \
    --thermoInterval "$thermoInterval" \
    --writeStateInterval "$writeStateInterval" \
    --writeEnergyInterval "$writeEnergyInterval" \
    --collisionFr "$collisionFr" \
    --integration "$integration" \
    --LJPhiOrder "$LJPhiOrder" \
    --momentI "$momentI" \
    --LJangularScale "$LJangularScale" \
    --saveFile "${outputDir}/run_${runIndex}.bp"



