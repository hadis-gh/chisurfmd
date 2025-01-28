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
collisionFr=${collisionFr:-0.016}

integration=${integration:-"VelocityVerlet"}

lowTemperature=${lowTemperature:-0.02}
highTemperature=${highTemperature:-1.00}
stepTemperature=${stepTemperature:-0.01}

maxRunIndex=($highTemperature - $lowTemperature) / $stepTemperature

MD_EXE=../test/testNParticleMD

run_md () {
    local runIndex="$1"
    local particleInit="$2"
    local temperature="$3"

    $MD_EXE --runIned "$runIndex" --particlesInit "$particlesInit" -T "$temperature" \
            --dt "$dt" -n "$particleNum" --seed "$seed" -t "$timeCooling" --areaL "$areaL" --exclusionRadius "$exclusionRadius" \
            --thermoInterval "$thermoInterval"  --writeStateInterval "$writeStateInterval" --writeEnergyInterval "$writeEnergyInterval" --collisionFr "$collisionFr" \
            --integration "$integration" --LJPhiOrder "$LJPhiOrder" --momentI "$momentI" --LJangularScale "$LJangularScale"
}

run_md "0" "RANDOM" "$highTemperature"

for ((i=1; i<=$maxRunIndex/2; i++)); do
    run_md "$i" "simulation_output_$((i-1)).dat" "$((highTemperature - i * stepTemperature))"
done

for ((i=$maxRunIndex/2+1; i<=$maxRunIndex; i++)); do
    run_md "$i" "simulation_output_$((i-1)).dat" "$((lowTemperature + (i - maxRunIndex/2) * stepTemperature))"
done