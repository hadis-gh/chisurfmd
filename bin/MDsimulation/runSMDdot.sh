#!/usr/bin/env bash

if [ -f config.sh ]; then
    . config.sh
fi

particleNum=${particleNum:-49}
areaL=${areaL:-20}
particlesDensity=${particlesDensity:-70.0}
particlesInit=${particlesInit:-"RANDOM"}
seed=${seed:-35}

timeSim=${timeSim:-100}
dt=${dt:-0.001}
writeStateInterval=${writeStateInterval:-0.05}
writeEnergyInterval=${writeEnergyInterval:-0.5}
thermoInterval=${thermoInterval:-100}

cutoff=${cutoff:-10.0}
exclusionRadius=${exclusionRadius:-0.8}
relaxationTime=${relaxationTime:-40.0}
collisionFr=$(echo "scale=4; 1 / 60" | bc -l)

integration=${integration:-"VelocityVerlet"}
targetT=${targetT:-1}

LJPhiOrder=${LJPhiOrder:-2}
momentI=${momentI:-1}
LJangularScale=${LJangularScale:-5}

saveParticles="N_Config"

MD_EXE=../test/testNParticleMD

CMD="$MD_EXE -T '$targetT' --particlesInit '$particlesInit' --saveParticles '${saveParticles}.dat' --dt '$dt' -n '$particleNum' \
    --seed '$seed' -t '$timeSim' --areaL '$areaL' --exclusionRadius '$exclusionRadius' \
    --thermoInterval '$thermoInterval'  --writeStateInterval '$writeStateInterval' --writeEnergyInterval '$writeEnergyInterval' \
    --integration '$integration' "
echo "$CMD"
eval $CMD

echo "done"
