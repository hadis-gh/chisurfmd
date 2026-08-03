#!/usr/bin/env bash

set -eu

config=${1:-"config.sh"}
if [ ! -f "$config" ]; then
    echo "Error: '$config' not found."
    exit 1
fi
source "$config"

MD_EXE=${MD_EXE:-"../test/testNParticleMD"}
outputDir=${outputDir:-"./outputs"}
mkdir -p "$outputDir"

logFile="log.txt"
exec > >(tee "$logFile") 2>&1
echo "Logging to: $logFile"

simulationType=${1:-"oriented"}
particleInit=${2:-"RANDOM"}
temperature=${3:-"$highTemperature"}
saveFile=${4:-"${outputDir}/final_run.bp"}

echo "============================================="
echo "    Single Molecular Dynamics Simulation     "
echo "============================================="

args=(
    --particlesInit "$particleInit"
    --particlesType "$particlesType"
    --temperature "$temperature"
    --time "$timeCooling" --dt "$dt"
    --saveFile "$saveFile"
    --particlesNum "$particleNum"
    --seed "$seed"
    --areaL "$areaL"
    --exclusionRadius "$exclusionRadius"
    --thermoInterval "$thermoInterval"
    --writeStateInterval "$writeStateInterval"
    --writeEnergyInterval "$writeEnergyInterval"
    --integration "$integration"
    --collisionFr "$collisionFr"
)

if [ "$simulationType" == "rotation" ]; then
    args+=(--LJangularScale "$LJangularScale"
           --LJPhiOrder "$LJPhiOrder"
           --LJalpha "$LJalpha"
           --momentI "$momentI")
fi

"$MD_EXE" "${args[@]}" || { 
    printf "Error: Simulation failed.\n";
    exit 1; 
}

printf "\nSimulation completed successfully.\n";
echo "------------------------------------------------------------"
