#!/usr/bin/env bash

start_time=$(date +%s)

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

LJalpha=${LJalpha:-"LJalpha"}
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

logFile="log.txt"
exec > >(tee "$logFile") 2>&1
echo "Logging to: $logFile"

run_md() {
    local runIndex="$1"
    local particleInit="$2"
    local temperature="$3"

    echo "Starting simulation #${runIndex} | Temperature: ${temperature}"

    $MD_EXE --particlesInit "$particleInit" \
        --temperature "$temperature" \
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
        --LJalpha "$LJalpha" \
        --LJPhiOrder "$LJPhiOrder" \
        --momentI "$momentI" \
        --LJangularScale "$LJangularScale" \
        --saveFile "${outputDir}/run_${runIndex}.bp"

    if [ $? -ne 0 ]; then
        echo "Error: Simulation #${runIndex} failed."
        exit 1
    fi

    echo -e "\nSimulation ${runIndex} of ${totalRuns} completed successfully.\n"
    echo "------------------------------------------------------------"
}

# Initial random configuration
initialRunIndex=0
run_md "$initialRunIndex" "RANDOM" "$highTemperature"

# Cooling loop
for ((i = 1; i <= numCoolingRuns; i++)); do
    temperature=$(echo "$highTemperature - $i * $stepTemperature" | bc)
    prevOutput="${outputDir}/run_$((i - 1)).bp"

    run_md "$i" "$prevOutput" "$temperature"
done

# Heating loop
for ((i = 1; i <= numHeatingRuns; i++)); do
    runIndex=$((numCoolingRuns + i))
    temperature=$(echo "$lowTemperature + $i * $stepTemperature" | bc)
    prevOutput="${outputDir}/run_$((runIndex - 1)).bp"
    run_md "$runIndex" "$prevOutput" "$temperature"
done

echo "All simulations completed. Total runs: ${totalRuns}"


# Record time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

printf "Total execution time: %02d:%02d:%02d (hh:mm:ss)\n" "$hours" "$minutes" "$seconds"