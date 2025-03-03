#!/usr/bin/env bash

if [ -f config.sh ]; then
    source config.sh
fi

numCoolingRuns=$(echo "($highTemperature - $lowTemperature) / $stepTemperature" | bc)
numHeatingRuns=$(echo "($highTemperature - $lowTemperature) / $stepTemperature" | bc)
totalRuns=$((numCoolingRuns + numHeatingRuns))

MD_EXE=${MD_EXE:-../test/testNParticleMD}
outputDir=${outputDir:-"./outputs"}
mkdir -p "$outputDir"

run_md() {
    local runIndex="$1"
    local particleInit="$2"
    local temperature="$3"
    local saveFile="${outputDir}/run_${runIndex}.bp"

    echo "Starting simulation #${runIndex} | Temperature: ${temperature}"

    ARGS=--particlesNum "$particleNum" \
        --seed "$seed" \
        --areaL "$areaL" \
        --exclusionRadius "$exclusionRadius" \
        --thermoInterval "$thermoInterval" \
        --writeStateInterval "$writeStateInterval" \
        --writeEnergyInterval "$writeEnergyInterval" \
        --integration "$integration" \
        --collisionFr "$collisionFr"

    $MD_EXE --runIndex "$runIndex" \
        --particlesInit "$particleInit" \
        --temperature "$temperature" \
        --time "$timeCooling" --dt "$dt" \
        --saveFile "$outputDir" \
        $ARGS

    if [ $? -ne 0 ]; then
        echo "Error: Simulation #${runIndex} failed."
        exit 1
    fi

    echo -e "\nSimulation ${runIndex} of ${totalRuns} completed successfully.\n"
    echo "------------------------------------------------------------"
}

run_temperature_loop() {
    local startTemp="$1"
    local stepTemp="$2"
    local numRuns="$3"
    local startIndex="$4"

    for ((i = 1; i <= numRuns; i++)); do
        runIndex=$((startIndex + i))
        temperature=$(echo "$startTemp + $i * $stepTemp" | bc)
        prevOutput="${outputDir}/run_$((runIndex - 1)).bp"

        run_md "$runIndex" "$prevOutput" "$temperature"
    done

}

initialRunIndex=0
run_md "$initialRunIndex" "RANDOM" "$highTemperature"

run_temperature_loop "$highTemperature" "$stepTemperature" "$numCoolingRuns" 0

run_temperature_loop "$lowTemperature" "$stepTemperature" "$numCoolingRuns" "$numCoolingRuns"

echo "All simulations completed. Total runs: ${totalRuns}"
