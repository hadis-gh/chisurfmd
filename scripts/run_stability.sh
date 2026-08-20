# --areaL 17.142857142857 \

logFile="log.txt"
exec > >(tee "$logFile") 2>&1
echo "Logging to: $logFile"

# ---------------------- Temperatures ----------------------

temps=(0.000001 0.00001 0.0001 0.001 0.01)
# ---------------------- General run function ----------------------
# --particlesInit "/home/hadis/chisurfmd/notebooks/best_configs/OP/bcluster_8.bp" \

# temps=(0.02 0.10 0.20 0.30 0.45 0.60)

run_temperature() {
    T=$1

    /home/hadis/chisurfmd/build_chimorse_brendsen/examples/testMD \
        --particlesInit "/home/hadis/chisurfmd/notebooks/best_configs/OP/bcluster_8.bp"\
        --resetInitialVelocity \
        --temperature "$T" \
        --thermoInterval 0.05\
        --collisionFr 0.01\
        --areaL 17.142857142857\
        --fixRadius 100 \
        --time 100 \
        --dt 0.001 \
        --seed 15 \
        --chiMorseModel "/home/hadis/chimorse/examples/models/chimorse_all.json" \
        --saveFile "outputs_MC_T_${T}.bp"
}

export -f run_temperature

parallel -j 6 run_temperature ::: "${temps[@]}"