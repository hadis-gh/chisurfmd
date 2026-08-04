MD_EXE="/home/hadis/gl_chisurfmd/vector/build/examples/testMD"
outputDir="/home/hadis/gl_chisurfmd/vector/runs/chimorse_annealing"

potentialType="ChiMorse"
chiMorseModel="/home/hadis/chimorse/examples/models/chimorse_all.json"

particleNum=25
particlesType="OP"
seed=15

areaL=100
exclusionRadius=2
fixRadius=100
momentI=1

highTemperature=0.6
lowTemperature=0.1
stepTemperature=0.1

timeCooling=12
timeHeating=12
dt=0.00001

writeStateInterval=0.01
writeEnergyInterval=0.1
thermoInterval=0.05

collisionFr=0.0166
integration="VelocityVerlet"