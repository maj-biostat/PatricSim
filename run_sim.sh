#!/sbin/bash

# /usr/bin/Rscript design_01.R -s 1 -n 1000 > log/sim1.log 2>&1 &
# /usr/bin/Rscript design_01.R -s 2 -n 1000 > log/sim2.log 2>&1 &

# s is scenario
# n is number of simulations
/usr/bin/Rscript design_02.R -s 1 -n 2000 > log/des02_scen1_sim.log 2>&1 &
/usr/bin/Rscript design_02.R -s 2 -n 2000 > log/des02_scen2_sim.log 2>&1 &
/usr/bin/Rscript design_02.R -s 3 -n 2000 > log/des02_scen3_sim.log 2>&1 &
/usr/bin/Rscript design_02.R -s 4 -n 2000 > log/des02_scen4_sim.log 2>&1 &
/usr/bin/Rscript design_02.R -s 5 -n 2000 > log/des02_scen5_sim.log 2>&1 &
/usr/bin/Rscript design_02.R -s 6 -n 2000 > log/des02_scen6_sim.log 2>&1 &

