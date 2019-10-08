#!/sbin/bash

# /usr/bin/Rscript design_01.R -s 1 -n 1000 > log/sim1.log 2>&1 &
# /usr/bin/Rscript design_01.R -s 2 -n 1000 > log/sim2.log 2>&1 &

# s is scenario
# n is number of simulations
# /usr/bin/Rscript design_04.R -s 5 -n 10 -m 10000 
/usr/bin/Rscript design_04.R -s 1  -n 1000 -m 10000 > log/des04_scen01_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 2  -n 1000 -m 10000 > log/des04_scen02_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 3  -n 1000 -m 10000 > log/des04_scen03_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 4  -n 1000 -m 10000 > log/des04_scen04_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 5  -n 1000 -m 10000 > log/des04_scen05_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 6  -n 1000 -m 10000 > log/des04_scen06_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 7  -n 1000 -m 10000 > log/des04_scen07_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 8  -n 1000 -m 10000 > log/des04_scen08_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 9  -n 1000 -m 10000 > log/des04_scen09_sim.log 2>&1 &
/usr/bin/Rscript design_04.R -s 10 -n 1000 -m 10000 > log/des04_scen10_sim.log 2>&1 &
