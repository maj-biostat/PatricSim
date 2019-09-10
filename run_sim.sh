#!/sbin/bash

/usr/bin/Rscript design_01.R -s 1 -n 1000 > log/sim1.log 2>&1 &
/usr/bin/Rscript design_01.R -s 2 -n 1000 > log/sim2.log 2>&1 &

