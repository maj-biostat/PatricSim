#!/sbin/bash

/usr/bin/Rscript main.R -s 1 -n 1000 > log/sim1.log 2>&1 &
/usr/bin/Rscript main.R -s 2 -n 1000 > log/sim1.log 2>&1 &

