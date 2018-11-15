#!/bin/sh

cat ndata_final.out | awk '{if(NR%5==3) print}' | cat -n | awk '{print $1,$8}' > tmp
awk -f cv.awk tmp | awk '{print 100, 0.10/(0.001*($2-$1))}'

