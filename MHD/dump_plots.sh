#!/bin/bash
triplot=../RUN/PLOTTING/triplot.py
mesh=../PRELOAD/inputformaxwell.msh
for file in plot_*.dat; do
    python $triplot $mesh $file 2 ${file%%.*}_ReR.pdf &
    python $triplot $mesh $file 3 ${file%%.*}_ImR.pdf &
    python $triplot $mesh $file 4 ${file%%.*}_ReZ.pdf &
    python $triplot $mesh $file 5 ${file%%.*}_ImZ.pdf &
    python $triplot $mesh $file 6 ${file%%.*}_Rephi.pdf &
    python $triplot $mesh $file 7 ${file%%.*}_Imphi.pdf &
done

for file in presn*.dat; do
    python $triplot $mesh $file 0 ${file%%.*}_Re.pdf &
    python $triplot $mesh $file 1 ${file%%.*}_Im.pdf &
done

python $triplot $mesh fluxvar.dat 0 psi.pdf &
python $triplot $mesh fluxvar.dat 1 q.pdf &
python $triplot $mesh fluxvar.dat 2 dens.pdf &
python $triplot $mesh fluxvar.dat 3 temp.pdf &
python $triplot $mesh fluxvar.dat 4 pres0.pdf &

wait
