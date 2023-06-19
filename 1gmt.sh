#!/bin/bash
### plot is1 canopy height model 
### 06/06/2023

##Interpolating to DEM, nearneighbor
#
# The two shorthands -Rg and -Rd stand for global domain (0/360 and -180/+180 in longitude respectively, 
# with -90/+90 in latitude)
#
# distance units 
# e Meter; k Km; d degree of arc; m minute of arc; 
#
#gmt nearneighbor is1_global.csv -Rd -I0.125  -Gis1.grd -S1d -N8/6
gmt makecpt -Cviridis -T0/40/1 -Z > elevation.cpt
#gmt grdgradient is1.grd -Ne0.8 -A100 -fg -Ggradient.nc
# -R-180/180/-85/85 
gmt grdimage is1.grd  -Igradient.nc -JN6i -R-180/180/-70/80  -P -Ba -B+t"IS1 canopy height (0.125 \260)" -Celevation.cpt -K -Q > fig.ps
gmt psscale -Dx8c/-1c+w8c/0.3c+jTC+h -Bxaf -By+lm  -Celevation.cpt  -O>> fig.ps
gmt psconvert fig.ps -A -Tj