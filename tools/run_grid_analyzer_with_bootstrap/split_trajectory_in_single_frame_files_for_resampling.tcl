# This script splits a trajectory in its constitutent frames in rst7 format.
#logfile console
set top {../../2217/T10/2217_con_eme.prmtop}
set traj {../../2217/T10/no_acqua_10_1_ALL_MD.nc}

mol load parm7 $top
mol addfile $traj step 1 first 0 last -1 waitfor all type netcdf

set numframes [molinfo 0 get numframes]
#for {set frame 0} {$frame < $numframes} {incr frame} {
#   animate write rst7 traj_frame${frame}.rst7 beg $frame end $frame skip 0 0
#}
   animate write rst7 ./frames/traj_frame0.rst7 beg 0 end 0 skip 0 0

quit

