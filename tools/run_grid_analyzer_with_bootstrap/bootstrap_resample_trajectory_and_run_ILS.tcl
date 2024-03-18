# This script splits a trajectory in its constitutent frames in rst7 format.
#logfile console
#set traj {../../2217/T10/no_acqua_10_1_ALL_MD.nc}

# generate random integer number in the range [0,max)
proc RandomInt {max} {
    return [expr {int(rand()*$max)}]
}


proc resample {maxframes id} {
   for {set frame 0} {$frame <= $maxframes} {incr frame} {
      set sample [RandomInt $maxframes]
      mol addfile ./frames/traj_frame${sample}.rst7 step 1 first 0 last -1 waitfor all type rst7 $id
   }
}



# Implicit Ligand Sampling run script
# ===================================

# Running ILS calculation:

# vmd -dispdev text -e <name_of_this_file>

# You will need VMD 1.8.7. or higher.

# Change the input parameters below to your liking.
# The filenames used in this script are relative to the directory
# for which it was generated but you can of course change then.

# If you have a CUDA enables Nvidia GPU VMD will use the GPU for
# the computation. Since the all the GPU resources will then be
# used for ILS your graphical display will freeze up, so don't be
# surprised. After finishing each frame the display will briefly
# be updated and freeze again.

# Comment this line out to prevent the use of the CUDA implementation:
#set env(VMDCUDAILS) 1

# You might want to do a quick test with 1 frame first to see if
# the syntax is correct and to determine the approximate runtime
# per frame.

# Adjustable parameters:
# ----------------------

# First and last frames to process
set first 0
set last  9999

# Resolution of final downsampled map in Angstrom
set res  1.0

# Subsampling of each dimension during computation
# i.e. each gridpoint of the final map will actually
# be downsampled from subres^3 points.
set subres 3

# Control of the angular spacing of probe orientation vectors,
# i.e. the number of probe conformers generated.
#
#   1: use 1 orientation only
#   2: use 6 orientations (vertices of octahedron)
#   3: use 8 orientations (vertices of hexahedron)
#   4: use 12 orientations (faces of dodecahedron)
#   5: use 20 orientations (vertices of dodecahedron)
#   6: use 32 orientations (faces+vert. of dodecahedron)
#  >6: geodesic subdivisions of icosahedral faces
#      with frequency 1, 2, ...
#
#   For each orientation a number of rotamers will be
#   generated. The angular spacing of the rotations
#   around the orientation vectors is chosen to be about
#   the same as the angular spacing of the orientation
#   vector itself.
#   If the probe ha at least one symmetry axis then the
#   rotations around the orientation vectors are reduced
#   accordingly. If there is an infinite oder axis (linear
#   molecule) the rotation will be omitted.
#   In case there is an additional perpendicular C2 axis
#   the half of the orientations will be ignored so that
#   there are no antiparallel pairs.
#
#   Probes with tetrahedral symmetry:
#   Here conf denotes the number of rotamers for each of
#   the 8 orientations defined by the vertices of the
#   tetrahedron and its dual tetrahedron.
set orient   2

# Cutoff energy above which the occupancy is regarded zero
# For GPUs energies of more than 87 always correspond to
# floating point values of zero for the occupancy. Hence
# there is no point going higher than that.
set maxen  30

# Temperature of the MD simulation
set T   283

# Nonbonded interaction cutoff
set cutoff 10.0

# The minmax box defining the free energy map
# (two opposite corners of the grid)
#set minmax {{11.424 9.051 10.685} {46.937 48.134 51.943}}
#set sel [atomselect top "protein"]

# The DX file containing the free energy map
set dxfile_prefix resamp_2217

# -------------------------------------------------------

# Set up the probe

# WARNING: The probe parameters have only been verified to
#          reproduce experimental solvation free energies for
#          xenon, oxygen, nitric oxide and carbon monoxide probes.
#          Use parameters for other probes as a starting point
#          for optimization.
# PARAMETROS OBTENIDOS DEL CAMPO DE FUERZAS DE CHARMM!!!! (JUAN Gay)
set pmol [mol new ils_CO.xyz]
set psel [atomselect $pmol all]
set n [atomselect $pmol "name C"]
$n set radius     2.00;  # VDW radius
$n set occupancy -0.11;   # VDW epsilon
set o [atomselect $pmol "name O"]
$o set radius 1.70;  # VDW radius
$o set occupancy -0.12;   # VDW epsilon

# -------------------------------------------------------

# Load  the topology file in mol number MOLID
set molid [mol new ../../2217/T10/2217_con_eme.prmtop type parm7]

# Load the trajectory coordinate file (we use ASCII written with box) 
mol addfile ../../2217/T10/no_acqua_10_1_ALL_MD.nc type netcdf first 0 last $last step 1 waitfor all 

#Align trajectory to the alignment template, based on alpha carbon locations
#obtain number of frames en NUMFRAMES
set numframes [molinfo $molid get numframes]
puts "frames: $numframes"

# The minmax box defining the free energy map
# (two opposite corners of the grid)
set minmax [measure minmax [atomselect top "protein"] ]

# Delete original trajectory
mol delete $molid

# Load Amber Parameters (requieresi:
# amber_vdw.tcl, amber_vdw.parms and a previous top loaded)
source amber_vdw.tcl
read_parm p 

# Run the analysis

#set sel [atomselect $molid "all"]
#set sel [atomselect $molid "all"]

# Run ILS
set top {../../2217/T10/2217_con_eme.prmtop}
set rs_molid [mol new $top type parm7]
set maxresamples 1000
set maxfr [expr {$numframes -1}]
for {set resample 1} {$resample < $maxresamples} {incr resample} {
   resample $last $rs_molid
   set numframes [molinfo $rs_molid get numframes]
   puts "frames: $numframes"
   set dxfile "${dxfile_prefix}_${resample}.dx"
   set maxfr [expr {$numframes -1}]
   assign_vdw $rs_molid
   puts $dxfile
   volmap ils $rs_molid  $minmax -cutoff $cutoff  \
       -res $res -subres $subres -probesel $psel -orient $orient \
        -maxenergy $maxen \
       -T $T -first $first -last $last \
       -o $dxfile
   animate delete  beg 0 end $last skip 0 $rs_molid
}

# Quit VMD when done with ILS calculation
#quit

