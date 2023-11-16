source /home/diegoa/inv/parma/daniele/analisis/ilsanalyzer.tcl

###################################################################
#           VARIABLE DEFINITION
###################################################################

# PDB file
set pdbfile_ini "2217_10C_frame0.pdb"

# DX file
set dxfile_ini "analisis/2217_10C_CO.dx"

# Define reference coordinates
set ref_coord_list [list {29.450000 35.230998 30.190000} {29.995001 24.528000 30.538000}]

###################################################################
#           MAIN PROGRAM
###################################################################

# Load pdb file with protein structure to VMD.
load_pdb_to_vmd $pdbfile_ini

# Load dx file with ILS grid results to VMD.
load_dx_to_vmd $dxfile_ini

# Read DX file grid to analyze.
read_DX $dxfile_ini

# Global Optimization
#set radius 7
#set outfile "global_opt.pdb"
#global_search_min $ref_coord_list $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta $radius $outfile

# Nudged Elastic Band
set start_indx_list [list  \
   {4 25} {25 21} {25 29} {25 35} {29 31} {29 55} {55 56} \
   {29 28} {28 39} {21 20} {20 36} {36 45} {45 77} {36 67} \
   {36 69} {69 61} {61 62} {67 69} {69 74} {74 73} {21 27} \
   {21 19} {27 26} {19 26} {26 37} {37 38} {37  7} { 4 13} \
   { 4  5} { 5  6} {6  11} { 4 10} {10 14} {14 11} \
]

mol new "global_opt.pdb" type {pdb} first 0 last -1 step 1 waitfor 1
mol modstyle 0 2 VDW 0.500000 20.000000
mol modcolor 0 2 Occupancy
mol modstyle 0 1 Isosurface 8.0 0 0 0 1 1

set out_fh1 [open "borders.pdb" w]
set out_fh2 [open "initial_path.pdb" w]
set out_fh3 [open "optimized_path.pdb" w]
set cnt 1
set ntot 0
set connect_list [list ]
global connect_list
foreach a $start_indx_list {
   puts "STARTING NEB ANALYSIS No $cnt/[llength $start_indx_list]"
   set ini [lindex $a 0]
   set end [lindex $a 1]
   set sel_ini [ atomselect top "index $ini" ]
   set sel_end [ atomselect top "index $end" ]
   set pos1 [ $sel_ini get {x y z} ]
   set pos2 [ $sel_end get {x y z} ]
   set pos1x [lindex $pos1 0 0]
   set pos1y [lindex $pos1 0 1]
   set pos1z [lindex $pos1 0 2]
   set pos2x [lindex $pos2 0 0]
   set pos2y [lindex $pos2 0 1]
   set pos2z [lindex $pos2 0 2]
   set npnts [grid_NEB $pos1x $pos1y $pos1z \
                       $pos2x $pos2y $pos2z \
                       $xOrigen $yOrigen $zOrigen \
                       $xdelta $ydelta $zdelta \
                       $out_fh3 $out_fh2 $out_fh1 \
                       $ntot ]
   set ntot [expr { $ntot + $npnts } ]
   set cnt [expr $cnt + 1]
}
close $out_fh3
close $out_fh2
close $out_fh1




#set start_points [list {23.995001 35.528000 31.538000 26.995001 34.528000 32.537998} {26.995001 34.528000 32.537998 28.995001 36.528000 36.537998} {28.995001 36.528000 36.537998 29.995001 32.528000 38.537998} {28.995001 36.528000 36.537998 31.995001 40.528000 36.537998} {26.995001 34.528000 32.537998 26.995001 29.528000 32.537998} ]
#set out_fh2 [open "initial_path.pdb" w]
#set out_fh3 [open "optimized_path.pdb" w]
#for {set i 0} { $i < [ llength $start_points ] } { incr i } {
#   puts "STARTING NEB ANALYSIS No $i"
#   set vector [lindex $start_points $i]
#   set pos1x [lindex $vector 0]
#   set pos1y [lindex $vector 1]
#   set pos1z [lindex $vector 2]
#   set pos2x [lindex $vector 3]
#   set pos2y [lindex $vector 4]
#   set pos2z [lindex $vector 5]
#   grid_NEB $pos1x $pos1y $pos1z \
#            $pos2x $pos2y $pos2z \
#            $xOrigen $yOrigen $zOrigen \
#            $xdelta $ydelta $zdelta \
#            $out_fh3 $out_fh2
#}
#close $out_fh3
#close $out_fh2



# Optimize from initial position.
#set opt [optimizer $ref_coord_x $ref_coord_y $ref_coord_z $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta]

# Now I have to separate the previous output into different variables because TCL is 
# that full of shit.
#set refx [lindex $opt 0]
#set refy [lindex $opt 1]
#set refz [lindex $opt 2]
#set Gmin [lindex $opt 3]

# Print PDB file with optimized geometries.
#set g 1
#set leaveopen 0
#set outfile "min.pdb"
#print_pdb $outfile $xOrigen $yOrigen $zOrigen $refx $refy $refz $xdelta $ydelta $zdelta $Gmin $g $leaveopen

#mol new $outfile type {pdb} first 0 last -1 step 1 waitfor 1
#mol new "global_opt.pdb" type {pdb} first 0 last -1 step 1 waitfor 1
#mol modstyle 0 2 VDW 0.500000 20.000000
#mol modcolor 0 2 Occupancy
#mol modstyle 0 1 Isosurface 8.0 0 0 0 1 1

mol new "initial_path.pdb" type {pdb} first 0 last -1 step 1 waitfor 1
mol modstyle 0 3 VDW 0.100000 20.000000
mol modcolor 0 3 ColorID 7

mol new "borders.pdb" type {pdb} first 0 last -1 step 1 waitfor 1
mol modstyle 0 4 VDW 0.100000 20.000000
mol modcolor 0 4 ColorID 3

mol new "optimized_path.pdb" type {pdb} first 0 last -1 step 1 waitfor 1
mol modstyle 0 5 VDW 0.100000 20.000000
mol modcolor 0 5 ColorID 4

set sel [atomselect 5 "all"]
$sel setbonds $connect_list 
