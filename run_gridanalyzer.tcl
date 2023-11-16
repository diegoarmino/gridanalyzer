source /home/diegoa/dev/gridanalyzer/gridanalyzer_lib.tcl
###################################################################
# GRID ANALYZER EXECUTION SCRIPT EXAMPLE
###################################################################
# This script is intended as an example for gridanalyzer library
# excecution. It is implemented to be used in three modes, 1, 2 or
# 3 by setting the "mode" variable, as follows.
#
# mode=0 ---> Gaussian kernel representation calculation (generates a .dx file).
# mode=1 ---> global minima search.
# mode=2 ---> Pathway search
#
# GAUSSIAN KERNEL REPRESENTATION (GKR):
# GKR is a method similar in purpuse to a histogram, so it can be used to obtain
# probabilities distribution functions (in this case in three dimensions), and
# so, by using the known formula G(x,y,z)=-RT ln( p(x,y,z) ) it is also possible to obtain
# the equivalent free energy. The advantage of this method compared to a simple
# histogram is that it can be shown that it converges much faster than the former
# (i.e. produces less noisy results), and doesnt depend on a given binning. A 2D histogram
# can be thought of as constructing a frequency plot by asigning a square to each
# event and just stacking squares in the corresponding bin. In GKR we use gaussian functions
# instead of squares and use no bins. So each gaussian function is evaluated in its right
# place, and are sumed toghether to get the final frequency plot (and normalized to get
# the probability distribution. There is one adjustable parameter though, the gaussian
# function width, sometimes called sigma, in reference to the normal distribution
# function tradicional nomenclature. There are ways to judiciously select the best sigma
# for each data set. But we have not yet implemented any of those methods.
# The script to runing GKR has two parts. Part 1 consist of a preparative step, in
# which a selection is made (for example all oxygen atoms of water molecules near the protein),
# and an intermediate file containing all the xyz coordinates of such water molecules for all
# frames of the trajectory is produced, alongside the coordinates of the corners of the
# smallest box containing all such selected water molecules.
# Part 2: The intermediate file is fed to the executable for the actual GKR analysis, which
# was written in fortran for performance reasons. This program calculates the GKR probability
# distribution function estimation as well as the equivalent free energy function, all evaluated
# onto a 3D grid and stored as a file in .DX format. The corresponding files are named
#
# free_energy.dx
# probability.dx
#
# Their names should be self-explanatory. The free energy file can then be procecessed using
# gridanalyze in order to obtain the minima and the minimum energy pathways between them (see below).
#
# GLOBAL MINIMA SEARCH
# This functionality gets as input a .dx file and searches for its minima. In this way it is possible
# to obtain, for example, CO docking sites in the case the dx file is the output of ILS analysis, or
# water sites, in the case the dx file is a free energy function obtained from GKR analysis of a
# trajectory. The program produces a pdb file contining all minima found with their corresponging
# energies in the "occupancy" slot of the pdb file, so they can be uploaded into VMD and colored
# using "occupancy" colour method in order to get a graphical representation of the positions and
# energies of the minima.
#
# MIGRATION PATHWAYS SEARCH
# This functionality searches for minimum energy pathways between the aforementioned minima using an
# algorithm similar to Nudged Elastic Band. This is, an initial pathway is drawn as an interpolation
# between the positions of two minima, and then they are subject to a restrained optimization, in
# such a way that each intermediate point can only move on a plane (or more accurately, a thin slab)
# perpendicular to the line passing through both minima.
# The program produces a pdb file containing all pathway points found in pdb format.
#
# IN SUMMARY
# Once the user has obtained the minima, he/she should select those the user whishes to probe for connecting
# pathways by setting the $start_indx_list variable (as a list of pairs of indexes). The indexes of each minima can
# be found by loading the optimization results pdb file into vmd, presing number 1, and then selecting the
# minima of interest. The index is reported at the console alongside other information such as coordinates, resname,
# etc. Once an index pairlist has been set up (see below, variable indx_list), the user may want to run
# the pathway finder, in order to do so,

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#            SET THE EXECUTION MODE HERE!
                     set mode 3

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###################################################################
#           VARIABLE DEFINITION
###################################################################

# PDB file
set pdbfile_ini "2217_10C_frame0.pdb"

# DX file
#set dxfile_ini "0030_10C_CO.dx"
set dxfile_ini "free_energy.dx"

# Define reference coordinates
set ref_coord_list [list {29.450000 35.230998 30.190000} {29.995001 24.528000 30.538000}]

# Set search radius. The search radius defines the volume to search for minima in, for the previoues
# reference coordinates. The optimizer won't consider grid points outside this cutoff.
set radius 12.0

# Define global optimization results file name (PDB). This file contains the coordinates and energy of
# each minima found.
set opt_filename "wat_gkr_global_opt.pdb"

# Define pathways filename (PDB). This file contains the pathways xyz coordinates as well as energies (except for
# minima points) for all the pathways found.
set pathways_filename "wat_gkr_optimized_path.pdb"

# Define the pruned .DX file name. This file is the same as the input .DX grid file, but with all values set to
# a very high energy value except for those close (up to a cutoff) to the migration pathways.
set pruned_grid_filename "wat_gkr_pruned_grid.dx"

# Set the index pair list for pathway search. It should be defined as a list of pairs with the following format:
#set start_indx_list [list  \
#    {4 0} \
#    {9 10} {10 5} {5 3} \
#]
set start_indx_list [list  \
    {22 10} {10 9} {9 14} {14 17} {14 13} {13 21} {21 66} {17 16} {16 24} {24 23} {23 20} {2 62} {2 19} \
    {58 57} {57 68} {68 73} {68 61} {61 62} {62 3} {3 19} {19 27} {61 69} {27 45} {27 20} {27 61}
]

###################################################################
#           SET GKR PARAMETERS AND INITIAL VARIABLES
###################################################################
#
# Selection for GKR analysis
set selection "noh and water and within 3 of protein or resname HEM"
# Path to GKR analysis executable
set path_GKR "/home/diegoa/dev/ilsanalyzer/gkranalysis.bin"
# Set the step of the grid for the sampling of the probability function in GKR analysis.
set gkr_delta 0.5
# Set the width of the gaussian kernel functions for GKR analysis
set gkr_sigma 0.5
# Set trajectory and prmtop filenames
set top {2217_con_eme.prmtop}
set traj {no_acqua_10_1_ALL_MD.nc}
# Set intermidiate file name.
set out_pdb "waters.pdb"
set max_frames 10000
# GKR analysis only stores a fraction of the trajectory on memory.
# Using the following variable set how big (in frames) this chunk
# of trajectory should be.
set batch 5000

# GKR analysis is performed by a code written in fortran. This code needs
# a file with a list of all the x y z coordinates of the probe (probably water)
# molecules for the whole trajectory, plus a heading for the coorinates of
# the bounding box containing all these molecules. The following variable defines
# if it is necessary to build this file. If the file has allready been produced
# then set do_build_input to 0, as the process of building this file is time-consuming.
#
# Build input file for gkranalyzer. 1=yes 0=no
#
set do_build_input 0

# Set the number of lines of input file for gkranalysis if known (produced in a previous run).
set gkr_line_cnt 3933787

###################################################################
#           SET LOAD RESULTS VARIABLES (MODE 3)
###################################################################
#
set ils_opt_filename "global_opt.pdb"
set ils_pathways_filename "optimized_path.pdb"
set ils_pruned_grid_filename "pruned_grid.dx"

set gkr_opt_filename "wat_gkr_global_opt.pdb"
set gkr_pathways_filename "wat_gkr_optimized_path.pdb"
set gkr_pruned_grid_filename "wat_gkr_pruned_grid.dx"

set ils_start_indx_list [list  \
   {25 4} { 4 13} {4 10} {10 14} {14 11} {11 6} {6 5} \
   {25 21} {21 36} {36 69} {69 74} {74 73} {69 61} {61 62} {36 45} {45 84} \
   {21 19} {19 27} {27 26} {26 37} {37 38} {37 7} \
   {25 29} {29 28} {29 39} {29 35} {35 31} {31 55} {55 56}
]

set ils_connect_list [list 1 {0 2} {1 3} 2 5 {4 6} {5 7} {6 8} {7 9} {8 10} {9 11} \
                           {10 12} 11 14 {13 15} {14 16} {15 17} 16 19 {18 20} \
                           19 22 {21 23} 22 25 {24 26} {25 27} {26 28} {27 29} \
                           28 31 {30 32} {31 33} {32 34} {33 35} 34 37 {36 38} \
                           {37 39} {38 40} {39 41} 40 43 {42 44} {43 45} {44 46} \
                           {45 47} 46 49 {48 50} {49 51} {50 52} {51 53} {52 54} \
                           53 56 {55 57} {56 58} 57 60 {59 61} 60 63 {62 64} {63 65} \
                           64 67 {66 68} {67 69} {68 70} 69 72 {71 73} {72 74} \
                           {73 75} {74 76} 75 78 {77 79} {78 80} {79 81} 80 83 \
                           {82 84} {83 85} {84 86} 85 88 {87 89} 88 91 {90 92} \
                           {91 93} 92 95 {94 96} {95 97} 96 99 {98 100} 99 102 \
                           {101 103} {102 104} {103 105} {104 106} 105 108 \
                           {107 109} {108 110} {109 111} 110 113 {112 114} \
                           {113 115} {114 116} 115 118 {117 119} 118 121 {120 122} \
                           121 124 {123 125} {124 126} {125 127} 126 129 {128 130} \
                           129 132 {131 133} {132 134} 133 ]


set gkr_start_indx_list [list  \
    {22 10} {10 9} {9 14} {14 17} {14 13} {13 21} {21 66} {17 16} {16 24} {24 23} {23 20} {2 62} {2 19} \
    {58 57} {57 68} {68 73} {68 61} {61 62} {62 3} {3 19} {19 27} {61 69} {27 45} {27 20} {27 61}
]

set gkr_connect_list [list 1 {0 2} {1 3} {2 4} 3 6 {5 7} {6 8} {7 9} 8 11 {10 12} {11 13} {12 14} {13 15} {14 16} {15 17} 16 19 {18 20} {19 21} {20 22} {21 23} {22 24} 23 26 {25 27} {26 28} {27 29} {28 30} {29 31} {30 32} {31 33} 32 35 {34 36} {35 37} 36 39 {38 40} {39 41} {40 42} {41 43} {42 44} {43 45} {44 46} {45 47} 46 49 {48 50} {49 51} {50 52} 51 54 {53 55} {54 56} {55 57} 56 59 {58 60} {59 61} {60 62} {61 63} {62 64} {63 65} 64 67 {66 68} {67 69} 68 71 {70 72} {71 73} {72 74} 73 76 {75 77} {76 78} {77 79} 78 81 {80 82} {81 83} {82 84} {83 85} {84 86} {85 87} 86 89 {88 90} {89 91} {90 92} {91 93} {92 94} 93 96 {95 97} {96 98} {97 99} {98 100} {99 101} {100 102} 101 104 {103 105} {104 106} {105 107} {106 108} {107 109} {108 110} {109 111} {110 112} 111 114 {113 115} {114 116} {115 117} {116 118} {117 119} {118 120} 119 122 {121 123} {122 124} {123 125} {124 126} {125 127} {126 128} 127 130 {129 131} {130 132} {131 133} {132 134} 133 136 {135 137} {136 138} {137 139} {138 140} {139 141} {140 142} {141 143} {142 144} {143 145} {144 146} {145 147} {146 148} {147 149} {148 150} {149 151} {150 152} 151 154 {153 155} {154 156} {155 157} {156 158} {157 159} {158 160} {159 161} {160 162} {161 163} 162 165 {164 166} {165 167} {166 168} {167 169} {168 170} {169 171} {170 172} {171 173} 172 175 {174 176} {175 177} {176 178} {177 179} {178 180} {179 181} {180 182} {181 183} {182 184} {183 185} {184 186} {185 187} {186 188} {187 189} 188 191 {190 192} {191 193} {192 194} {193 195} {194 196} {195 197} {196 198} {197 199} {198 200} 199]


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                   MAIN PROGRAM CODE STARTS HERE!
#     -----------------------------------------------------------------
#     DO NOT, I REPEAT, DO NOT EDIT THE FOLLOWING CODE UNLESS YOU KNOW
#                          WHAT YOU'RE DOING.
#     -----------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###################################################################
#  MODE 0: GAUSSIAN KERNEL REPRESENTATION ANALYSIS
###################################################################

if { $mode == 0 } {
   gkr_analysis_preparation $selection $batch $max_frames $top $traj $out_pdb $do_build_input
   exec $path_GKR "$out_pdb.xyz" $gkr_line_cnt $gkr_delta $gkr_sigma > "gkr.log"
}

###################################################################
#  MODE 1: DOCKING SITE SEACH (FREE ENERGY GRID OPTIMIZATION)
###################################################################

if { $mode == 1 } {

   # Load pdb file with protein structure to VMD.
   load_pdb_to_vmd $pdbfile_ini

   # Load dx file with ILS grid results to VMD.
   load_dx_to_vmd $dxfile_ini

   # Read DX file grid to analyze.
   read_DX $dxfile_ini

   set outfile $opt_filename
   global_search_min $ref_coord_list $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta $radius $outfile
}

###################################################################
#  MODE 2: PATHWAY SEARCH
###################################################################

if { $mode == 2 } {

   # Load pdb file with protein structure to VMD.
   load_pdb_to_vmd $pdbfile_ini

   # Load dx file with ILS grid results to VMD.
   load_dx_to_vmd $dxfile_ini

   # Read DX file grid to analyze.
   read_DX $dxfile_ini

   # Search for pathways.
   search_pathways $start_indx_list $opt_filename $pathways_filename

   # Created pruned .DX grid for better visualization of the isosurfaces related to the pathways found.
   write_pruned_DX $pruned_grid_filename $pathways_filename 4.0

   # Write connectivity list. This should be copied and pasted into the wat_connect_list
   # variable of mode 3 in order to be able to correctly visualize pathways results.
   puts $connect_list
}


###################################################################
#  MODE 3: LOAD RESULTS TO VMD
###################################################################


if { $mode == 3 } {
   # Load pdb file with protein structure to VMD.
   load_pdb_to_vmd $pdbfile_ini

   set ils_minlist [list ]
   foreach a $ils_start_indx_list {
      foreach b $a {
         if { [ lsearch $ils_minlist $b ] < 0 } {
            lappend ils_minlist $b }
      }
   }

   set gkr_minlist [list ]
   foreach a $gkr_start_indx_list {
      foreach b $a {
         if { [ lsearch $gkr_minlist $b ] < 0 } {
            lappend gkr_minlist $b }
      }
   }

   color Display Background white
   display cuedensity 0.050000
   display shadows on
   display ambientocclusion on

   mol modstyle    0 0 NewCartoon 0.300000 10.000000 4.100000 0
   mol modmaterial 0 0 AOChalky
   mol modcolor    0 0 ColorID 3
   mol modselect   1 0 resname HEM and noh
   mol modcolor    1 0 ColorID 6
   mol modmaterial 1 0 AOChalky
   mol modmaterial 2 0 AOChalky
   mol modselect   2 0 noh and resid 21 48 88 37 75

   mol new $ils_opt_filename type {pdb} first 0 last -1 step 1 waitfor 1
   mol modstyle 0 top VDW 0.500000 20.000000
   mol modcolor 0 top Occupancy
   mol modselect 0 top "index $ils_minlist"
   mol modmaterial 0 top AOChalky

   mol new $ils_pathways_filename type {pdb} first 0 last -1 step 1 waitfor 1
   mol modstyle 0 top Licorice 0.200000 12.000000 12.000000
   mol modcolor 0 top Occupancy
   set sel [atomselect top "all"]
   $sel setbonds $ils_connect_list
   mol modmaterial 0 top AOChalky

   mol new $ils_pruned_grid_filename type dx first 0 last -1 step 1 waitfor 1 volsets 0
   mol modstyle 0 top Isosurface 12.00 0 0 0 1 1
   mol modmaterial 0 top Transparent
   mol modcolor 0 top ColorID 7

   mol new $gkr_opt_filename type {pdb} first 0 last -1 step 1 waitfor 1
   mol modstyle 0 top VDW 0.500000 20.000000
   mol modcolor 0 top Occupancy
   mol modselect 0 top "index $gkr_minlist"
   mol modmaterial 0 top AOChalky

   mol new $gkr_pathways_filename type {pdb} first 0 last -1 step 1 waitfor 1
   mol modstyle 0 top Licorice 0.200000 12.000000 12.000000
   mol modcolor 0 top Occupancy
   set sel [atomselect top "all"]
   $sel setbonds $gkr_connect_list
   mol modmaterial 0 top AOChalky
   mol modselect 0 top "all not index 96 to 200 42 to 44"

   mol new $gkr_pruned_grid_filename type dx first 0 last -1 step 1 waitfor 1 volsets 0
   mol modstyle 0 top Isosurface 12.00 0 0 0 1 1
   mol modmaterial 0 top Transparent
   mol modcolor 0 top ColorID 7

   # EXTRA REPRESENTATIONS

   # ILS

   mol color ColorID 7
   mol representation VDW 0.500000 20.000000
   mol selection index 25 4 13 10 14 11 6 5
   mol material AOChalky
   mol addrep 1
   mol color ColorID 7
   mol representation Licorice 0.200000 12.000000 12.000000
   mol selection index 0 to 36
   mol material AOChalky
   mol addrep 2

   mol color ColorID 4
   mol representation VDW 0.500000 20.000000
   mol selection index 21 36 19 27 45 84 69 74 73 61 62 26 37 38 7
   mol material AOChalky
   mol addrep 1
   mol color ColorID 4
   mol representation Licorice 0.200000 12.000000 12.000000
   mol selection index 36 to 106
   mol material AOChalky
   mol addrep 2

   mol color ColorID 11
   mol representation VDW 0.500000 20.000000
   mol selection index 29 35 28 39 31 56
   mol material AOChalky
   mol addrep 1
   mol color ColorID 11
   mol representation Licorice 0.200000 12.000000 12.000000
   mol selection index 106 to 112 114 to 150
   mol material AOChalky
   mol addrep 2

   mol color ColorID 1
   mol representation VDW 0.700000 20.000000
   mol selection index 25
   mol material AOChalky
   mol addrep 1

   # GKR

   display resetview

   mol off 3
   mol off 6

}
