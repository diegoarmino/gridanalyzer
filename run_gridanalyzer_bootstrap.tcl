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
# mode=3 ---> Results Visualization mode.
# mode=4 ---> Bootstrap error analysis.
# mode=5 ---> KDE of Bootstrap results.
#
# GAUSSIAN KERNEL DENSITY ESTIMATION (KDE):
# KDE is a method similar in purpuse to a histogram, so it can be used to obtain
# probability distribution functions (in this case in three dimensions), and
# so, by using the known formula G(x,y,z)=-RT ln( p(x,y,z) ) it is also possible to obtain
# the equivalent free energy. The advantage of this method compared to a simple
# histogram is that it can be shown that it converges much faster than the former
# (i.e. produces less noisy results), and doesnt depend on a given binning. A 2D histogram
# can be thought of as constructing a frequency plot by asigning a square to each
# event and just stacking squares in the corresponding bin. In KDE we use gaussian functions
# instead of squares and use no bins. So each gaussian function is evaluated in its right
# place, and are sumed toghether to get the final frequency plot (and normalized to get
# the probability distribution. There is one adjustable parameter though, the gaussian
# function width, sometimes called sigma, in reference to the normal distribution
# function tradicional nomenclature. There are ways to judiciously select the best sigma
# for each data set. But we have not yet implemented any of those methods.
# The script to runing KDE has two parts. Part 1 consist of a preparative step, in
# which a selection is made (for example all oxygen atoms of water molecules near the protein),
# and an intermediate file containing all the xyz coordinates of such water molecules for all
# frames of the trajectory is produced, alongside the coordinates of the corners of the
# smallest box containing all such selected water molecules.
# Part 2: The intermediate file is fed to the executable for the actual KDE analysis, which
# was written in fortran for performance reasons. This program calculates the KDE probability
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
# water sites, in the case the dx file is a free energy function obtained from KDE analysis of a
# trajectory. The program produces a pdb file contining all minima found with their corresponging
# energies in the "occupancy" slot of the pdb file, so they can be uploaded into VMD and colored
# using "occupancy" colour method in order to get a graphical representation of the positions and
# energies of the minima.
#
# MIGRATION PATHWAYS SEARCH
# This functionality searches for minimum energy pathways between the aforementioned minima using an
# algorithm similar to the Elastic Band methods. This is, an initial pathway is drawn as an interpolation
# between the positions of two minima, and then they are subject to a restrained optimization, in
# such a way that each intermediate point can only move on a plane (or more accurately, a thin slab)
# perpendicular to the line passing through both minima.
# The program produces a pdb file containing all pathway points found in pdb format.
#
# BOOTSTRAP ERROR ANALYSIS
# In this mode a resampling (with replacement) of the original trajectory must previously been done. To so a companion script is provided. Given an original MD trajectory, a resampled trajectory is one that is obtained from it by taking the same number of frames at random from the original includding repeated frames as a possibility. About 1000 resamplings are recomended.
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
                     set mode 5

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###################################################################
#           VARIABLE DEFINITION
###################################################################

# PDB file
set pdbfile_ini "2217_10C_frame0.pdb"

# DX file
set dxfile_ini "resamp_2217_1.dx"
#set dxfile_ini "free_energy.dx"

# Define reference coordinates
#set ref_coord_list [list {29.450000 35.230998 30.190000} {29.995001 24.528000 30.538000}]
#set ref_coord_list [list {28.209072 35.922890 30.807768} {33.612045 35.186863 32.503319} {24.272484 37.532246 37.002148} {29.128218 32.173157 26.666540} {24.325396 34.537834 26.891190} {30.933334 19.255215 28.709063} {30.018442 20.611382 35.833599} {23.951830 22.161655 27.178566} {34.786106 33.545963 40.964973}]
set ref_coord_list [list {29.43 35.23 30.19}]

#{27.416798 20.137243 42.200703}] 

# Set search radius. The search radius defines the volume to search for minima in, for the previoues
# reference coordinates. The optimizer won't consider grid points outside this cutoff.
set radius 20.0

# Define global optimization results file name (PDB). This file contains the coordinates and energy of
# each minima found.
#set opt_filename "wat_kde_global_opt.pdb"
set opt_filename "resamp_2217_1_opt.pdb"

# Define pathways filename (PDB). This file contains the pathways xyz coordinates as well as energies (except for
# minima points) for all the pathways found.
set pathways_filename "resamp_2217_1_optimized_path.pdb"

# Define the pruned .DX file name. This file is the same as the input .DX grid file, but with all values set to
# a very high energy value except for those close (up to a cutoff) to the migration pathways.
set pruned_grid_filename "resamp_2217_1_pruned_grid.dx"

# Set the index pair list for pathway search. It should be defined as a list of pairs with the following format:
#set start_indx_list [list  \
#    {4 0} \
#    {9 10} {10 5} {5 3} \
#]
#set start_indx_list [list \
#	]
###################################################################
#           VARIABLES FOR BOOTSTRAPING (MODE 4)
###################################################################
set rsmp_n 1153
set rsmp_dx_prefix "resamp_2217"

###################################################################
#           VARIABLES FOR BOOTSTRAPING (MODE 5)
###################################################################
# Selection for KDE analysis
set selection "all"
# Set path to kde fortran binary executable.
set path_KDE "/home/diegoa/dev/gridanalyzer/kdeanalysis.bin"
# Set the step of the grid for the sampling of the probability function in KDE analysis.
set kde_delta 0.5
# Set the width of the gaussian kernel functions for KDE analysis.
# It is best to compute this number using a nonparametric method suitable
# for multimodal distributions such as Improved Sheater Jones.
# ISJ method is implemented in KDEpy library.
set kde_sigma 2.2 
# Define file name of input pdb file.
set in_pdb "all.pdb"
# Define file name of the input for KDE binary.
set out "kde_input"

###################################################################
#           SET KDE PARAMETERS AND INITIAL VARIABLES
###################################################################
#
# Selection for KDE analysis
#set selection "noh and water and within 3 of protein or resname HEM"
# Path to KDE analysis executable
#set path_KDE "/home/diegoa/dev/ilsanalyzer/kdeanalysis.bin"
# Set the step of the grid for the sampling of the probability function in KDE analysis.
#set kde_delta 0.5
# Set the width of the gaussian kernel functions for KDE analysis
#set kde_sigma 2.2
# Set trajectory and prmtop filenames
#set top {2217_con_eme.prmtop}
#set traj {no_acqua_10_1_ALL_MD.nc}
#set pdb_in "all.pdb"
# Set intermidiate file name.
#set max_frames 10000
# KDE analysis only stores a fraction of the trajectory on memory.
# Using the following variable set how big (in frames) this chunk
# of trajectory should be.
#set batch 5000

# KDE analysis is performed by a code written in fortran. This code needs
# a file with a list of all the x y z coordinates of the probe (probably water)
# molecules for the whole trajectory, plus a heading for the coorinates of
# the bounding box containing all these molecules. The following variable defines
# if it is necessary to build this file. If the file has allready been produced
# then set do_build_input to 0, as the process of building this file is time-consuming.
#
# Build input file for kdeanalyzer. 1=yes 0=no
#
#set do_build_input 0

# Set the number of lines of input file for kdeanalysis if known (produced in a previous run).
#set kde_line_cnt 3933787

###################################################################
#           SET LOAD RESULTS VARIABLES (MODE 3)
###################################################################
#
#set ils_opt_filename "global_opt.pdb"
#set ils_opt_filename "2217_10C_ils_opt.pdb"
#set ils_pathways_filename "2217_10C_ils_optimized_path.pdb"
#set ils_pruned_grid_filename "2217_10C_ils_pruned_grid.dx"
#
#set kde_opt_filename "wat_kde_global_opt.pdb"
#set kde_pathways_filename "wat_kde_optimized_path.pdb"
#set kde_pruned_grid_filename "wat_kde_pruned_grid.dx"
#
#
#set ils_start_indx_list [list \
#	{217 111} {111 112} {112 47} {112 134} {134 87} {87 88} {88 121} {121 120} {120 85} {85 902} {902 906} {906 933} {933 938} {85 98} {98 133} {133 111} \
#	{217 248} {248 274} {274 287} {287 286} {248 249} {249 286} {286 268} {268 312} {312 338} {268 276} {276 283} {283 316} {316 359} {359 314} {314 394} {394 393} {394 408} {408 352} {352 333} {333 285} {285 273} {273 257} {257 244} {244 216} {216 212} {212 207} {207 205} {205 206} {206 204} {206 1453} {206 1420} {1420 1446} {1446 1413} {1413 1411} {1411 1397} {1397 1399} {1453 1469} {1469 1486} {1453 1494} {1494 1512} {1512 1530} {1530 1541} {1541 281} {281 204} {248 240} {240 247} {247 232} {232 275} {283 269} {269 224} {224 195} {195 157} {157 285} {157 1095} {1095 193} {193 184} {184 1073} {1073 1100} {1100 1093} {1093 238} {238 1101} {238 155} {1101 246} {246 284} {284 273} {284 323} {157 233} {233 186} {186 168} {168 178} {178 187} {187 196} {196 195} {196 181} {181 150} {150 102} {324 349} {349 350} {350 357} {357 365} {365 396} {396 359} {207 217} \
#	]

#set ils_connect_list [list 1 {0 2} {1 3} {2 4} {3 5} {4 6} {5 7} 6 9 {8 10} {9 11} {10 12} {11 13} {12 14} 13 16 {15 17} {16 18} {17 19} {18 20} {19 21} 20 23 {22 24} {23 25} {24 26} 25 28 {27 29} {28 30} {29 31} 30 33 {32 34} {33 35} 34 37 {36 38} {37 39} {38 40} 39 42 {41 43} 42 45 {44 46} {45 47} {46 48} {47 49} {48 50} {49 51} {50 52} 51 54 {53 55} {54 56} {55 57} {56 58} {57 59} {58 60} {59 61} {60 62} 61 64 {63 65} {64 66} 65 68 {67 69} {68 70} 69 72 {71 73} {72 74} 73 76 {75 77} {76 78} {77 79} {78 80} {79 81} {80 82} {81 83} {82 84} {83 85} 84 87 {86 88} {87 89} {88 90} {89 91} 90 93 {92 94} {93 95} {94 96} 95 98 {97 99} {98 100} {99 101} {100 102} 101 104 {103 105} 104 107 {106 108} 107 110 {109 111} 110 113 {112 114} {113 115} {114 116} 115 118 {117 119} 118 121 {120 122} {121 123} {122 124} {123 125} {124 126} {125 127} 126 129 {128 130} 129 132 {131 133} {132 134} 133 136 {135 137} {136 138} {137 139} {138 140} {139 141} 140 143 {142 144} {143 145} 144 147 {146 148} {147 149} 148 151 {150 152} {151 153} {152 154} 153 156 {155 157} {156 158} {157 159} {158 160} {159 161} {160 162} 161 164 {163 165} {164 166} {165 167} {166 168} {167 169} 168 171 {170 172} {171 173} {172 174} {173 175} 174 177 {176 178} {177 179} {178 180} 179 182 {181 183} {182 184} {183 185} {184 186} 185 188 {187 189} {188 190} {189 191} 190 193 {192 194} {193 195} {194 196} 195 198 {197 199} {198 200} {199 201} {200 202} 201 204 {203 205} {204 206} {205 207} {206 208} {207 209} 208 211 {210 212} {211 213} {212 214} {213 215} 214 217 {216 218} {217 219} 218 221 {220 222} {221 223} {222 224} 223 226 {225 227} {226 228} {227 229} {228 230} 229 232 {231 233} {232 234} {233 235} {234 236} {235 237} 236 239 {238 240} {239 241} {240 242} {241 243} {242 244} 243 246 {245 247} {246 248} {247 249} {248 250} {249 251} 250 253 {252 254} {253 255} {254 256} {255 257} {256 258} {257 259} {258 260} {259 261} {260 262} {261 263} {262 264} 263 266 {265 267} {266 268} {267 269} {268 270} {269 271} 270 273 {272 274} 273 276 {275 277} {276 278} 277 280 {279 281} {280 282} {281 283} {282 284} {283 285} {284 286} {285 287} 286 289 {288 290} {289 291} {290 292} {291 293} 292 295 {294 296} 295 298 {297 299} {298 300} {299 301} {300 302} 301 304 {303 305} 304 307 {306 308} {307 309} {308 310} {309 311} {310 312} 311 314 {313 315} {314 316} {315 317} 316 319 {318 320} 319 322 {321 323} {322 324} {323 325} {324 326} 325 328 {327 329} {328 330} {329 331} {330 332} {331 333} 332 335 {334 336} {335 337} {336 338} {337 339} {338 340} {339 341} {340 342} {341 343} {342 344} 343 346 {345 347} 346 349 {348 350} {349 351} 350 353 {352 354} {353 355} 354 357 {356 358} {357 359} {358 360} {359 361} 360 363 {362 364} {363 365} {364 366} 365 368 {367 369} {368 370} {369 371} {370 372} {371 373} {372 374} 373 376 {375 377} 376 379 {378 380} {379 381} {380 382} {381 383} {382 384} {383 385} {384 386} {385 387} {386 388} {387 389} 388 391 {390 392} {391 393} {392 394} {393 395} {394 396} {395 397} {396 398} {397 399} {398 400} {399 401} {400 402} {401 403} {402 404} 403 406 {405 407} {406 408} {407 409} {408 410} {409 411} {410 412} {411 413} {412 414} 413 416 {415 417} 416 419 {418 420} {419 421} 420 423 {422 424} {423 425} {424 426} {425 427} {426 428} 427 430 {429 431} {430 432} 431 434 {433 435} 434 437 {436 438} 437 440 {439 441} 440 443 {442 444} {443 445} 444 447 {446 448} {447 449} {448 450} 449 452 {451 453} {452 454} {453 455} 454 457 {456 458} {457 459} {458 460} {459 461} {460 462} {461 463} 462 465 {464 466} {465 467} {466 468} {467 469} {468 470} {469 471} {470 472} 471 474 {473 475} {474 476} {475 477} {476 478} 477 480 {479 481} {480 482} {481 483} 482 485 {484 486} {485 487} {486 488} {487 489} 488 491 {490 492} {491 493} {492 494} 493 496 {495 497} 496 499 {498 500} {499 501} {500 502} 501 504 {503 505} {504 506} 505 508 {507 509} {508 510} 509 512 {511 513} 512 515 {514 516} {515 517} {516 518} 517 520 {519 521} {520 522} {521 523} {522 524} 523 526 {525 527} {526 528} {527 529} 528 531 {530 532} {531 533} {532 534} {533 535} {534 536} {535 537} {536 538} {537 539} 538 541 {540 542} 541 544 {543 545} {544 546} {545 547} {546 548} {547 549} {548 550} 549 552 {551 553} {552 554} 553 556 {555 557} {556 558} {557 559} {558 560} {559 561} {560 562} {561 563} {562 564} {563 565} 564]
#	
#
#
#set kde_start_indx_list [list  \
#    {22 10} {10 9} {9 14} {14 17} {14 13} {13 21} {21 66} {17 16} {16 24} {24 23} {23 20} {2 62} {2 19} \
#    {58 57} {57 68} {68 73} {68 61} {61 62} {62 3} {3 19} {19 27} {61 69} {27 45} {27 20} {27 61}
#]
#
#set kde_connect_list [list 1 {0 2} {1 3} {2 4} 3 6 {5 7} {6 8} {7 9} 8 11 {10 12} {11 13} {12 14} {13 15} {14 16} {15 17} 16 19 {18 20} {19 21} {20 22} {21 23} {22 24} 23 26 {25 27} {26 28} {27 29} {28 30} {29 31} {30 32} {31 33} 32 35 {34 36} {35 37} 36 39 {38 40} {39 41} {40 42} {41 43} {42 44} {43 45} {44 46} {45 47} 46 49 {48 50} {49 51} {50 52} 51 54 {53 55} {54 56} {55 57} 56 59 {58 60} {59 61} {60 62} {61 63} {62 64} {63 65} 64 67 {66 68} {67 69} 68 71 {70 72} {71 73} {72 74} 73 76 {75 77} {76 78} {77 79} 78 81 {80 82} {81 83} {82 84} {83 85} {84 86} {85 87} 86 89 {88 90} {89 91} {90 92} {91 93} {92 94} 93 96 {95 97} {96 98} {97 99} {98 100} {99 101} {100 102} 101 104 {103 105} {104 106} {105 107} {106 108} {107 109} {108 110} {109 111} {110 112} 111 114 {113 115} {114 116} {115 117} {116 118} {117 119} {118 120} 119 122 {121 123} {122 124} {123 125} {124 126} {125 127} {126 128} 127 130 {129 131} {130 132} {131 133} {132 134} 133 136 {135 137} {136 138} {137 139} {138 140} {139 141} {140 142} {141 143} {142 144} {143 145} {144 146} {145 147} {146 148} {147 149} {148 150} {149 151} {150 152} 151 154 {153 155} {154 156} {155 157} {156 158} {157 159} {158 160} {159 161} {160 162} {161 163} 162 165 {164 166} {165 167} {166 168} {167 169} {168 170} {169 171} {170 172} {171 173} 172 175 {174 176} {175 177} {176 178} {177 179} {178 180} {179 181} {180 182} {181 183} {182 184} {183 185} {184 186} {185 187} {186 188} {187 189} 188 191 {190 192} {191 193} {192 194} {193 195} {194 196} {195 197} {196 198} {197 199} {198 200} 199]
#

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
#  MODE 0: GAUSSIAN KERNEL DENSITY ESTIMATION ANALYSIS
###################################################################

if { $mode == 0 } {
   kde_analysis_preparation $selection $batch $max_frames $top $traj $out_pdb $do_build_input
   exec $path_KDE "$out_pdb.xyz" $kde_line_cnt $kde_delta $kde_sigma > "kde.log"
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
   global_search_min $ref_coord_list $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta $endgridX $endgridY $endgridZ $outfile
   quit
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

   set kde_minlist [list ]
   foreach a $kde_start_indx_list {
      foreach b $a {
         if { [ lsearch $kde_minlist $b ] < 0 } {
            lappend kde_minlist $b }
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
#   mol modselect 0 top "not index 104 111 113 14 6 74"
   

   mol new $ils_pruned_grid_filename type dx first 0 last -1 step 1 waitfor 1 volsets 0
   mol modstyle 0 top Isosurface 12.00 0 0 0 1 1
   mol modmaterial 0 top Transparent
   mol modcolor 0 top ColorID 7

   set top {/home/diegoa/inv/parma/2217/2217_con_eme.prmtop}
   set traj {/home/diegoa/inv/parma/2217/no_acqua_10_1_ALL_MD.nc}
   mol load parm7 $top
   mol addfile $traj step 20 first 0 last -1 waitfor all type netcdf
   mol selection "protein"
   mol representation NewCartoon
   mol color colorID 3
   mol modstyle 0 0 NewCartoon 0.300000 12.000000 2.120000 0
   mol addrep top
#
   mol selection "resname HEM"
   mol representation Licorice
   mol color Name
   mol addrep top
#
   mol selection "resid 75 37 88 48 21 125"
   mol representation Licorice
   mol color Name
   mol addrep top


#   mol new $kde_opt_filename type {pdb} first 0 last -1 step 1 waitfor 1
#   mol modstyle 0 top VDW 0.500000 20.000000
#   mol modcolor 0 top Occupancy
#   mol modselect 0 top "index $kde_minlist"
#   mol modmaterial 0 top AOChalky

#   mol new $kde_pathways_filename type {pdb} first 0 last -1 step 1 waitfor 1
#   mol modstyle 0 top Licorice 0.200000 12.000000 12.000000
#   mol modcolor 0 top Occupancy
#   set sel [atomselect top "all"]
#   $sel setbonds $kde_connect_list
#   mol modmaterial 0 top AOChalky
#   mol modselect 0 top "all not index 96 to 200 42 to 44"

#   mol new $kde_pruned_grid_filename type dx first 0 last -1 step 1 waitfor 1 volsets 0
#   mol modstyle 0 top Isosurface 12.00 0 0 0 1 1
#   mol modmaterial 0 top Transparent
#   mol modcolor 0 top ColorID 7
#
   # EXTRA REPRESENTATIONS

   # ILS
#   set top {/home/diegoa/inv/parma/2217/2217_con_eme.prmtop}
#   set traj {/home/diegoa/inv/parma/2217/no_acqua_10_1_ALL_MD.nc}
#   mol new parm7 $top
#   mol addfile $traj step 20 first 0 last -1 waitfor all type netcdf
#   mol selection "protein"
#   mol representation NewCartoon
#   mol color colorID 3
#   mol modstyle 0 0 NewCartoon 0.300000 12.000000 2.120000 0
#   mol color ColorID 7
#   mol representation VDW 0.500000 20.000000
#   mol selection index 25 4 13 10 14 11 6 5
#   mol material AOChalky
#   mol addrep 1
#   mol color ColorID 7
#   mol representation Licorice 0.200000 12.000000 12.000000
#   mol selection index 0 to 36
#   mol material AOChalky
#   mol addrep 2
#
#   mol color ColorID 4
#   mol representation VDW 0.500000 20.000000
#   mol selection index 21 36 19 27 45 84 69 74 73 61 62 26 37 38 7
#   mol material AOChalky
#   mol addrep 1
#   mol color ColorID 4
#   mol representation Licorice 0.200000 12.000000 12.000000
#   mol selection index 36 to 106
#   mol material AOChalky
#   mol addrep 2
#
#   mol color ColorID 11
#   mol representation VDW 0.500000 20.000000
#   mol selection index 29 35 28 39 31 56
#   mol material AOChalky
#   mol addrep 1
#   mol color ColorID 11
#   mol representation Licorice 0.200000 12.000000 12.000000
#   mol selection index 106 to 112 114 to 150
#   mol material AOChalky
#   mol addrep 2
#
#   mol color ColorID 1
#   mol representation VDW 0.700000 20.000000
#   mol selection index 25
#   mol material AOChalky
#   mol addrep 1

   # KDE

   display resetview

   mol off 3
   mol off 6

}

###################################################################
#  MODE 4: BOOTSTRAP ERROR ANALYSIS
###################################################################

if { $mode == 4 } {

   # Load pdb file with protein structure to VMD.
   load_pdb_to_vmd $pdbfile_ini

   for {set rsmp_i 982} { $rsmp_i <= $rsmp_n } { incr rsmp_i } {
      # Load dx file with ILS grid results to VMD.
      puts "ANALYZING RESAMPLE #$rsmp_i"
      set dxfile "${rsmp_dx_prefix}_${rsmp_i}.dx"
      set outfile "${rsmp_dx_prefix}_opt_${rsmp_i}.pdb"
      load_dx_to_vmd $dxfile

      # Read DX file grid to analyze.
      read_DX $dxfile

      global_search_min $ref_coord_list $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta $endgridX $endgridY $endgridZ $outfile
      puts "FINISHED RESAMPLE #$rsmp_i"

   }
}

###################################################################
#  MODE 5: GAUSSIAN KERNEL DENSITY ESTIMATION ANALYSIS OF BOOTSTRAP
###################################################################

if { $mode == 5 } {
   set kde_line_cnt [ bootstrap_kde_preparation $selection $in_pdb $out ]
   puts "Line count of input file"
   puts $kde_line_cnt
   exec $path_KDE "$out.xyz" $kde_line_cnt $kde_delta $kde_sigma > "kde.log"
}

