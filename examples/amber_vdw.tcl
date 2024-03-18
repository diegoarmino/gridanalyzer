
# READAMBERVDWARAMS: This routine parses all the AMBER NONBONDED (van der Waals/
# Lennard-Jones) parameters from a list of AMBER parameter files.
#
# Example: read_parm customparams1.par customparams2.par


proc read_parm {arg} {
global typelist rminlist epslist

## Read the file
  set opened [open "amber_vdw.parms"]
  set content [read $opened]
  close $opened

## Split into records on newlines
  set records [split $content "\n"]
#  puts Registros_del_archivo	
#  puts $records	
#  puts fin_del_archivo	

## Iterate over the records
   foreach line $records {
    
#     puts analizando_fila 
#     puts $line	
     set list [split $line]
     set data {}
     foreach var $list {
	if {$var != "" } {lappend data $var}
     }
     set type [lindex $data 0]
     set rmin [lindex $data 1]
     set eps [lindex $data 2]
    
#     puts $type
#     puts $rmin
#     puts $eps
     lappend typelist $type
     lappend rminlist $rmin  
     lappend epslist $eps
   }
   puts fin_del_ciclo_sobre_parms
#   puts $typelist
#   puts $rminlist
#   puts $epslist	
}


# ASSIGNAMBEPARMS: This routine sets the beta and occupancy fields of every atom
# in the specified molecule to what is needed forrunning the implicit ligand sampling 
# analysis. "read_parms" must have been previously run on the desired parameter 
# files.
#
# Example: assign_vdw 0

proc assign_vdw {molid} {
global typelist rminlist epslist

#  puts $molid
  set atomtypes [[atomselect $molid all] get type]
#  puts Atom_types	
#  puts $atomtypes 

  set atomradii {}
  set atomepsilon {}

  foreach tp $atomtypes { 
      set index [lsearch $typelist $tp]
      lappend atomradii [lindex $rminlist $index]
      lappend atomepsilon  [expr {-1 * [lindex $epslist $index]}]
      if {$index == -1 } {
      set index 0
      puts $tp not_found }
       	
  }   		

#  puts $atomradii
#  puts $atomepsilon  
  [atomselect $molid all] set radius $atomradii
  [atomselect $molid all] set occupancy $atomepsilon
#  [atomselect $molid all] set beta $atomepsilon
}


