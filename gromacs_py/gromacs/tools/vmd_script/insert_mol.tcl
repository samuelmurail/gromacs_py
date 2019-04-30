#/Applications/VMD\ 1.8.6.app/Contents/vmd/vmd_MACOSXX86 -dispdev text -e add_alcool.tcl -args chain pdb_out 

proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}



[atomselect top "all"] set beta 0
set mol_insert_chain [lindex $argv 0]
set mol_length [lindex $argv 1]

set cutoff_water_clash 1.0
set cutoff_prot_off 7.0
set cutoff_mol_off 9.0
set cutoff_prot_in 10.0


set mol [atomselect top "chain $mol_insert_chain and noh" ]
set mol_in [$mol get residue]
set uniqueResid_mol [lsort -unique $mol_in]

set mol_num [ $mol num ]
puts $mol_num

for { set i 0 } { $i < $mol_num } { incr i $mol_length} {
	
	puts " $i / $mol_num"
	
    set water_residue [[atomselect top "beta 0 and water and (within $cutoff_prot_in of protein and not chain $mol_insert_chain) and not (within $cutoff_prot_off of protein and not chain $mol_insert_chain) and not (within $cutoff_mol_off of chain $mol_insert_chain) "] get residue]
	set water_mv [atomselect top "water and residue [lindex $water_residue 0]"]
    puts " [lindex $water_residue 0]"
	$water_mv set beta 1.0
    set res [lindex $uniqueResid_mol $i]
	set mol_mv [atomselect top "chain $mol_insert_chain and residue [expr $res] to [expr $res + $mol_length-1] "]
	#puts " [$water_mv num] [$mol_mv num]"

	set coor_water [center_of_mass $water_mv]
	set coor_mol [center_of_mass $mol_mv]
	
	#puts [lindex $coor_water 0]

	set vec_1 [list  [expr [lindex $coor_mol 0] -[lindex $coor_water 0]] [expr [lindex $coor_mol 1]- [lindex $coor_water 1]] [expr [lindex $coor_mol 2]- [lindex $coor_water 2]] ]
	set vec_2 [list  [expr [lindex $coor_water 0]- [lindex $coor_mol 0]]  [expr [lindex $coor_water 1] -[lindex $coor_mol 1]]  [expr [lindex $coor_water 2]- [lindex $coor_mol 2]] ]

	$water_mv moveby $vec_1
	$mol_mv moveby $vec_2

	$water_mv delete
	$mol_mv delete
}


puts  $uniqueResid_mol

[atomselect top "not beta = 1.0 and not same residue as (water within $cutoff_water_clash of chain $mol_insert_chain)" ]  writepdb [lindex $argv 2]

quit
