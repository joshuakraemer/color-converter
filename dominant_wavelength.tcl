package require math::constants
::math::constants::constants pi

set x [lindex $argv 0]
set y [lindex $argv 1]

set whitepoint_x 0.312721
set whitepoint_y 0.329126

proc angle {x y} {
	variable pi
	variable whitepoint_x
	variable whitepoint_y
	set diff_x [expr {$x - $whitepoint_x}]
	set diff_y [expr {$y - $whitepoint_y}]
	set angle [expr {atan2($diff_y, $diff_x)*-1}]; # clockwise = positive
	if {$angle < 0} {set angle [expr $angle + 2*$pi]}
	return $angle
}

# angles are expressed relative to the minimal angle
proc relative_angle {x y reference_angle} {
	variable pi
	set relative_angle [expr {[angle $x $y] - $reference_angle}]
	if {$relative_angle < 0} {set relative_angle [expr $relative_angle + 2*$pi]}
	return $relative_angle
}

set coordinates [list]
set file_channel [open cc2012xyz2_fine_5dp.csv r]
while {-1 != [gets $file_channel line]} {
	lappend coordinates [split $line ,]
}
close $file_channel

set reference_angle [angle [lindex $coordinates 0 1] [lindex $coordinates 0 2]]
set max_angle [relative_angle [lindex $coordinates end 1] [lindex $coordinates end 2] $reference_angle]
set wanted_angle [relative_angle $x $y $reference_angle]

set sign {}
if {$wanted_angle > $max_angle} {
	set sign -
	set wanted_angle [expr {$wanted_angle - $pi}]
}


# search for the nearest angle

set lower_limit_index 0
set lower_limit_angle 0
set upper_limit_index [expr {[llength $coordinates] - 1}]
set upper_limit_angle $max_angle

while {($upper_limit_index - $lower_limit_index) > 1} {
	set test_index [expr {($upper_limit_index + $lower_limit_index)/2}]
	set test_angle [relative_angle [lindex $coordinates $test_index 1] [lindex $coordinates $test_index 2] $reference_angle]

	if {($test_angle > $wanted_angle)} {
		set upper_limit_index $test_index
		set upper_limit_angle $test_angle
	} elseif {($test_angle < $wanted_angle)} {
		set lower_limit_index $test_index
		set lower_limit_angle $test_angle
	} else {
		puts "$sign[lindex $coordinates $test_index 0]"
		exit
	}
}

set lower_limit_difference [expr {$wanted_angle - $lower_limit_angle}]
set upper_limit_difference [expr {$upper_limit_angle - $wanted_angle}]

if {$upper_limit_difference < $lower_limit_difference} {
	puts "$sign[lindex $coordinates $upper_limit_index 0]"
} else {
	puts "$sign[lindex $coordinates $lower_limit_index 0]"
}
