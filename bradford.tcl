set source_X [lindex $argv 0]
set source_Y [lindex $argv 1]
set source_Z [lindex $argv 2]

# Liverpool (Xiao & al 2011) dark room
set source_white_X 97.31273
set source_white_Y 100
set source_white_Z 138.59590

# D65 reference white (2° observer)
set target_white_X 95.0156
set target_white_Y 100
set target_white_Z 108.8199

proc XYZ_to_RGB {X Y Z} {
	set R [expr {0.8951*$X/$Y + 0.2664 - 0.1614*$Z/$Y}]
	set G [expr {-0.7502*$X/$Y + 1.7135 + 0.0367*$Z/$Y}]
	set B [expr {0.0389*$X/$Y - 0.0685 + 1.0296*$Z/$Y}]
	return [list $R $G $B]
}


set source_RGB [XYZ_to_RGB $source_X $source_Y $source_Z]
set source_R [lindex $source_RGB 0]
set source_G [lindex $source_RGB 1]
set source_B [lindex $source_RGB 2]

set source_white_RGB [XYZ_to_RGB $source_white_X $source_white_Y $source_white_Z]
set source_white_R [lindex $source_white_RGB 0]
set source_white_G [lindex $source_white_RGB 1]
set source_white_B [lindex $source_white_RGB 2]

set target_white_RGB [XYZ_to_RGB $target_white_X $target_white_Y $target_white_Z]
set target_white_R [lindex $target_white_RGB 0]
set target_white_G [lindex $target_white_RGB 1]
set target_white_B [lindex $target_white_RGB 2]

set adapted_R [expr {$target_white_R*$source_R/$source_white_R}]
set adapted_G [expr {$target_white_G*$source_G/$source_white_G}]
set p [expr {($source_white_B/$target_white_B)**0.0834}]
set adapted_B [expr {$target_white_B*($source_B/$source_white_B)**$p}]

set adapted_X [expr {0.9870*$adapted_R*$source_Y - 0.1471*$adapted_G*$source_Y + 0.1600*$adapted_B*$source_Y}]
set adapted_Y [expr {0.4323*$adapted_R*$source_Y + 0.5184*$adapted_G*$source_Y + 0.0493*$adapted_B*$source_Y}]
set adapted_Z [expr {-0.0085*$adapted_R*$source_Y + 0.0400*$adapted_G*$source_Y + 0.9685*$adapted_B*$source_Y}]

puts "$adapted_X $adapted_Y $adapted_Z"
