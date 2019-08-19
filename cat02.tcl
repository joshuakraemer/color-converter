package require math::constants
::math::constants::constants e

set source_X [lindex $argv 0]
set source_Y [lindex $argv 1]
set source_Z [lindex $argv 2]

# Liverpool dark room (Xiao & al 2011, Excel file)
#set source_white_X 97.31273
#set source_white_Y 100
#set source_white_Z 138.59590
#set ambient_L 22.92
#set F 0.9

# Liverpool dark room (Xiao & al 2015)
set source_white_X 98.0
set source_white_Y 100.0
set source_white_Z 139.7
#set ambient_L 22.92
#set F 0.9

# Liverpool D65 illumination (Xiao & al 2011)
#set source_white_X 97.01096
#set source_white_Y 100
#set source_white_Z 135.10462
#set ambient_L 41.3
#set F 1.0

# D65 reference white (2° observer)
set target_white_X 95.0156
set target_white_Y 100.0
set target_white_Z 108.8199


# worked example
# https://books.google.de/books?id=0LVMLOSEeqoC
# target X = 18.25, Y = 24.66, Z = 27.94

#set ambient_L 200
#set F 1
#set source_white_X 98.88
#set source_white_Y 90
#set source_white_Z 32.03
#set target_white_X 100
#set target_white_Y 100
#set target_white_Z 100
#set source_X 19.31
#set source_Y 23.93
#set source_Z 10.14


# calculate degree of adaptation

#set D [expr {$F*(1 - (1/3.6) * $e ** ((-$ambient_L - 42)/92.0))}]
set D 1
#puts $D

# convert to the spectrally sharpened CAT02 LMS space
# matrices: http://rit-mcsl.org/fairchild/PDFs/AppearanceLec.pdf

proc XYZ_to_LMS {X Y Z} {
	set L [expr {0.7328*$X + 0.4296*$Y - 0.1624*$Z}]
	set M [expr {-0.7036*$X + 1.6975*$Y + 0.0061*$Z}]
	set S [expr {0.0030*$X + 0.0136*$Y + 0.9834*$Z}]
	return [list $L $M $S]
}

proc LMS_to_XYZ {L M S} {
	set X [expr {1.0961*$L - 0.2789*$M + 0.1827*$S}]
	set Y [expr {0.4544*$L + 0.4735*$M + 0.0721*$S}]
	set Z [expr {-0.0096*$L - 0.0057*$M + 1.0153*$S}]
	return [list $X $Y $Z]
}

set source_LMS [XYZ_to_LMS $source_X $source_Y $source_Z]
set source_L [lindex $source_LMS 0]
set source_M [lindex $source_LMS 1]
set source_S [lindex $source_LMS 2]

set source_white_LMS [XYZ_to_LMS $source_white_X $source_white_Y $source_white_Z]
set source_white_L [lindex $source_white_LMS 0]
set source_white_M [lindex $source_white_LMS 1]
set source_white_S [lindex $source_white_LMS 2]

set target_white_LMS [XYZ_to_LMS $target_white_X $target_white_Y $target_white_Z]
set target_white_L [lindex $target_white_LMS 0]
set target_white_M [lindex $target_white_LMS 1]
set target_white_S [lindex $target_white_LMS 2]


# CAT02 transformation

set target_L [expr {(($D*$source_white_Y*$target_white_L)/($target_white_Y*$source_white_L) + 1 - $D) * $source_L}]
set target_M [expr {(($D*$source_white_Y*$target_white_M)/($target_white_Y*$source_white_M) + 1 - $D) * $source_M}]
set target_S [expr {(($D*$source_white_Y*$target_white_S)/($target_white_Y*$source_white_S) + 1 - $D) * $source_S}]


# convert to XYZ space

set target_XYZ [LMS_to_XYZ $target_L $target_M $target_S]

puts $target_XYZ
