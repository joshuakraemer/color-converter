# source color under C illuminant

set x 0.102776
set y 0.102864
set Y 100

set X [expr {$x*$Y/double($y)}]
set Z [expr {(1-$x-$y)*$Y/$y}]

set X 96.8652
set Y 100
set Z 116.6144

puts "$X $Y $Z"

set X_D65 [expr {$X*0.9904476 - $Y*0.0071683 - $Z*0.0116156}]
set Y_D65 [expr {$X*-0.0123712 + $Y*1.0155950 - $Z*0.0029282}]
set Z_D65 [expr {$X*-0.0035635 + $Y*0.0067697 + $Z*0.9181569}]

puts "$X_D65 $Y_D65 $Z_D65"
