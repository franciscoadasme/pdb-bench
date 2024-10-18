proc distance {} {
    set sel [atomselect top "all"]
    set mat [measure fit $sel $sel]
    $sel move $mat
}

set pdbfile [lindex $argv 0]
set repeats [lindex $argv 1]

mol new $pdbfile

set output [time {distance} $repeats]
puts "elapsed = [expr {double([lindex $output 0]) * 1e-6}]"

exit
