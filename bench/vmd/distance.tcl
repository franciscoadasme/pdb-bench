proc distance {} {
    set r1 [atomselect top "resid 50"]
    set r2 [atomselect top "resid 60"]
    set min_d 999999
    foreach p1 [$r1 get {x y z}] {
        foreach p2 [$r2 get {x y z}] {
            set d [veclength [vecsub $p1 $p2]]
            if {$d < $min_d} {
                set min_d $d
            }
        }
    }
}

set pdbfile [lindex $argv 0]
set repeats [lindex $argv 1]

mol new $pdbfile

set output [time {distance} $repeats]
puts "elapsed = [expr {double([lindex $output 0]) * 1e-6}]"

exit
