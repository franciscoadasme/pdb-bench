proc ramachandran {} {
    [atomselect top "all"] get {phi psi}
}

set pdbfile [lindex $argv 0]
set repeats [lindex $argv 1]

mol new $pdbfile

set output [time {ramachandran} $repeats]
puts "elapsed = [expr {double([lindex $output 0]) * 1e-6}]"

exit
