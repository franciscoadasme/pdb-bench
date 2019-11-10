proc count {} {
    [atomselect top "resname ALA and name CA"] num
}

set pdbfile [lindex $argv 0]
set repeats [lindex $argv 1]

mol new $pdbfile

set output [time {count} $repeats]
puts "elapsed = [expr {double([lindex $output 0]) * 1e-6}]"

exit
