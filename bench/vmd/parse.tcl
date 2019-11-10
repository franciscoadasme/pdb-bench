set pdbfile [lindex $argv 0]
set repeats [lindex $argv 1]

set output [time {mol new $pdbfile} $repeats]
puts "elapsed = [expr {double([lindex $output 0]) * 1e-6}]"

exit
