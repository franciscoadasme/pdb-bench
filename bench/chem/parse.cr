require "chem"

pdbfile = ARGV[0]
repeats = ARGV[1].to_i

total_elapsed = (0...repeats).sum do
  Time.measure do
    Array(Chem::Structure).from_pdb Path[pdbfile]
  end
end

puts total_elapsed.total_seconds / repeats
