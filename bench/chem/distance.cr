require "chem"

pdbfile = ARGV[0]
repeats = ARGV[1].to_i

struc = Chem::Structure.read pdbfile

total_elapsed = (0...repeats).sum do
  Time.measure do
    r1 = struc.dig 'A', 50
    r2 = struc.dig 'A', 60
    min_dist = r1.atoms.min_of do |a1|
      r2.atoms.min_of do |a2|
        a1.coords.distance a2.coords
      end
    end
  end
end

puts total_elapsed.total_seconds / repeats
