require "chem"

pdbfile = ARGV[0]
repeats = ARGV[1].to_i

st = Chem::Structure.read pdbfile

total_elapsed = (0...repeats).sum do
  Time.measure do
    r1 = st.dig 'A', 50
    r2 = st.dig 'A', 60
    min_dist = r1.each_atom.min_of do |a1|
      r2.each_atom.min_of { |a2| Chem::Spatial.distance a1, a2 }
    end
  end
end

puts total_elapsed.total_seconds / repeats
