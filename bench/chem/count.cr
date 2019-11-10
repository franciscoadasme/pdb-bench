require "chem"

pdbfile = ARGV[0]
repeats = ARGV[1].to_i

st = Chem::Structure.read pdbfile

total_elapsed = (0...repeats).sum do
  Time.measure do
    st.each_residue.count &.name.==("ALA")
  end
end

puts total_elapsed.total_seconds / repeats
