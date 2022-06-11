
from protein_docking import load_protein, create_system, add_restraints, \
                            output_coordinates, create_moves, angstrom, \
                            degrees, celsius, lam_clj, lam_restraint, \
                            wrap

# Load the first protein - the file(s) should contain a single molecule
# that contains the parameterised protein. Lots of file formats are supported.
protein0 = load_protein(["example/thrombin.prmtop", "example/thrombin.rst"])

# Load the second protein - the file(s) should contain a single molecule
# that contains the parameterised protein. Lots of file formats are supported.
protein1 = load_protein(["example/thrombin.prmtop", "example/thrombin.rst"])

# Now create a Sire system that contains the two proteins, and that
# can calculate the (soft core) intermolecular coulomb and LJ energy
# between those proteins
system = create_system(protein0, protein1)

# Now add some restraints as defined in the passed restraints file
system = add_restraints(system, "example/restraints.txt")

# Print out details of the system
print(system)
print(system.energies())

# Save the initial coordinates to a PDB file
output_coordinates(system, "output/0000.pdb")

# Now create some Monte Carlo moves that operate at the 
# specified temperature, and that perform maximum moves
# as specified. This will just rotate the proteins...
# (we use a high temperature to force a lot of sampling)
moves = create_moves(system,
                     max_translate=0*angstrom,
                     max_rotate=5*degrees,
                     temperature=50*celsius)

print(system.property("alpha"))
print(system.property("shiftDelta"))

# First start off with rotating the protein to minimise
# the restraint force. Do this by turning off the 
# clj forcefield
system.setComponent(lam_restraint, 1.0)
system.setComponent(lam_clj, 0.0)

# Do 10 cycles of Monte Carlo...
for i in range(1, 11):
    print(i, moves)
    system = moves.move(system, 100, False)
    print(system.energy())
    output_coordinates(system, "output/%04d.pdb" % i)

# Now let the proteins translate so that they can get
# closer to each other. Turn on the CLJ forcefield
# with 50% of the strength of the restraints. Also 
# turn up the soft core parameters so that the proteins
# are 'fuzzy', preventing us getting locked in 
# local minima

moves = create_moves(system,
                     max_translate=0.5*angstrom,
                     max_rotate=5*degrees,
                     temperature=50*celsius)

system.setComponent(lam_clj, 0.5)
system.setProperty("alpha", wrap(0.5))
system.setProperty("shiftDelta", wrap(2.0))

# Do 10 more cycles of Monte Carlo
for i in range(11, 21):
    print(i, moves)
    system = moves.move(system, 100, False)
    print(system.energy())
    output_coordinates(system, "output/%04d.pdb" % i)

# turn down the sampling for a more refined fit?
moves = create_moves(system,
                     max_translate=0.5*angstrom,
                     max_rotate=1.0*degrees,
                     temperature=35*celsius)

system.setComponent(lam_clj, 1.0)
system.setProperty("alpha", wrap(0.5))
system.setProperty("shiftDelta", wrap(2.0))

# Do 10 more cycles of Monte Carlo
for i in range(21, 31):
    print(i, moves)
    system = moves.move(system, 100, False)
    print(system.energy())
    output_coordinates(system, "output/%04d.pdb" % i)

# turn down fully and use a fully hard CLJ interaction
# for final refinement?
moves = create_moves(system,
                     max_translate=0.5*angstrom,
                     max_rotate=0.25*degrees,
                     temperature=25*celsius)

system.setComponent(lam_clj, 1.0)
system.setProperty("alpha", wrap(0.0))
system.setProperty("shiftDelta", wrap(0.0))

# Do 20 more cycles of Monte Carlo
for i in range(31, 51):
    print(i, moves)
    system = moves.move(system, 100, False)
    print(system.energy())
    output_coordinates(system, "output/%04d.pdb" % i)
