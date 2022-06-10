
from protein_docking import load_protein, create_system, add_restraints, \
                            output_coordinates, create_moves

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

# Now create some Monte Carlo moves
moves = create_moves(system)

for i in range(1, 51):
    print(i, moves)
    system = moves.move(system, 100, False)
    print(system.energy())
    output_coordinates(system, "output/%04d.pdb" % i)
