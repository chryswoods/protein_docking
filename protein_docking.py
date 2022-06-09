
# Try to load the new Sire
# (this ensures forwards compatibility)
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass


from Sire.Units import angstrom


def load_protein(files):
    """Load the protein from the passed files. This returns 
       the protein as a single molecule
    """
    from Sire.IO import MoleculeParser

    if type(files) is str:
        files = [files]

    mols = MoleculeParser.load(files).molecules()

    if mols.nMolecules() == 1:
        return mols[mols.molNums()[0]]

    # join all of the molecules into a single molecule
    print("NEED TO ADD FUNCTIONALITY TO JOIN MOLECULES")
    raise ProgramError("This doesn't support multi-molecule proteins")


def run():
    print("RUN")


def create_system(protein0, protein1, 
                  coulomb_cutoff=7.5*angstrom,
                  lj_cutoff=7.5*angstrom,
                  delta=15*angstrom):
    """Build a system for the two proteins"""

    from Sire.System import System
    from Sire.MM import InterGroupFF, CLJSoftShiftFunction
    from Sire.Mol import MGIdx    
    from Sire.Maths import Vector

    system = System("protein_docking")

    # Create a forcefield to calculate the intermolecular
    # interaction between the proteins. This is probably
    # best being a soft-core forcefield, so that it can
    # start soft, and then become progressively harder
    ff = InterGroupFF("cljff")

    func = CLJSoftShiftFunction(coulomb_cutoff, lj_cutoff)
    ff.setCLJFunction(func)


    # now translate the proteins so that they are separated
    # from each other
    box0 = protein0.evaluate().aaBox()
    box1 = protein1.evaluate().aaBox()

    # move the first protein to (0,0,0)
    protein0 = protein0.move().translate(-box0.center()).commit()

    # move the second protein to (X, 0, 0) where 
    # X is far enough away that the proteins don't overlap
    # (proteins are separated by 'delta'
    X = box0.radius() + box1.radius() + delta.to(angstrom)

    protein1 = protein1.move().translate((Vector(X,0,0)-box1.center())).commit()

    # now add the proteins to the forcefield
    ff.add(protein0, MGIdx(0))
    ff.add(protein1, MGIdx(1))

    # and add this to the system
    system.add(ff)

    return system


def add_restraints(system, restraints_file):
    """Add the restraints specified in 'restraints_file'
       to the system
    """

    from Sire.System import System
    from Sire.Units import kcal_per_mol
    from Sire.Mol import MolIdx, MGIdx, AtomName, ResNum
    from Sire.MM import RestraintFF, DistanceRestraint

    # do everything in a copy of the System
    system = System(system)

    # forcefield to hold all of the restraints
    ff = RestraintFF("restraints")

    # symbol that represents the distance between the atoms
    r = DistanceRestraint.r()

    with open(restraints_file, "r") as FILE:
        for line in FILE.readlines():
            line = line.lstrip().rstrip()

            if line.startswith("#"):  # comment
                continue

            words = line.split()

            # format is
            # atomname_0  residuenum_0  atomname_1  residuenum_1  R0/angstrom  k/kcal mol-1 angstrom-2
            if len(words) < 6:
                continue

            atomname0 = words[0]
            resnum0 = int(words[1])
            atomname1 = words[2]
            resnum1 = int(words[3])
            R0 = float(words[4]) * angstrom
            k = float(words[5]) * kcal_per_mol / (angstrom*angstrom)

            # define the restraint function - we will use a simple harmonic
            # bond, but anything could in theory be used
            func = k.value() * (r - R0.value())**2

            # find the two atoms
            protein0 = system[MGIdx(0)][MolIdx(0)].molecule()
            protein1 = system[MGIdx(1)][MolIdx(0)].molecule()

            atom0 = protein0[ResNum(resnum0) + AtomName(atomname0)]
            atom1 = protein1[ResNum(resnum1) + AtomName(atomname1)]

            restraint = DistanceRestraint(atom0, atom1, func)

            ff.add(restraint)

    system.add(ff)

    return system


def output_coordinates(system, filename):
    """Write the passed system to the passed file"""

    from Sire.IO import MoleculeParser
    import os

    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

    MoleculeParser.save(system, filename)

