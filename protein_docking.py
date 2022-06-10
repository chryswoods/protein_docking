
# Try to load the new Sire
# (this ensures forwards compatibility)
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass


from Sire.Units import angstrom, degrees, celsius
from Sire.CAS import Symbol as _Symbol

lam_clj = _Symbol("lambda_clj")
lam_restraint = _Symbol("lambda_restraint")


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
    return mols


def create_system(protein0, protein1, 
                  coulomb_cutoff=7.5*angstrom,
                  lj_cutoff=7.5*angstrom,
                  delta=15*angstrom):
    """Build a system for the two proteins"""

    from Sire.System import System
    from Sire.MM import InterGroupFF, CLJSoftShiftFunction
    from Sire.Mol import MGIdx, MoleculeGroup
    from Sire.Maths import Vector

    group0 = MoleculeGroup("protein0")
    group1 = MoleculeGroup("protein1")

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

    group0.add(protein0)
    group1.add(protein1)

    # and add this to the system
    system.add(ff)
    system.add(group0)
    system.add(group1)

    system.setProperty("alpha", 0.0)
    system.setProperty("shiftDelta", 0.0)

    return system


def add_restraints(system, restraints_file):
    """Add the restraints specified in 'restraints_file'
       to the system
    """

    from Sire.System import System
    from Sire.Units import kcal_per_mol
    from Sire.Mol import MolIdx, MGIdx, AtomName, ResNum, MGName
    from Sire.MM import RestraintFF, DistanceRestraint
    from Sire.FF import FFIdx

    # do everything in a copy of the System
    system = System(system)

    # forcefield to hold all of the restraints
    ff = RestraintFF("restraints")

    # Get a handle on the intermolecular forcefield
    cljff = system.forceFields()[FFIdx(0)]

    # symbol that represents the distance between the atoms
    r = DistanceRestraint.r()

    protein0 = system[MGName("protein0")][MolIdx(0)]
    protein1 = system[MGName("protein1")][MolIdx(0)]

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
            atom0 = protein0[ResNum(resnum0) + AtomName(atomname0)]
            atom1 = protein1[ResNum(resnum1) + AtomName(atomname1)]

            restraint = DistanceRestraint(atom0, atom1, func)

            ff.add(restraint)

    system.add(ff)

    system.setComponent(lam_clj, 1.0)
    system.setComponent(lam_restraint, 1.0)

    nrg = lam_clj * cljff.components().total() + lam_restraint * ff.components().total()

    system.setComponent(system.totalComponent(), nrg)

    return system


def output_coordinates(system, filename):
    """Write the passed system to the passed file"""

    from Sire.IO import MoleculeParser
    import os

    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

    MoleculeParser.save(system, filename)


def create_moves(system,
                 max_translate=0.1*angstrom,
                 max_rotate=5*degrees,
                 temperature=25*celsius):
    """Create the Monte Carlo moves for the passed
       system, which translates and rotates the two
       proteins by a maximum of 'translate_delta'
       and 'rotate_delta' respectively
    """

    from Sire.Move import RigidBodyMC, WeightedMoves
    from Sire.Mol import MGName

    protein0 = system[MGName("protein0")]
    protein1 = system[MGName("protein1")]

    move0 = RigidBodyMC(protein0)
    move1 = RigidBodyMC(protein1)

    for move in [move0, move1]:
       move.setMaximumTranslation(max_translate)
       move.setMaximumRotation(max_rotate)
       move.setTemperature(temperature)

       # if the protein is composed of multiple molecules,
       # then ensure that these are all rotated and translated
       # as one
       move.setSynchronisedTranslation(True)
       move.setSynchronisedRotation(True)
    
    moves = WeightedMoves()
    moves.setCheckRunningTotal(False)
    moves.add(move0, 1)
    moves.add(move1, 1)

    return moves


