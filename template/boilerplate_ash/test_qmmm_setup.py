"""
Test ash read Amber Parmtop and inpcord in QM/MM setting


Notes: atom indices as obtained from previous attempt to set up QM/MM with CP2K

&QM_KIND N
    MM_INDEX 1667 2039 2051 2103 2161 2169 2172 3119 3613 3671 3679 3682
&END QM_KIND
&QM_KIND H
    MM_INDEX 1668 1670 1672 1673 1675 1676 2040 2042 2044 2045 2047 2048 2052 2053 2104 2106 2108 2109 2111 2112 2162 2164 2166 2167 2171 2173 2175 3120 3122 3124 3125 3127 3128 3614 3616 3618 3619 3621 3622 3672 3674 3676 3677 3681 3683 3685 8158 8159 8161 8162 8164 8165
&END QM_KIND
&QM_KIND C
    MM_INDEX 1669 1671 1674 1677 1680 2041 2043 2046 2049 2054 2105 2107 2110 2113 2116 2163 2165 2168 2170 2174 2176 3121 3123 3126 3129 3132 3615 3617 3620 3623 3626 3673 3675 3678 3680 3684 3686
&END QM_KIND
&QM_KIND O
    MM_INDEX 1678 1679 1681 2050 2055 2114 2115 2117 2177 3130 3131 3133 3624 3625 3627 3687 8157 8160 8163
&END QM_KIND
&QM_KIND FE
    MM_INDEX 8155 8156
&END QM_KIND

"""
from ash import *

#===============================================================================
# Global variables
#===============================================================================

INPCRDFILE = r"/scratch/09069/dhp563/ash/test/amber/prmtop_inpcrd/prmtop_inpcrd_relaxed/6YD0_solv_tip3.inpcrd"
PRMTOPFILE = r"/scratch/09069/dhp563/ash/test/amber/prmtop_inpcrd/prmtop_inpcrd_relaxed/noflag.prmtop"

# QM atoms from old file
# QM_ATOMS = [1667, 2039, 2051, 2103, 2161, 2169, 2172, 3119, 3613, 3671, 3679, 3682, # Nitrogen
#             1668, 1670, 1672, 1673, 1675, 1676, 2040, 2042, 2044, 2045, 2047, 2048, #Hydrogen
#             2052, 2053, 2104, 2106, 2108, 2109, 2111, 2112, 2162, 2164, 2166, 2167,
#             2171, 2173, 2175, 3120, 3122, 3124, 3125, 3127, 3128, 3614, 3616, 3618,
#             3619, 3621, 3622, 3672, 3674, 3676, 3677, 3681, 3683, 3685, 8158, 8159,
#             8161, 8162, 8164,  #fix by removing 8165
#             1669, 1671, 1674, 1677, 1680, 2041, 2043, 2046, 2049, 2054, 2105, 2107, # Carbon
#             2110, 2113, 2116, 2163, 2165, 2168, 2170, 2174, 2176, 3121, 3123, 3126,
#             3129, 3132, 3615, 3617, 3620, 3623, 3626, 3673, 3675, 3678, 3680, 3684, 
#             3686, 
#             1678, 1679, 1681, 2050, 2055, 2114, 2115, 2117, 2177, 3130, 3131, 3133, # Oxygen
#             3624, 3625, 3627, 3687, 8157, 8160, 8163, 8156, #fix by adding 8156
#             8155, 8154] # Iron


#QM atoms from prep_qmmm_active_site.py
QM_ATOMS = [8156, 8157, 8158, 8159, 8160, 8161, 8162, 8163, 8164, 1670, 1671, 1672, 1673, 1674, 1675, 1676, 1677, 1678, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050, 2051, 2052, 2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2164, 2165, 2166, 2167, 2168, 2169, 2170, 2171, 2172, 2173, 2174, 3122, 3123, 3124, 3125, 3126, 3127, 3128, 3129, 3130, 3559, 3560, 3561, 3562, 3563, 3564, 3565, 3566, 3567, 3568, 3616, 3617, 3618, 3619, 3620, 3621, 3622, 3623, 3624, 3674, 3675, 3676, 3677, 3678, 3679, 3680, 3681, 3682, 3683, 3684, 8154, 8155, 8157, 8158, 8160, 8161, 8163, 8164]

#===============================================================================
# Functions 
#===============================================================================

# Test ran successfully with the provided files
def test_single_point(prmtopfile, inpcrdfile):
    """
    Read amber parmtop and inpcrd file as 
    an ASH fragment and run a single point calculation 
    with OpenMMTheory
    """
    #Amber files

    #Read coordinates from Amber INPCRD and PRMTOP FILES
    frag = Fragment(amber_prmtopfile=prmtopfile, amber_inpcrdfile=inpcrdfile)

    #Creating OpenMMobject using AMBER forcefield files
    #Note: Here we need to turn autoconstraints and rigidwater constraints off.
    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile=prmtopfile, printlevel=1,
        periodic=True, autoconstraints=None, rigidwater=False, platform='CUDA')

    #Run a simple energy+gradient job at the MM level to test whether everything is correct.
    Singlepoint(theory=openmmobject, fragment=frag)


# Test to set up sample QM/MM calculation
def test_QMMM(prmtopfile: str, inpcrdfile: str, qmatoms: list):
    #create ash fragment
    frag = Fragment(amber_prmtopfile = prmtopfile, 
                    amber_inpcrdfile = inpcrdfile)
    
    #create openmm object
    openmmobject = OpenMMTheory(Amberfiles = True, 
                                amberprmtopfile = prmtopfile,
                                printlevel = 1,
                                periodic = True,
                                autoconstraints = None,
                                rigidwater = False,
                                platform = 'CUDA')
    
    #Define QM region with atom indices, starting from 0
    qmatomlist = qmatoms

    #default from tutorial, suitable for an expensive system with transition metals
    #TPSSh functional, RIJCOSX scf strategy, D3BJ dispersion, SARC/J ZORA relativistic effects, ZORA tightscf slowconv
    ORCAinpline = "! TPSSh RIJCOSX  D3BJ SARC/J ZORA-def2-SVP ZORA tightscf slowconv"
    ORCAblocklines = """
    %maxcore 2000
    %scf
    MaxIter 1000
    end
    """

    # QM-region: Charge and multiplicity
    # from git discussion, this might not work as intended since
    # there might be a neutral hard setting to maintain consistency between QM and MM region
    # test this with the default p-multiplicity, charge = 2 multiplicity = 1
    charge = 0
    multiplicity = 1

    print("Creating ORCA QM object")
    # Create ORCA QM object
    orcaobject = ORCATheory(orcasimpleinput = ORCAinpline,
                            orcablocks = ORCAblocklines,
                            numcores = 8, 
                            orcadir = r"/work/09069/dhp563/ls6/orcaa/all_three")
    
    print("Creating QM/MM object")
    # Create QM/MM object
    qmmmobject = QMMMTheory(qm_theory = orcaobject,
                            mm_theory = openmmobject,
                            fragment = frag,
                            embedding = "Elstat",
                            qmatoms = qmatomlist,
                            printlevel = 2)                                          

    print("Running single point calculation")
    # Single-point job to test QM/MM setup
    Singlepoint(theory = qmmmobject,
                fragment = frag,
                charge = charge,
                mult = multiplicity)


#=============================
# Main
#=============================

if __name__ == "__main__":
    #test_single_point(PRMTOPFILE, INPCRDFILE)
    test_QMMM(PRMTOPFILE, INPCRDFILE, QM_ATOMS)




# #Read coordinates from either an XYZ-file, a PDB-file, or an ASH-file (.ygg)
# frag = Fragment(xyzfile="system.xyz", conncalc=False)

# #Creating OpenMMobject using CHARMM forcefield files
# #Note: Here we need to turn autoconstraints and rigidwater constraints off.
# openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
#     charmmprmfile=parfile, autoconstraints=None, rigidwater=False)

# #Forcefield files
# forcefielddir="/home/bjornsson/ASH-vs-chemshell-protein/QM-MM/FeMoco-test1/forcefielddir/"
# topfile=forcefielddir+"top_all36_prot.rtf"
# parfile=forcefielddir+"par_all36_prot.prm"
# psffile=forcefielddir+"new-XPLOR-psffile.psf"

# #Define QM region
# #IMPORTANT: Atom indices start at 0 in ASH.
# # Define either as lists in script:
# #qmatoms = [0, 5, 6, 7, 8]
# #Or read in list from file called: qmatoms (atom indices separated by space)
# qmatomlist = read_intlist_from_file("qmatoms")

# #Define QM-theory. Here ORCA
# ORCAinpline="! TPSSh RIJCOSX  D3BJ SARC/J ZORA-def2-SVP ZORA tightscf slowconv"
# ORCAblocklines="""
# %maxcore 2000
# %scf
# MaxIter 500
# end
# """

# #QM-region: Charge and multiplicity
# charge=-5
# mult=4

# #Create ORCA QM object
# orcaobject = ORCATheory(orcasimpleinput=ORCAinpline,
#                         orcablocks=ORCAblocklines, numcores=8)

# # Create QM/MM OBJECT
# qmmmobject = QMMMTheory(qm_theory=orcaobject, mm_theory=openmmobject,
#     fragment=frag, embedding="Elstat", qmatoms=qmatomlist, printlevel=2)

# # Single-point job to test QM/MM setup
# Singlepoint(theory=qmmmobject, fragment=frag, charge=charge,mult=mult)
