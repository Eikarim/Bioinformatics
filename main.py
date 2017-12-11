
# FIRST WE IMPORT ALL THE LIBRARIES

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from ForceField import VdwParamset
from ResLib import  ResiduesDataLib
from simpleinteractions import simpleinteractions
from plot import plot
from energyinteractions import energies

import sys
import argparse
import os


############################################

# Then we need to characterize some parameters

COVLNK = 2.0
HBLNK  = 3.5

# We define all the atoms that can make an hydrogen bond
all_polars = [
    'N', 'ND1', 'ND2', 'NE',  'NE1', 'NE2', 'NH1', 'NH2', 'NZ',
    'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG',  'OG1', 'OH',
    'S', 'SD',  'SG'
]
backbone_polars =  ['N','O'] #And the ones that are in the backbone
waternames = ['WAT','HOH']

##########################################

# Now we use the Parser to create some extra arguments that can help using our function.
parser = argparse.ArgumentParser(
        prog='polarContacts',
        description='Polar contacts detector'
        )

parser.add_argument(
    '--backonly',
    action='store_true',
    dest='backonly',
    help='Restrict to backbone'  #If this argument == True, only backbone atoms are loaded
)

parser.add_argument(
        '--nowats',
        action='store_true',
        dest='nowats',   
        help='Exclude water molecules'
    )
    
parser.add_argument(
        '--diel',
        type= float,
        action='store',
        dest='diel',
        default = 1.0,  
        help='Relative dielectric constant'
    )
    
parser.add_argument(
        '--vdw',
        action='store',
        dest='vdwprm',
        help='VDW Parameters file'
    )
    
parser.add_argument(
        '--rlib',
        action='store',
        dest='reslib',
        help='AminoAcid library'
    )

parser.add_argument( ### This is a modification of the initial script, since this was wrong,
		'--pdb_path',
		action='store',
		dest="pdb_path",
		help = "PDB file")

parser.add_argument( ### We add an extra argument for if somebody want to use another MSMS. Although default MSMS file is attached to the folder.
        '--msms_path',
        action='store',
        dest="msms_path",
        default='./msms.MacOSX.2.6.1',  
        help= "MSMS file"
        )



args = parser.parse_args()

  

backonly = args.backonly  #We are assigning the args to a variable with the same name
nowats =args.nowats
pdb_path = args.pdb_path
vdwprm = args.vdwprm
reslib = args.reslib
diel = args.diel
msms_path = args.msms_path
    
# Load VDW parameters
vdwParams = VdwParamset(vdwprm)


# Load AA Library
aaLib = ResiduesDataLib(reslib)

    
if not pdb_path:
    parser.print_help()
    sys.exit(2)


### Until here we only loaded things.

parser = PDBParser(PERMISSIVE=1)   #We create a PDB Parser which creates the structure. Permissive means that the errors can be overlooked.
    
try:
    st = parser.get_structure('st', pdb_path)  #Here we put the argument pdb_path, for making it more dynamic.
except OSError:
    print ("#ERROR: loading PDB")
    sys.exit(2)

# Checking for models
if len(st) > 1:
    print ("#WARNING: Several Models found, using only first")

# We told the script to use the structure 1.
st = st[0]

# And we make a list of polar atoms
polats = []
if backonly:   #If backonly were in the arguments, only backbone atoms are used

    selected_atoms = backbone_polars   
else:
    selected_atoms = all_polars

for at in st.get_atoms():   #And we append each atom to the list
    if at.id in selected_atoms:   
        polats.append(at)   
            
      
#Now we will search for hydrogen-nitrogen interactions.
nbsearch = NeighborSearch(polats)   
hblist = []
for at1, at2 in nbsearch.search_all(HBLNK):  
    if at1.get_parent() == at2.get_parent():   #This is Kruskal, we are looking for not taking atoms of the same residue
        continue

    if (at1-at2) < COVLNK:		#Here we discard the atoms that make a covalent bond
        continue
    if abs(at2.get_parent().id[1] - at1.get_parent().id[1]) == 1:  #Here we see if one residue is neighbour to the other.
        continue 

    if nowats:  #And here if the argument nowats it's true, we get the molecule of each atoms, and if its water, we don't take it.
        if at1.get_parent().get_resname() in waternames \
            or at2.get_parent().get_resname() in waternames:
            continue
            

    if at1.get_serial_number() < at2.get_serial_number():      
        hblist.append([at1, at2]) 
    else:
        hblist.append([at2, at1])

print ()
print ("All parameters have been charged")


######### FROM HERE, WE WILL CALL TO EXTERNAL SCRIPT DEPENDING ON THE FUNCTION WE WANT.


while True:
    os.system("clear")

    mode = input("""
        Disponible functions:

        0) See the charged parameters
        1) Load a list with residue interactions
        2) Load a list with residue interactions and energies
        3) Make a plot of the type of residues interacting. 
        4) Ala Scanning (In progress)
        5) Exit
       
       Please, make your choice: """)

    if mode == '0':

        print()
        print ("Settings")
        print ("--------")
        for k,v in vars(args).items():   #K is the variable name, and V the value
            print ('{:10}:'.format(k),v) 
        print ("{} atom types loaded".format(vdwParams.ntypes))
        print ("{} amino acid atoms loaded".format(aaLib.nres))
        print ()
        enter = input("""
        Press any key to continue: """)    



    if mode == "1":
       
        simpleinteractions(hblist)
        enter = input("""
        Press any key to continue: """)


    if mode == '2':
       
        energies(hblist,aaLib,vdwParams,st, msms_path)
        enter = input("""
        Press any key to continue: """)


    if mode == '3':
       
        print (""" 
            Plot legend:

            If the interaction is only made in the Backbone == 1
            If the interaction is between side chains == 2
            If the interaction is between Backbone and Side chain == 3
    """)
        plot(hblist, backbone_polars)
        enter = input("""
        Press any key to continue: """)

    if mode == '5':
     
        print()
        print ("Thanks for using this program. Come again")
        print()
        break

    else:
        print()
        





